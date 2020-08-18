

library(survival)
library(limma)

# read RNA file 
rna <- read.table('RNA/BLCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt',nrows=20533, header=T,row.names=1,sep='\t')
# and take off first row cause we don't need it
rna <- rna[-1,]

#reading clinical data
clinical <- t(read.table('Clinical/BLCA.merged_only_clinical_clin_format.txt',header=T, row.names=1, sep='\t', fill = T))
clinical_keep_safe <- clinical



# first I remove genes whose expression is == 0 in more than 50% of the samples:
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}
remove <- rem(rna)
rna <- rna[-remove,]

#Now I need to identify normal and tumor samples. 
# see the values
table(substr(colnames(rna),14,14))

# get the index of the normal/control samples
n_index <- which(substr(colnames(rna),14,14) == '1')
t_index <- which(substr(colnames(rna),14,14) == '0')

# apply voom function from limma package to normalize the data
vm <- function(x){
  cond <- factor(ifelse(seq(1,dim(x)[2],1) %in% t_index, 1,  0))
  d <- model.matrix(~1+cond)
  x <- t(apply(x,1,as.numeric))
  ex <- voom(x,d,plot=F)
  return(ex$E)
}

rna_vm  <- vm(rna)
#colnames(rna_vm) <- gsub('\\.','-',substr(colnames(rna),1,12))
colnames(rna_vm) <- colnames(rna)
colnames(rna_vm) <- gsub("\\.", "-",colnames(rna_vm))
# and check how data look, they should look normally-ish distributed
hist(rna_vm)
# we can remove the old "rna" cause we don't need it anymor
rm(rna)

#Now we can finally scale the data. the reason to do so is because we don't want to 
#use ONLY the fold changes. if we use FC then we average the expression values across 
#all samples, losing the heterogeity that is characteristic of those data. we therefore 
#transform data to z-score so that per each patient for each gene we will 
#have a measure of how many SD away from the mean that is and we will consider those
#with Z > +/- 1.96 (roughly p=0.05 or 2 SD away) to be differentially expressed
#To obtain z-scores for the RNASeq data, we use following formula:
# zscore scaling against normals
scal <- function(x,y){
  mean_n <- rowMeans(y)  # mean of normal
  sd_n <- apply(y,1,sd)  # SD of normal
  # z score as (value - mean normal)/SD normal
  res <- matrix(nrow=nrow(x), ncol=ncol(x))
  colnames(res) <- colnames(x)
  rownames(res) <- rownames(x)
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      res[i,j] <- (x[i,j]-mean_n[i])/sd_n[i]
    }
  }
  return(res)
}
z_rna <- scal(rna_vm[,t_index],rna_vm[,n_index])



# set the rownames keeping only gene name
rownames(z_rna) <- sapply(rownames(z_rna), function(x) unlist(strsplit(x,'\\|'))[[1]]) 

rm(rna_vm)#we don't need it anymore

#
# match the patient ID in clinical data with the colnames of z_rna
#clinical <- data.frame(clinical)
#clinical$IDs <- toupper(clinical$patient.bcr_patient_barcode)
#clin_for_sur come from TCGA biolibks and ad.expression
sum(clin_for_surv$barcode %in% colnames(z_rna)) # we have 408 patients that we could use

# get the columns that contain data we can use: days to death, new tumor event, last day contact to....
#ind_keep <- grep('days_to_new_tumor_event_after_initial_treatment',colnames(clinical))
#
# this is a bit tedious, since there are numerous follow ups, let's collapse them together and keep the first value (the higher one) if more than one is available
#new_tum <- as.matrix(clinical[,ind_keep])
#new_tum_collapsed <- c()
#for (i in 1:dim(new_tum)[1]){
#  if ( sum ( is.na(new_tum[i,])) < dim(new_tum)[2]){
#    m <- min(new_tum[i,],na.rm=T)
#    new_tum_collapsed <- c(new_tum_collapsed,m)
#  } else {
#    new_tum_collapsed <- c(new_tum_collapsed,'NA')
#  }
#}

#
# do the same to death
ind_keep <- grep("paper_Days.until.death",colnames(clin_for_surv))
death <- as.matrix(clin_for_surv[,ind_keep])
death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if ( sum ( is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,'NA')
  }
}

# and days last follow up here we take the most recent which is the max number
ind_keep <- grep("paper_Days.to.last.followup",colnames(clin_for_surv))
fl <- as.matrix(clin_for_surv[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]){
  if ( sum (is.na(fl[i,])) < dim(fl)[2]){
    m <- max(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,'NA')
  }
}

# and put everything together
new_clin <- data.frame(death_collapsed,fl_collapsed)
colnames(new_clin) <- c('death_days', 'followUp_days')

#Since we want to do censored analysis, we need to have something to censor the data with.
#For example, if a patient has no death data BUT there is a date to last followup it means 
#that after that day we know nothing about the patient, therefore after that day it cannot 
#be used for calculations/Kaplan Meier plot anymore, therefore we censor it.
#so now we need to create vectors for both 'time to new tumor' and 'time to death' 
#that contain also the data from censored individuals.

# create vector with time to new tumor containing data to censor for new_tumor
new_clin$new_time <- c()
for (i in 1:length(as.numeric(as.character(new_clin$new_tumor_days)))){
  new_clin$new_time[i] <- ifelse ( is.na(as.numeric(as.character(new_clin$new_tumor_days))[i]),
                                   as.numeric(as.character(new_clin$followUp_days))[i],as.numeric(as.character(new_clin$new_tumor_days))[i])
}

# create vector time to death containing values to censor for death
new_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(new_clin$death_days)))){
  new_clin$new_death[i] <- ifelse ( is.na(as.numeric(as.character(new_clin$death_days))[i]),
                                    as.numeric(as.character(new_clin$followUp_days))[i],as.numeric(as.character(new_clin$death_days))[i])
}

# create vector for death censoring
table(clin_for_surv$paper_Vital.status)

#create death event
new_clin$death_event <- ifelse(clin_for_surv$paper_Vital.status == "Alive", 0,1)

#
#finally add row.names to clinical
rownames(new_clin) <- clin_for_surv$barcode

#Now, one more thing to do is to use the z-score values we obtained from 
#the RNASeq data and define which samples are altered or do not change in this 
#case I just look at mRNA deregulation, you can divide the data into up- and 
#down- regulated too if you wish.
# create event vector for RNASeq data
event_rna <- t(apply(z_rna, 1, function(x) ifelse(abs(x) > 1.96,1,0)))

# since we need the same number of patients in both clinical and RNASeq data take 
#the indices for the matching samples
ind_tum <- which(unique(colnames(z_rna)) %in% rownames(new_clin))
ind_clin <- which(rownames(new_clin) %in% colnames(z_rna))

#make row name more readable
df <- data.frame(rownames(z_rna))
rownames(z_rna) <- df$rownames.z_rna.

# pick your gene of interest
ind_gene <- which(rownames(z_rna) == "MAX")
# check how many altered samples we have
table(event_rna[ind_gene,])

##
# run survival analysis
s <- survfit(Surv(as.numeric(as.character(new_clin$new_death))[ind_clin],new_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum])
s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(new_clin$new_death))[ind_clin],new_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum]), error = function(e) return(NA))

# extraect the p.value
pv <- ifelse ( is.na(s1),next,(round(1 - pchisq(s1$chisq, length(s1$n) - 1),3)))[[1]]

print(pv)
# plot the data
plot(survfit(Surv(as.numeric(as.character(new_clin$new_death))[ind_clin],new_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum]),
     col=c(1:3), frame=F, lwd=2,main=paste('BLCA',rownames(z_rna)[ind_gene],sep='\n'), xlab = "days to death" , ylab = "survival probability")

# add lines for the median survival
# add lines for the median survival
x1 <- ifelse(is.na(as.numeric(summary(s)$table[,'median'][1])), "NA", 
             as.numeric(summary(s)$table[,'median'][1]))
x2 <- as.numeric(summary(s)$table[,'median'][2])
if(x1 != "NA" || x2 != "NA"){
  lines(c(0, x1), c(0.5, 0.5), col="blue")
  lines(c(x1, x1), c(0, 0.5), col="black")
  lines(c(x2, x2), c(0, 0.5), col="red")
}

# add legend
legend(4000, 1.150, 
       legend=paste("p.value = ", pv[[1]], sep=""), 
       bty="n", 
       cex=1.1)
legend(max(as.numeric(as.character(new_clin$death_days)[ind_clin]), na.rm = T) * 1.1, 1.0, 
       legend=c(paste("NotAltered=",x1), paste("Altered=",x2)), 
       bty="n",
       cex=1.1, 
       lwd=1, 
       col=c("black","red"))

