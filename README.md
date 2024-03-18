# RNA seq survival analysis in R
Originally this was the method I used to do survival analysis on gene expression (RNA-seq) in bladder cancer TCGA data. For publishing here I decided to add more details and steps in a way that helps everybody who needs to get to know the basics and codes needed for cancer survival analysis on RNA-seq data.

## When is this needed?

In some cases when you have a list of differentially expressed genes/genes belonging to a specific GO term/somatically mutated genes it would be very powerful to find a correlation between the status of dysregulation in these genes and the survival time of patients.

In this project, I performed  survival analysis on bladder cancer data (RNA-seq and DNA-seq) from TCGA. We were interested in finding out whether: 

a) Do dysregulation of epigenetic-related genes is associated with bladder cancer patient survival? 

b) Do mutations in the epigenetic-related genes is associated with bladder cancer patient survival?  

Here we will focus only RNA-seq data. Before diving into the analysis, I am going to share some basics needed to know to better understand the analysis steps.

## Survival Analysis intro
In fact, survival analysis corresponds to a number of statistical methods employed to find the time it takes for an event of interest to occur in a group of patients. To be more specific in cancer research survival we may use survival analysis to answer questions about “time” from operation(surgery) to patient death, “time” from starting a treatment regime to cancer progression, and “time” from response to a drug to disease recurrence. These are somehow classic uses of survival analysis. However, here we will use “gene expression” and “mutation data” to assess their impact on patient survival. Indeed, we want to know the correlation (if any) of specific gene dysregulation on bladder cancer patient survival “time”.

## Some  concepts:

Two basic concepts in doing survival analysis are survival time and the event in a study. Here the time from disease diagnosis to the occurrence of the event of interest (death) is referred to as survival time. In addition to “death”, there are other events in cancer studies: relapse and progression. One may consider “relapse” to do a relapse-free survival analysis. Relapse is defined as the time between response to treatment and recurrence of the disease.

In real-world scenarios in a cohort of patients, it is fairly common to have patients who we fail to follow up, withdraw from the study or may show no event in the defined time of the study. These situations result in observation which we call censored observation. We will include censored observations in our analysis for the sake of having more and more data points, however, this needs to treat such data points as censored. For example, if a subject has no death data but has a date to last follow-up, we can use this data until a certain time point (last follow-up update). In the survival plot (Kaplan Meier plot) these data point is indicated with a + sign on the survival line. To describe survival data, we will use survival probability which corresponds to the probability that a patient survives from the beginning of the time (for example diagnosis of cancer) to a particular future time (end of the study).


## Steps toward performing survival analysis
There are three general steps toward doing a survival analysis: 1) Downloading data, 2)Data cleaning, recoding, and transformation, and 3)Performing survival analysis.

### 1-Downloading  data
We can use two approaches to retrieve data

#### Approach A:

Direct data retrieval, no need to install any packages. I have adapted this approach from a Biostar [post](https://www.biostars.org/p/153013/) by @Tris. 

1-	RNA data: Go to the FireBrowse ( http://gdac.broadinstitute.org/ ), select your dataset (we were interested in “Bladder urothelial carcinoma”) under the “Data” column click "Browse". In the new webpage popup window scroll down to "mRNASeq" and then select ```illuminahiseq_rnaseqv2-RSEM_genes_normalized```. Download it, and extract the file ```BLCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt``` to your working directory. 

2-	Clinical data: In the popped window scroll down to the section “Clinical”, find ```Merge_Clinical``` and download it. Extract ```BLCA.merged_only_clinical_clin_format.txt``` into your working directory.

To read data in R:
```R
# loading libraries:
require(TCGAbiolinks)
require(SummarizedExperiment)
require(limma)
require(survival)
require(survminer)

###reading expression matrix
rna <- read.table('BLCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt', sep = "\t", header = T, row.names = 1)
# looking at first rows:
head(rna) 
# as you can see, we have to remove the first row of the dataset
# removing unwanted row:
rna <- rna[-1, ]
# row names should be gene names but as you can see, it is a combination of gene symbol and gene Entrez id. For example "A1BG" gene is indicated as "A1BG|1" . 
#We should polish the row name to only contain the gene symbols.
df <- data.frame(name = row.names(rna)) # keeping rownames as a temporary data frame
df <- data.frame(do.call('rbind', strsplit(as.character(df$name),'|',fixed=TRUE))) # This do magic like "text to column" in Excel!
df$X1[df$X1 == "?"] <- df$X2 # some genes are only presented by Entrez gene number, to keep these genes
rowName <- df$X1
# find duplicates in rowName, if any
table(duplicated(rowName))
#FALSE  TRUE 
#20530     1 
# in order to resolve the duplication issue
rowName[duplicated(rowName) == TRUE]
#[1] "SLC35E2"
grep("SLC35E2", rowName)
#[1] 16301 16302
rowName[16302] <- "SLC35E2_2"
# Setting RNA row names 
row.names(rna) <- rowName
rm(df, rowName) # removing datasets that we do not need anymore
###reading clinical data
clinical <- read.table('BLCA.merged_only_clinical_clin_format.txt',header=T, row.names=1, sep='\t', fill = TRUE) 
View(clinical)# it is better to transpose the data set
clinical <- t(clinical)
```

#### Approach B: 

Alternatively, it is possible to download data using third-party package like ```TCGABiolink```:

1-	RNA data
```R

query <- GDCquery(project = "TCGA-BLCA",
                           data.category = "Gene expression",
                           data.type = "Gene expression quantification",
                           platform = "Illumina HiSeq", 
                           file.type  = "normalized_results",
                           experimental.strategy = "RNA-Seq",
                           legacy = TRUE)
GDCdownload(query, method = "api")
dat <- GDCprepare(query = query, save = TRUE, save.filename = "exp.rda")
rna <- as.data.frame(SummarizedExperiment::assay(dat))
```

2-	Clinical data: following the step one, 
```R
clinical <- data.frame(dat@colData)
```
I'd rather to use Approach B, since it is returning you the most updated clinical data. However, both should work fine. 

### 2-Data cleaning, recoding, and transformation
However almost all genes included in the ```rna``` matrix, it is quite logical to have genes that show no expression (show 0 expression in all samples) or very low expression or uneven expression pattern ( having 0 value in >= 50% of cases). We need to keep these types of genes out from our analysis.

```R
#to find how many genes show no expression in the rna matrix
table(rowSums(rna) == 0)
#FALSE  TRUE 
#19677   270
#visualizing RNA read counts distribution
hist(log10(rowSums(rna)), main = "log10-RNA read count dist")
```

![alt text](https://github.com/hamid-gen/RNA_seq_survival_analysis_in_R/blob/master/log10-RNA%20read%20count.png "log10-RNA read count")

```R
#to remove genes with no-expression in more 50% of samples:
fif <- dim(rna)[2]/2
no_exp <- data.frame(count = apply(rna, 1, function(x) length(which(x== 0))))
no_exp$gene <- row.names(no_exp)
no_exp <- no_exp[no_exp$count > fif, ]$gene
rna <- rna[- which(row.names(rna) %in% no_exp), ]
```

#### RNA-seq data normalization
In order to normalize RNA-seq data, we use ```voom``` function from ```limma``` package.  As stated in the manual this fuction "transform count data to 
log2-counts per million (logCPM), estimate the mean-variance relationship and use this to compute appropriate observational-level weights".
```voom``` needs to be supplied by a count matrix and a "design matrix with rows corresponding to samples and columns to coefficients to be estimated". 

Before diving into the normalization step, we need to familiarize ourselves with the TCGA barcode structure and meaning. A full description can be found [here](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/). 
A typical barcode is something like “TCGA-CF-A1HS-01A-11R-A13Y-07 “. As detailed by the TCGA working group letters 14 to 15 – here 01 denotes sample type: 
Tumor types range from 01 - 09, normal types from 10 - 19, and control samples from 20 - 29. So the barcode in our example is a tumoral sample barcode. 
To identify how many tumor and normal samples we have in our data we can do: 

```R
table(substr(colnames(rna),14,15))
 #01  11  # so we have 408 tumors and 19 normal samples.
#408  19 
```
Now using the barcode we can make indexes for tumor and normal samples

```R
normal_index <- which(substr(colnames(rna),14,14) == '1')
tumor_index <- which(substr(colnames(rna),14,14) == '0')

# apply voom function from limma package to normalize the data
vm <- function(x){
  cond <- factor(ifelse(seq(1,dim(x)[2],1) %in% tumor_index, 1,  0))
  d <- model.matrix(~1+cond)
  x <- t(apply(x,1,as.numeric))
  ex <- voom(x,d,plot=F)
  return(ex$E)
}

rna_vm  <- vm(rna)

# restoring column names
colnames(rna_vm) <- colnames(rna)
# After these steps we expect that data to have a somehow gaussian distribution.
hist(rna_vm)
```
![alt text](https://github.com/hamid-gen/RNA_seq_survival_analysis_in_R/blob/master/rna_vm.png)

#### RNA-seq data scaling and encoding

*Credits of function code in this section go to @Tris Biostar

To use gene expression matrix in survival analysis usually we encode genes as high or low expressed genes. To do so both fold change and z-score are fine.
However, due to retaining heterogeneity in data, the latter is preferred.
The general formula for calculating the z-score is as 

```z = [(gene X expression value in tumor)-(mean gene X expression value in normal)]/(standard deviation gene X expression in normal)```.

Using z-score we will have a measure of how many SD away from the mean a gene is. Also we will consider those genes  with |Z| > 1.96 to be differentially expressed: + 1.96 (up-regulated) and -1.96 (down-regulated)

```R
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
z_rna <- scal(rna_vm[,tumor_index],rna_vm[,normal_index])

rm(rna_vm)
```
Now it is time to define ```survival time``` and ```event``` in ```clinical``` dataset. 
```survival time``` could be found under columns that contain "days" in the clinical dataset.

```R
colnames(clinical)[grep("days", colnames(clinical))]
#[1] "days_to_collection"                           
#[2] "days_to_last_follow_up"                       
#[3] "days_to_diagnosis"                            
#[4] "days_to_birth"                                
#[5] "days_to_death"                                
#[6] "paper_Combined.days.to.last.followup.or.death"
```
We will consider ```"days_to_death"``` as survival time. Also in cases we have no ```"days_to_death"``` data, ```"days_to_last_follow_up"``` would be considered.
```R
clinical$new_death <- ifelse(is.na(clinical$days_to_death), clinical$days_to_last_follow_up, clinical$days_to_death)
clinical$new_death[clinical$new_death == 0] <- NA
```
Data for ```event``` could be found under column ```vital_status``` . We may want to recode this column.

```R
# to see what we have in th vital_status column
table(clinical$vital_status)

#       Alive         Dead Not Reported 
#         235          191            1 
# exclude patient with not reported vital_status
clinical[clinical$vital_status == "Not Reported", ]$barcode
#[1] "TCGA-K4-A4AB-01B-12R-A28M-07"

clinical <- clinical[-which(row.names(clinical) == "TCGA-K4-A4AB-01B-12R-A28M-07"), ]
#recoding vital_status
clinical$event <- ifelse(clinical$vital_status == "Alive", 0,1)
# create a subset from original clinical data
new_clin <- clinical[, c("new_death", "event")]


# remove sample with "not reported" vital status from expression matrix
z_rna <- z_rna[, - grep("TCGA-K4-A4AB-01B", colnames(z_rna))]
```
The final step before doing survival analysis is to encode RNA-seq data to be dysregulated and intact. by dysregulated we mean genes with |z-score| > 1.96.
```R
dys_rna <- t(apply(z_rna, 1, function(x) ifelse(abs(x) > 1.96,"dysregulated","intact")))

```
#### Data encoding by maximally selected rank statistics

As an alternative, we can use maximally selected rank statistics to identify samples with high (low) expression values for a given gene of interest. In fact, maximally selected rank statistics look for the optimal cutpoint that maximizes the separation between samples based on their expression levels of the gene. By ranking the samples and testing different cutpoints, this statistical method helps determine the most significant division of samples into groups with distinct expression patterns. This approach allows researchers to identify samples with high (low) expression values, enabling further investigation into the biological implications of these expression levels and their association with specific conditions or phenotypes.

There is an R package that allows for using maximally selected rank statistics implementation, `maxtstat`, see [here](https://cran.r-project.org/web/packages/maxstat/vignettes/maxstat.pdf) for more information and also there is a function `surv_cutpoint` from `survminer` package, ([link]((https://www.rdocumentation.org/packages/survminer/versions/0.4.9/topics/surv_cutpoint))) which uses this method to cut samples into low and high expression group while considering the survival data. 


### 3-Performing survival analysis
We will use packages ```survival``` and ```survminer``` to do analysis. Suppose we are interested in the *EMP1* gene. It has been suggested that this gene is 
a survival gene for [Bladder cancer](https://www.nature.com/articles/s41420-020-00295-x). Further, in this paper *TPM1*, *NRP2*, *FGFR1*, *CAVIN1*, and *LATS2* were 
identified as bladder cancer survival-related genes. 
```R
fin_dat <- data.frame(gene = dys_rna[row.names(dys_rna) == "EMP1", ])
fin_dat <- merge( fin_dat, new_clin, by = 0)
#table(fin_dat$gene)
# fitting model
fit1 <- survfit(Surv(new_death, event) ~ gene, data = fin_dat)
print(fit1)

# calculating pvalue
fit2 <- survdiff(Surv(new_death, event) ~ gene, data = fin_dat)
pv <- ifelse ( is.na(fit2),next,(round(1 - pchisq(fit2$chisq, length(fit2$n) - 1),3)))[[1]]
print(pv)
```
One of the most intresting aspect of survival analysis is to have survival probability in a graph (Kaplan–Meier curve). 
To draw KM curve:
```R
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit1,
          pval = TRUE, conf.int = TRUE,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          surv.median.line = "hv", # Specify median survival
          ggtheme = theme_bw(),
          palette = c("#990000", "#000099"))
```
![alt text](https://github.com/hamid-gen/RNA_seq_survival_analysis_in_R/blob/master/surv.PNG)

#### using maxstat for survival analysis:
Below is the code that can be used to perform survival analysis using the `maxstat` method. 
**Note** I did not run the code, so the following should be considered as a guide rather than a working code chunk and it's adapted from [this post](http://www.sthda.com/english/wiki/survminer-0-2-4):

```r
# 1. Determine the optimal cutpoint of EMP1
res.cut <- surv_cutpoint(myeloma, time = "time", event = "event",
   variables = c("DEPDC1", "WHSC1", "CRIM1"))
summary(res.cut)
##        cutpoint statistic
## DEPDC1    279.8  4.275452
## WHSC1    3205.6  3.361330
## CRIM1      82.3  1.968317
# 2. Plot cutpoint for DEPDC1
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(res.cut, "DEPDC1", palette = "npg")
## $DEPDC1
survminer
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
##           time event DEPDC1 WHSC1 CRIM1
## GSM50986 69.24     0   high   low  high
## GSM50988 66.43     0    low   low   low
## GSM50989 66.50     0    low  high  high
## GSM50990 42.67     1    low  high   low
## GSM50991 65.00     0    low   low   low
## GSM50992 65.20     0   high   low   low
# 4. Fit survival curves and visualize
fit <- survfit(Surv(time, event) ~DEPDC1, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE)
```



### Performing survival analysis for all genes
To this aim, we can use a for loop.
```R
all_gene <- row.names(dys_rna)
result = data.frame( gene=character(0), pval=numeric(0), dysregulated=numeric(0), intact=numeric(0))

for (i in all_gene){
    fin_dat <- data.frame(gene = dys_rna[row.names(dys_rna) == i, ])
    fin_dat <- merge( fin_dat, new_clin, by = 0)
    if (dim(table(fin_dat$gene)) > 1){
    fit2 <- survdiff(Surv(new_death, event) ~ gene, data = fin_dat)
    pv <- ifelse ( is.na(fit2),next,(round(1 - pchisq(fit2$chisq, length(fit2$n) - 1),3)))[[1]]

   gene <- i
   dysregulated <- table(fin_dat$gene)[1]
   intact <- table(fin_dat$gene)[2]
   pval = pv
   result[i, ] = c(gene, pval, dysregulated, intact)
    }
}

```
Inspect the result file carefully, p values should be interpreted in the context of having a balanced sample size in both dysregulated and intact groups. When one group - here dysregulated group is more likely to have a low sample number, it is more likely to get you a significant p-value while this would not be true in most cases.

The following table represents a sub-set from the ```result``` table for six survival-related genes (mentioned above). As you can see, in our analysis some of these genes show significant association and two of them show non-significant p-value.

|gene    | pval    |dysregulated    |intact |
|------- |:-------:|:--------------:|------:|
|EMP1    | 0.005   | 117            |290    |
|FGFR1   | 0.051   | 273            |134    | 
|TPM1    | 0.056   | 54             |353    | 
|NRP2    | 0.147   | 88             |319    |
|LATS2   | 0.186   | 89             |318    |

### Finding all survival-related gene 

In some cases for example for the ARID1A gene, we can see three distinct classes of expression: low (z-score <= -1.96), normal (-1.96 < z-score > 1.96), and high (z-score >= 1.96) In contrast to simply have dysregulated and intact expression, it is possible to have genes with a high, normal and low expression for survival analysis. To do so please consider this code:


```R
################ Finding genes which have values under the three categories High, Norm, Low
gene_table <- data.frame( gene=character(0), High=numeric(0), Norm = numeric(0), Low=numeric(0), lab = character(0))
gene <- row.names(dys_rna)
for (i in gene) {
  df <- data.frame(gene = dys_rna[row.names(dys_rna) == i, ])
  tb <- table(df)
  if (dim(tb) == 3){
    gene <- i
    High <- tb[1]
    Low <- tb[2]
    Norm <- tb[3]
    lab <- paste(names(tb)[1],names(tb)[2],names(tb)[3], sep = "_" )
    gene_table[i, ] = c(gene, High, Norm, Low, lab)
  }
  
}

############## finding genes which have data only in two of the three states (High,Norm,Low)
gene_table_2x2 <- data.frame( gene=character(0), state1=numeric(0), state2= numeric(0), lab1 = character(0), lab2 = character(0))

gene <- row.names(dys_rna)
for (i in gene) {
  df <- data.frame(gene = dys_rna[row.names(dys_rna) == i, ])
  tb <- table(df)
  if (dim(tb) == 2){
    gene <- i
    state1 <- tb[1]
    state2 <- tb[2]
    lab1 <- paste(names(tb)[1])
    lab2 <- paste(names(tb)[2])
    gene_table_2x2[i, ] = c(gene, state1, state2, lab1, lab2)
  }
  
}

#h.table
h.tab <- gene_table_2x2[gene_table_2x2$lab1 == "High", ]
names(h.tab)[2] <- names(table(h.tab$lab1))[1]
names(h.tab)[3] <- names(table(h.tab$lab2))[1]
h.tab$Low <- NA
h.tab <- h.tab[, colnames(gene_table)[1:4]]

#l.table
l.tab <- gene_table_2x2[gene_table_2x2$lab1 == "Low", ]
names(l.tab)[2] <- names(table(l.tab$lab1))[1]
names(l.tab)[3] <- names(table(l.tab$lab2))[1]
l.tab$High <- NA
l.tab <- l.tab[, colnames(gene_table)[1:4]]
#
gene_table_2x2 <- rbind(h.tab, l.tab)
gene_table_2x2$lab <- NA
#
gene_table <- rbind(gene_table, gene_table_2x2)
gene_table <- gene_table[,-5]

# geting the result table from the previous analysis:
gene_table[,2:4] <- sapply(gene_table[, 2:4], as.numeric) # setting type of columns for numbers as numeric
gene_table[is.na(gene_table)] <- 0

# defining classes for each gene
gene_table$state <- ifelse(gene_table$Norm < 15 & gene_table$High >= 15 & gene_table$Low >= 15, "HL", 
                        ifelse(gene_table$Low < 15 & gene_table$High >= 15 & gene_table$Norm >= 15, "HN", 
                               ifelse(gene_table$High < 15 & gene_table$Low >= 15 & gene_table$Norm >= 15, "LN",
                                      ifelse(gene_table$High >= 15 & gene_table$Low >= 15 & gene_table$Norm >= 15, "HNL", "flag"))))
# see what we get
table(gene_table$state)

#hnl
hnl.dys_rna <- t(apply(z_rna[row.names(z_rna) %in% gene_table[gene_table$state == "HNL", ]$gene, ], 1, function(x) ifelse(x >= 1.96,"High",ifelse(x <= -1.96, "Low", "Norm"))))

gene <- row.names(hnl.dys_rna)
hnl.result = data.frame( gene=character(0), pval=numeric(0), High=numeric(0), Norm=numeric(0), Low=numeric(0))

for (i in gene){
  fin_dat <- data.frame(gene = dys_rna[row.names(dys_rna) == i, ])
  fin_dat <- merge( fin_dat, new_clin, by = 0)
  if (dim(table(fin_dat$gene)) > 1){
    fit2 <- survdiff(Surv(new_death, event) ~ gene, data = fin_dat)
    pv <- ifelse ( is.na(fit2),next,(round(1 - pchisq(fit2$chisq, length(fit2$n) - 1),3)))[[1]]
    
    gene <- i
    High <- table(fin_dat$gene)[1]
    Norm <- table(fin_dat$gene)[3]
    Low <- table(fin_dat$gene)[2]
    pval = pv
    hnl.result[i, ] = c(gene, pval, High,Norm, Low)
  }
}

#hn
hn.dys_rna <- t(apply(z_rna[row.names(z_rna) %in% gene_table[gene_table$state == "HN", ]$gene, ], 1, function(x) ifelse(x >= 1.96,"High",ifelse(x <= -1.96, "Low", "Norm"))))
gene <- gene_table[gene_table$newcol == "HN", ]$gene
hn.result = data.frame( gene=character(0), pval=numeric(0), High=numeric(0), Norm=numeric(0))

for (i in all_gene){
  fin_dat <- data.frame(gene = hn.dys_rna[row.names(hn.dys_rna) == i, ])
  fin_dat <- merge( fin_dat, new_clin, by = 0)
  fin_dat <- fin_dat[-which(fin_dat$gene == "Low"), ]
  if (dim(table(fin_dat$gene)) > 1){
  fit2 <- survdiff(Surv(new_death, event) ~ gene, data = fin_dat)
  pv <- ifelse ( is.na(fit2),next,(round(1 - pchisq(fit2$chisq, length(fit2$n) - 1),3)))[[1]]
  #
  gene <- i
  High <- table(fin_dat$gene)[1]
  Norm <- table(fin_dat$gene)[2]
  pval = pv
  hn.result[i, ] = c(gene, pval, High, Norm)
  }
}

#ln
ln.dys_rna <- t(apply(z_rna[row.names(z_rna) %in% gene_table[gene_table$state == "LN", ]$gene, ], 1, function(x) ifelse(x >= 1.96,"High",ifelse(x <= -1.96, "Low", "Norm"))))

gene <- gene_table[gene_table$newcol == "LN", ]$gene
ln.result = data.frame( gene=character(0), pval=numeric(0), Low=numeric(0), Norm=numeric(0))

for (i in all_gene){
  fin_dat <- data.frame(gene = ln.dys_rna[row.names(ln.dys_rna) == i, ])
  fin_dat <- merge( fin_dat, new_clin, by = 0)
  fin_dat <- fin_dat[-which(fin_dat$gene == "High"), ]
  if (dim(table(fin_dat$gene)) > 1){
  fit2 <- survdiff(Surv(new_death, event) ~ gene, data = fin_dat)
  pv <- ifelse ( is.na(fit2),next,(round(1 - pchisq(fit2$chisq, length(fit2$n) - 1),3)))[[1]]
  #
  gene <- i
  Low <- table(fin_dat$gene)[1]
  Norm <- table(fin_dat$gene)[2]
  pval = pv
  ln.result[i, ] = c(gene, pval, High, Norm)
  }
}

# Combining results into one dataset
hn.result$Low <- NA
ln.result$High <- NA
hn.result <- hn.result[, colnames(hnl.result)]
ln.result <- ln.result[, colnames(hnl.result)]
#
result <- rbind(hnl.result, hn.result, ln.result)
# Selecting significantly associated genes
result$pval <- as.numeric(result$pval)
sig.result <- result[result$pval <= 0.05 ,]

write.table(sig.result, file = "sig.surv.associated.gene.csv", row.names = T, quote = F)
```


