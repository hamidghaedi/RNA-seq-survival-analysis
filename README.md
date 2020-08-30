# RNA seq survival analysis in R
Originally this was the method I used to do survival analysis on gene expression (RNA-seq) in bladder cancer TCGA data. For publishing here I decided to add more details and steps in a way that helps everybody who needs to get to know the basics and codes needed for cancer survival analysis on RNA-seq data.

## When is this needed?

In some cases when you have a list of differentially expressed genes/genes belongs to a specific GO term/somatically mutated genes it would be very powerful to find a correlation between the status of dysregulation in these genes and survival time of patients.

In this project, I performed  survival analysis on bladder cancer data (RNA-seq and DNA-seq) from TCGA. We were interested in finding whether: 

a) Do dysregulation of epigenetic-related genes is associated with bladder cancer patient survival? 

b) Do mutations in the epigenetic-related genes is associated with bladder cancer patient survival?  

Here we will focus only RNA-seq data. Before diving into analysis, I am going to share some basic needed to know to better understand the analysis steps.

## Survival Analysis intro
In fact, survival analysis corresponds to a number of statistical methods employed to find the time it takes for an event of interest to occur in a group of patients. To be more specific in cancer research survival we may use survival analysis to answer questions about “time” from operation(surgery) to patient death, “time” from starting a treatment regime to cancer progression, and “time” from response to a drug to disease recurrence. These are somehow classic use of survival analysis. However, here we will use “gene expression” and “mutation data” to assess their impact on patient survival. Indeed, we want to know correlation (if any) of specific gene dysregulation on bladder cancer patient survival “time”.

## Some  concepts:

Two basic concepts in doing survival analysis are survival time and the event in a study. Here the time from disease diagnosis to the occurrence of the event of interest (death) is referred to as survival time. In addition to “death”, there are other events in cancer studies: relapse and progression. One may consider “relapse” to do relapse-free survival analysis. Relapse is defined as the time between response to treatment and recurrence of the disease.

In the real-world scenarios in a cohort of patients, it is fairly common to have patients who we fail to follow them up, withdraw from the study, or may show no event in the defined time of the study. These situations result in observation which we called them censored observation. We will include censored observations in our analysis for the sake of having more and more data points, however, this needs to treat such data points as censored. For example, if a subject has no death data but has a date to last follow up, we can use this data until a certain time point (last follow update). In survival plot (Kaplan Meier plot) these data point is indicated with a + sign on survival line. To describe survival data, we will use survival probability which corresponds to the probability that a patient survives from the beginning of the time (for example diagnosis of cancer) to a particular future time (end of the study).


## Steps toward performing survival analysis

From here on, we will focus on gene expression data. For steps needed for survival analysis using mutation data please refer here. 

### Downloading  data
We can use two approaches to retrive data

#### Approach A:

Direct data retrieval, no need to install any packages. I have adapted this approach from a Biostar [post](https://www.biostars.org/p/153013/) by @Tris. 

1-	RNA data : Go to the FireBrowse ( http://gdac.broadinstitute.org/ ), select your dataset (we were interested on “Bladder urothelial carcinoma”) under “Data” column click "Browse". In the new webpage popup window scroll down to "mRNASeq" and then select ```illuminahiseq_rnaseqv2-RSEM_genes_normalized```. Download it, and extract the file ```BLCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt``` to your working directory. 

2-	Clinical data: In the popped window scroll down to the section “Clinical”, find ```Merge_Clinical``` and download it. Extract ```BLCA.merged_only_clinical_clin_format.txt``` into your working directory.

To read data in R:
```R
# loading librariies:
require(TCGABiolinks)
require(limma)
require(survival)
require(survminer)

###reading exptression matrix
rna <- read.table('BLCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt', sep = "\t", header = T, row.names = 1)
# looking at first rows:
head(rna) 
# as you can see, we have to remove the first row of the datset
# removing unwanted row:
rna <- rna[-1, ]
# row names should be gene name but as you can see, it is a combination of gene symbol and gene Entrez id. for example "A1BG" gene is indicated as "A1BG|1" . 
#We should polish row name to only contain gene symbol.
df <- data.frame(name = row.names(rna)) # keeping rownames as a temporary data frame
df <- data.frame(do.call('rbind', strsplit(as.character(df$name),'|',fixed=TRUE))) # this do magic like "text to column" in Excel!
df$X1[df$X1 == "?"] <- df$X2 # some genes are only presented by Entrez gene number, to keep these gene
rowName <- df$X1
# find duplicates in rowName, if any
table(duplicated(rowName))
#FALSE  TRUE 
#20530     1 
# in order to resilve duplucation issue
rowName[duplicated(rowName) == TRUE]
#[1] "SLC35E2"
grep("SLC35E2", rowName)
#[1] 16301 16302
rowName[16302] <- "SLC35E2_2"
#setting rna row names 
row.names(rna) <- rowName
rm(df, rowName) # removing datasets that we do not need anymore
###reading clinical data
clinical <- read.table('BLCA.merged_only_clinical_clin_format.txt',header=T, row.names=1, sep='\t', fill = TRUE) 
View(clinical)# it is better to transpose the data set
clinical <- t(clinical)
```

#### Approach B: 

Alternatively, it is possible to download data using third party package like ```TCGABiolink```:

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
I'd rather to use Approach B, since it is returning you most updated clinical data. However, both should work fine. 

### Data cleaning, recoding and transformation
However almost all genes included in the ```rna``` matrix, it is quite logical to have genes which show no expression (show 0 expression in all samples) or very low expression or uneven expression pattern ( having 0 value in >= 50% of cases). We need to keep these types of genes out from our analysis.

```R
#to find howmany genes show no expression in the rna matrix
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
In order to normalize RNA-seq data , we use ```voom``` function from ```limma``` package.  As stated in the manual this fuction "transform count data to 
log2-counts per million (logCPM), estimate the mean-variance relationship and use this to compute appropriate observational-level weights".
```voom``` needs to be supplied by a count matrix and a "design matrix with rows corresponding to samples and columns to coefficients to be estimated". 

Before diving into normalization step, we need to get familiarize ourselves with the TCGA barcode structure and meaning. Full description can be found [here](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/). 
A typical barcode is something like “TCGA-CF-A1HS-01A-11R-A13Y-07 “. As detailed by the TCGA working group letter 14 to 15 – here 01 denote sample type: 
Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. So the barcode in our example is a tumoral sample barcode. 
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

*Credits of function codes in this section goes to @Tris Biostar

To use gene expression matrix in survival analysis usually we encode genes as high or low expressed genes. To do so both fold change and z-score are fine.
However, due to retaining heterogeneity in data, the latter is preferred.
The general formula for calculating z-score is as 

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
```survival time``` could be find under columns that contains "days" in the clinical dataset.

```R
colnames(clinical)[grep("days", colnames(clinical))]
#[1] "days_to_collection"                           
#[2] "days_to_last_follow_up"                       
#[3] "days_to_diagnosis"                            
#[4] "days_to_birth"                                
#[5] "days_to_death"                                
#[6] "paper_Combined.days.to.last.followup.or.death"
```
We will consider ```"days_to_death"``` as survival time.Also in cases we have no ```"days_to_death"``` data, ```"days_to_last_follow_up"``` would be considered.
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
z_rna <- z_rna[, -"TCGA-K4-A4AB-01B-12R-A28M-07"]
```
Final steps before doing servival analysis is to encode RNA-seq data to dysregulated and intact. by dysregulated we mean genes with |z-score| > 1.96.
```R
dys_rna <- t(apply(z_rna, 1, function(x) ifelse(abs(x) > 1.96,"dysregulated","intact")))


```
### Performing survival analysis
We will use packages ```survival``` and ```survminer``` to do analysis. Suppose we are intrested in *EMP1* gene. It has been suggested that this gene is 
a survival gene for [Bladder cancer](https://www.nature.com/articles/s41420-020-00295-x). Further, in this paper *TPM1*, *NRP2*, *FGFR1*, *CAVIN1*, and *LATS2* were 
identified as bladder cancer survival-related genes. 
```R
fin_dat <- data.frame(gene = dys_rna[row.names(dys_rna) == "EMP1", ])
fin_dat <- merge( fin_dat, new_clin, by = 0)
#table(fin_dat$gene)
# fitting model
fit1 <- survfit(Surv(new_death, event) ~ gene, data = fin_dat)
print(fit1)

# calcilating pvalue
fit2 <- survdiff(Surv(new_death, event) ~ gene, data = fin_dat)
pv <- ifelse ( is.na(fit2),next,(round(1 - pchisq(fit2$chisq, length(fit2$n) - 1),3)))[[1]]
print(pv)
```
One of the most intresting aspect of survival analysis is to have survival probability in a graph (Kaplan–Meier curve). 
To draw KM curve:
```R
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
          pval = TRUE, conf.int = TRUE,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          surv.median.line = "hv", # Specify median survival
          ggtheme = theme_bw(),
          palette = c("#990000", "#000099"))
```
![alt text](https://github.com/hamid-gen/RNA_seq_survival_analysis_in_R/blob/master/surv.PNG)

##### performing survival analysis for all genes
To this aim we can use a for loop.
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
Inspect the result file carefully, p values should be interpreted in the context of having a balance sample size in both dysregulated and intact group. When one group - here dysregulated group is more likely has a low sample number, it is more likely to get you a significant p-value while this would not be true in most cases.

The following table represents a sub-set from the ```result``` table for six survival-related genes (mentioned above). As you can see, in our analysis some of these genes show significant association and two of them show non-significant p-value.

|gene    | pval    |dysregulated    |intact |
|------- |:-------:|:--------------:|------:|
|EMP1    | 0.005   | 117            |290    |
|FGFR1   | 0.051   | 273            |134    | 
|TPM1    | 0.056   | 54             |353    | 
|NRP2    | 0.147   | 88             |319    |
|LATS2   | 0.186   | 89             |318    |

