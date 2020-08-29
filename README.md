# RNA seq survival analysis in R
Originally this was the method I used to do survival analysis on gene expression (RNA-seq) in bladder cancer TCGA data. 
For publishing here I decided to add more details and steps in a way that helps everybody who need to get to know basics and codes needed for survival.

## When is this needed?

In some cases when you have a list of differentially expressed genes/genes belongs to a specific GO term/somatically mutated genes it would be very powerful to find correlation between status of dysregulation in these genes and survival time of patients. 

In this project I performed survival analysis on bladder cancer data (RNA-seq and DNA-seq) from TCGA. We were interested in to find whether :
a) Do dysregulation of epigenetic-related genes is associated with bladder cancer patient survival?
b) Do mutations in the epigenetic-related genes is associated with bladder cancer patient survival?
Before diving into analysis, I am going to share some basic needed to know to better understand the analysis steps. 

## Introduction to Survival Analysis

In fact, survival analysis corresponds to a number of statistical methods employed to find  the time it takes for an event of interest to occur in a group of patients. 
To be more specific in cancer research survival we may use survival analysis to answer question about “time” from operation(surgery) to patient death, “time” from starting a treatment regime to cancer progression and “time” from response to a drug to disease recurrence. These are somehow classic use of survival analysis. However, here we will use “gene expression” and “mutation data” to assess their impact on patient survival. Indeed, we want to know correlation (if any) of specific gene dysregulation on bladder cancer patient survival “time”. 

## Some  concepts:

Two basic concepts in doing survival analysis are survival time and event in a study. Here the time from disease diagnosis to the occurrence of the event of interest (death) is referred as survival time. In addition to “death”, there are other events in cancer studies: relapse and progression. One may consider “relapse” to do relapse-free survival analysis. Relapse is defined as the time between response to treatment and recurrence of the disease. 

In the real-world scenarios in a cohort of patients it is fairly common to have patients who we fail to follow them up, withdraw from the study or may show no event in the defined time of the study. These situations result in observation which we called them censored observation. We will include censored observations in our analysis for the sake of having more and more data points, however this needs to treat such data point as censored. For example,  if a subject have no death data but have a date to last follow up, we can use this data until to the certain time point (last follow up date). In survival plot (Kaplan Meier plot) these data point is indicated with a + sign on survival line. 
In order to describe survival data, we will use survival probability which corresponds to the probability that a patient survives from the beginning of the time (for example diagnosis of cancer) to a particular future time (end of the study). 


## Steps toward performing survival analysis

From here on, we will focus on gene expression data. For steps needed for survival analysis using mutation data please refer here. 

### Downloading  data
We can use two approaches to retrive data

#### Approach A:

Direct data retrieval, no need to install any packages. I have adapted this approach from a Biostar [post](https://www.biostars.org/p/153013/) by @Tris. Also some code chuncks for data preparation and analysis are adopted from this post. 

1-	RNA data : Go to the FireBrowse ( http://gdac.broadinstitute.org/ ), select your dataset (we were interested on “Bladder urothelial carcinoma”) under “Data” column click "Browse". In the new webpage popup window scroll down to "mRNASeq" and then select ```illuminahiseq_rnaseqv2-RSEM_genes_normalized```. Download it, and extract the file ```BLCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt``` to your working directory. 

2-	Clinical data: In the popped window scroll down to the section “Clinical”, find ```Merge_Clinical``` and download it. Extract ```BLCA.merged_only_clinical_clin_format.txt``` into your working directory.

To read data in R:
```R
# loading librariies:
require(TCGABiolink)
require(limma)

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
# make column name shorter to only have sample name (1,12 letter)
colnames(rna_vm) <- substr(colnames(rna_vm),1,12)
# Since no longer "rna" dataset is needed, we can remove it. 
rm(rna)

# After these steps we expect that data to have a somehow gaussian distribution.
hist(rna_vm)
```
![alt text](https://github.com/hamid-gen/RNA_seq_survival_analysis_in_R/blob/master/rna_vm.png)

#### RNA-seq data scaling and encoding

To use gene expression matrix in survival analysis usually we encode genes as high or low expressed genes. To do so both fold change and z-score are fine.
However due to retaining heterogeneity in data the latter is preferred. Using z-score we will have a measure of how many SD away from the mean a gene is. Also we will consider those with |Z| > 1.96 to be differentially expressed: + 1.96 (up-regulated) and -1.96 (down-regulated)

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
z_rna <- scal(rna_vm[,t_index],rna_vm[,n_index])

rm(rna_vm)
```
