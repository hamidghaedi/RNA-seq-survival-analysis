# RNA_seq_survival_analysis_in_R
Survival analysis on gene expression (RNA-seq) in bladder cancer TCGA data

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

Direct data retrieval, no need to install any packages.

1-	RNA data : Go to the FireBrowse ( http://gdac.broadinstitute.org/ ), select your dataset (we were interested on “Bladder urothelial carcinoma”) under “Data” column click "Browse". In the new webpage popup window scroll down to "mRNASeq" and then select ```illuminahiseq_rnaseqv2-RSEM_genes_normalized```. Download it, and extract the file ```BLCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt``` to your working directory. 

2-	Clinical data: In the popped window scroll down to the section “Clinical”, find ```Merge_Clinical``` and download it. Extract ```BLCA.merged_only_clinical_clin_format.txt``` into your working directory.

To read data in R:
```R
###reading exptression matrix
rna <- read.table('BLCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt', sep = "\t", header = T, row.names = 1)
# looking at first rows:
head(rna) 
# as you can see, we have to remove the first row of the datset
# removing unwanted row:
rna <- rna[-1, ]
# row names should be gene name but as you can see, it is a combination of gene symbol and gene Entrez id. for example _"A1BG"_ gene is indicated as _"A1BG|1"_ . 
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

###reading clinical data
clinical <- read.table('BLCA.merged_only_clinical_clin_format.txt',header=T, row.names=1, sep='\t', fill = TRUE) 
View(clinical)# it is better to transpose the data set
clinical <- t(clinical)
```

#### Approach B: 

Alternatively, it is possible to download data using third party package ```TCGABiolink```:

1-	RNA data
```R
library(TCGAbiolink)
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
```
clinical <- data.frame(dat@colData)
```
### Downloading  data
