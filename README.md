# RNA_seq_survival_analysis_in_R
Survival analysis on gene expression (RNA-seq) in bladder cancer TCGA data

# When is this needed?

In some cases when you have a list of differentially expressed genes/genes belongs to a specific GO term/somatically mutated genes it would be very powerful to find correlation between status of dysregulation in these genes and survival time of patients. 

In this project I performed survival analysis on bladder cancer data (RNA-seq and DNA-seq) from TCGA. We were interested in to find whether :
a) Do dysregulation of epigenetic-related genes is associated with bladder cancer patient survival?
b) Do mutations in the epigenetic-related genes is associated with bladder cancer patient survival?
Before diving into analysis, I am going to share some basic needed to know to better understand the analysis steps. 

# Introduction to Survival Analysis

In fact, survival analysis corresponds to a number of statistical methods employed to find  the time it takes for an event of interest to occur in a group of patients. 
To be more specific in cancer research survival we may use survival analysis to answer question about “time” from operation(surgery) to patient death, “time” from starting a treatment regime to cancer progression and “time” from response to a drug to disease recurrence. These are somehow classic use of survival analysis. However, here we will use “gene expression” and “mutation data” to assess their impact on patient survival. Indeed, we want to know correlation (if any) of specific gene dysregulation on bladder cancer patient survival “time”. 

# Some  concepts:

Two basic concepts in doing survival analysis are survival time and event in a study. Here the time from disease diagnosis to the occurrence of the event of interest (death) is referred as survival time. In addition to “death”, there are other events in cancer studies: relapse and progression. One may consider “relapse” to do relapse-free survival analysis. Relapse is defined as the time between response to treatment and recurrence of the disease. 

In the real-world scenarios in a cohort of patients it is fairly common to have patients who we fail to follow them up, withdraw from the study or may show no event in the defined time of the study. These situations result in observation which we called them censored observation. We will include censored observations in our analysis for the sake of having more and more data points, however this needs to treat such data point as censored. For example,  if a subject have no death data but have a date to last follow up, we can use this data until to the certain time point (last follow up date). In survival plot (Kaplan Meier plot) these data point is indicated with a + sign on survival line. 
In order to describe survival data, we will use survival probability which corresponds to the probability that a patient survives from the beginning of the time (for example diagnosis of cancer) to a particular future time (end of the study). 


# Steps toward performing survival analysis

From here on, we will focus on gene expression data. For steps needed for survival analysis using mutation data please refer here. 

1-	Downloading intended data

2-	Data wrangling

3-	Performing survival analysis
