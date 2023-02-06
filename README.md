# Windows of Susceptibility Analysis for Zika Virus

We just include the required R-package installation, introduction, conclusions, and references parts in the `README.md`
A complete report is avaliable at [Notion](https://tianrui-qi.notion.site/WOS-Analysis-for-Zika-Virus-d4bfad0386634b29b4c815b8a26db7f9), same as the HTML version knitted from `report.Rmd` through [RStudio](https://www.rstudio.com). 

## Required R-package Installation

```r
r = getOption("repos")
r["CRAN"] = "http://cran.rstudio.com"
options(repos = r)

if (!require("gplots")) {
  install.packages("gplots")
  library(gplots)
}
if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}
if (!require("ggbiplot")) {
   devtools::install_git("https://github.com/vqv/ggbiplot.git")
   library(ggbiplot)
}
if (!require("ggbiplot")) {
   install.packages("ggbiplot")
   library(ggbiplot)
}
if (!require("ComplexHeatmap")) {
  library(devtools)
  install_github("jokergoo/ComplexHeatmap")
  library(ComplexHeatmap)
}
if (!require("fields")) {
  install.packages("fields")
  library(fields)
}
if (!require("hablar")) {
  install.packages("hablar")
  library(hablar)
}
if (!require("knitr")) {
  install.packages("knitr")
  library(knitr)
}
if (!require("matlab")) {
   install.packages("matlab")
   library(matlab)
}
if (!require("tibble")) {
  install.packages("tibble")
  library(tibble)
}
if (!require("tidyr")) {
  install.packages("tidyr")
  library(tidyr)
}
knitr::opts_chunk$set(echo = TRUE)
```

## Introduction

From [WHO](https://www.who.int/news-room/fact-sheets/detail/zika-virus), Zika virus disease is caused by a virus transmitted primarily by Aedes mosquitoes, which bite during the day. 
Symptoms are generally mild and include fever, rash, conjunctivitis, muscle and joint pain, malaise or headache.
Symptoms typically last for 2–7 days. Most people with Zika virus infection do not develop symptoms. 
Zika virus infection during pregnancy can cause infants to be born with Microcephaly and other congenital malformations, known as congenital Zika syndrome. 
Infection with Zika virus is also associated with other complications of pregnancy including preterm birth and miscarriage.
An increased risk of neurologic complications is associated with Zika virus infection in adults and children, including Guillain-Barré syndrome, neuropathy and myelitis.

Researchers led by Dr. Bennett used existing data from in-vitro models of human brain development based on embryonic stem cells to understand when embryos are particularly susceptible to Zika infection resulting in Microcephaly.
These models produce the layers of the cerebral cortex in a petri dish.
The study is uses data from the [Cortecon](http://cortecon.neuralsci.org/) that consists of RNA transcription data measured as the layers of the cerebral cortex develop at nine different time points between 0 and 77 days.
An analysis of this data was published in a 2019 paper in the journal [Nature](https://www.nature.com/articles/s41598-019-39318-8).
The prior study found that in humans that genes associated with microcephaly are enriched for the Neuroectoderm Stage (Cluster 2).
Genes associated with microcephaly and changed by Zika infection have an even stronger enrichment. 
This suggests that human embryos are particularly susceptible to Zika-induced microcephaly in the first trimester, possibly before the mother even know she is pregnant.

Clustering was used to identify the stages of development in humans.
Principal component analysis (PCA), data visualization, and log odds ratio analysis were used to
show when these stages occur.
An advanced visualization, called a Susceptibility Window Ontological Transcription (SWOT) Clock, is available at [here](https://semnext.tw.rpi.edu/swotclock/).
Given a set of genes associated with a disease, the SWOT clock can be used to help identify time periods of cortical development that may be susceptible to changes in those genes.
Recall each cluster of gene represents a distinct stage of development, if a set of genes occurs more than expected (a.k.a enriched) for a given stages then it is likely that brain development may be impacted by changes in the activities of those genes during that stage.
The SWOT Clock provides a visualization and a statistical analysis of the log odds ratio of each cluster along with its p-value.
Since the set of disease genes are significantly enriched in a cluster, then the time period associated with that cluster is a window of susceptibility. 
A cluster is significantly enriched if it has a positive log odds ratio and a p-value < 0.1.

Perform the windows of susceptibility (WOS) analysis based on mouse data from a similar brain-in-a-dish model for mice instead of human data used by Dr. Bennett. 
Analyze the same sets of microcephaly-associated genes and Zika-associated genes to see if we can detect a similar WOS for Microcephaly and Zika-induced microcephaly in mice as in humans. 
Same technique also applied to cognitive disorders, anther disease that impacts mental function, to see if we get a similar result. 

## Conclusions

We got similar WOS result in human and mouse model for both Zika and cognitive disorder.
For Zika virus, analysis in both models show the second stage is most susceptible. 
The last stage, the upper layers stage, is most susceptible for the cognitive disorders in human model where the stages D and E are most susceptible in mouse model, which is also end stage. 

## References

This is a mini project for the MATP 4400 Data Introduction to Data Mathematics by Dr. Kristin Bennett and Dr. John Erickson, Spring 2022, at Rensselaer Polytechnic Institute, Troy, NY. The background and goal of the project is given in *IDM Lab 5: Mini Project*. 
The structure of this report is revise on the basis of *PreLab 5* and *Lab 5* templates given by instructors. 
Dr. Bennett and Dr. Erickson gave perfect instruction about all the methods including K-meaning clustering, PCA, and log odds ratio analysis and R visualization techniques in MATP 4400. 
