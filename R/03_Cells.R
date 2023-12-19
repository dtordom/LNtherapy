##########################
## Patterns behind Cellcept response
## R version 4.3.1 (2023-06-16)
##########################
## Get associations between cell proportions and response/non respose to LN drugs

library("ggplot2")
linrary("ggprubr")
set.seed(123456788)

load("DEGS.RData")

## ······································································· Step 1 
## Load data
