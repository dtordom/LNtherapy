##########################
## Get DEG genes beetwen responders/ non-responders
## R version 4.3.1 (2023-06-16)
##########################

## Load Gene expression data and metadata from:
## https://drive.google.com/drive/folders/15KVn3tckVPZKieCsSZizdi8bAQFZkOIC?usp=sharing

load("Dataset.RData")

DEGS<-list()
DATA<-list()

## MMF
