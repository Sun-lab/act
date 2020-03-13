# Permanova for TCGA survival time, 
# by comparing patients with long survival time sersus short survival time

rm(list=ls())
library(vegan)
library(e1071)
library(pROC)
source("../R_batch1/Unifrac_dist.R")
set.seed(2020)

## --------------------------------------------------------
## read in sample tcga
## --------------------------------------------------------
cf = "../data/COAD_composition.txt"
datC = read.table(cf, header = T, as.is = T, row.names = 1, sep='\t',
                  check.names=FALSE)
rownames(datC) = sub("^X", "", rownames(datC))
dim(datC)

clinic_dat = read.table("../data/expression_v2_sample.txt", sep="\t",
                        as.is=TRUE, header=TRUE)

clinic_dat$surv_time = ifelse(is.na(clinic_dat$diagnoses.days_to_death),
                              clinic_dat$diagnoses.days_to_last_follow_up,
                              clinic_dat$diagnoses.days_to_death)

clinic_dat$stage = gsub("stage ", "", clinic_dat$diagnoses.tumor_stage)
clinic_dat$stage = gsub('[a-c]', "", clinic_dat$stage)
clinic_dat$stage[which(clinic_dat$stage=="not reported")] = NA
clinic_dat$death = ifelse(clinic_dat$diagnoses.vital_status == 'dead', 1, 0)
clinic_dat$gender = clinic_dat$demographic.gender
clinic_dat$age_at_diagnosis = clinic_dat$diagnoses.age_at_diagnosis

sam = match(rownames(datC),clinic_dat$participant)
if(anyNA(sam)){
  stop('clinical data does not match with CIBERSOR output')
}

dA = cbind(datC, clinic_dat[sam, c("surv_time", "death",'stage',
                                   "diagnoses.age_at_diagnosis",
                                   'demographic.gender', 'participant',
                                   "diagnoses.days_to_death",
                                   "diagnoses.days_to_last_follow_up",
                                   "diagnoses.vital_status")])

dA = dA[!is.na(dA$surv_time) &!is.na(dA$death), ]
n = nrow(dA)
n

# long survival vs short 0 = short, 1 = 'long'
dA$surv_time = ifelse(is.na(dA$diagnoses.days_to_death), 
                              dA$diagnoses.days_to_last_follow_up, 
                              dA$diagnoses.days_to_death)
# 1 years 
dA$surv_1year = NA
dA$surv_1year[dA$diagnoses.days_to_death < 365.25] = 0
dA$surv_1year[dA$diagnoses.days_to_death > 365.25 | 
                      dA$diagnoses.days_to_last_follow_up > 365.25] = 1 
table(dA$surv_1year, useNA = 'ifany')
table(dA$surv_1year, 
      dA$surv_time >365.25, useNA = 'ifany')


# 2 years
dA$surv_2year = NA
dA$surv_2year[dA$diagnoses.days_to_death < 365.25*2] = 0
dA$surv_2year[dA$diagnoses.days_to_death > 365.25*2 | 
                      dA$diagnoses.days_to_last_follow_up > 365.25*2] = 1 
table(dA$surv_2year, useNA = 'ifany')
table(dA$surv_2year, 
      dA$surv_time >365.25*2, useNA = 'ifany')

# 5 years
dA$surv_5year = NA
dA$surv_5year[dA$diagnoses.days_to_death < 365.25*5] = 0
dA$surv_5year[dA$diagnoses.days_to_death > 365.25*5 | 
                      dA$diagnoses.days_to_last_follow_up > 365.25*5] = 1 
table(dA$surv_5year, useNA = 'ifany')
table(dA$surv_5year, 
      dA$surv_time >365.25*5, useNA = 'ifany')

table(dA$diagnoses.vital_status, 
      is.na(dA$diagnoses.days_to_death), useNA = 'ifany')


## --------------------------------------------------------
## distance matrices
## --------------------------------------------------------
datC = datC[dA$participant, ]

## 22 cell types
tree22 = read.tree(file="../data/immune_cell_lineage_tree_434_gene_22ct.tre")
tree22$edge.length
datP = NodeProp22(datC, tree22)  #for 22 cell types
D22ct = GUniFrac(datP, tree22, alpha = 0.5)
dim(D22ct)

## 12 cell types
tree12 = read.tree(file="/fh/fast/sun_w/licai/COAD/data_lineage_dist/immune_cell_lineage_tree_cor_434_gene.tre")
tree12$edge.length
datC2 = CellType12(datC)
dim(datC2)
datP  = NodeProp(datC2, tree12)  #for 12 cell types
D12ct = GUniFrac(datP, tree12, alpha = 0.5)
dim(D12ct)


table(colnames(D12ct) == dA$participant)
table(rownames(D22ct) == dA$participant)

## --------------------------------------------------------
## load dist & permanova for binary 
## --------------------------------------------------------

# source("/fh/fast/sun_w/licai/COAD/R_batch5/Permanova_cov_test.R")

# 1 year 
sam2kp1 = which(!is.na(dA$surv_1year))
length(sam2kp1)
pval12_year1 = adonis(D12ct[sam2kp1,sam2kp1] ~ as.factor(dA$surv_1year[sam2kp1]))
pval12_year1
pval22_year1 = adonis(D22ct[sam2kp1,sam2kp1] ~ as.factor(dA$surv_1year[sam2kp1]))
pval22_year1

# 2 year 
sam2kp1 = which(!is.na(dA$surv_2year))
length(sam2kp1)
pval12_year2 = adonis(D12ct[sam2kp1,sam2kp1] ~ as.factor(dA$surv_2year[sam2kp1]))
pval12_year2
pval22_year2 = adonis(D22ct[sam2kp1,sam2kp1] ~ as.factor(dA$surv_2year[sam2kp1]))
pval22_year2


# 5 year 
sam2kp1 = which(!is.na(dA$surv_5year))
length(sam2kp1)
pval12_year5 = adonis(D12ct[sam2kp1,sam2kp1] ~ as.factor(dA$surv_5year[sam2kp1]))
pval12_year5
pval22_year5 = adonis(D22ct[sam2kp1,sam2kp1] ~ as.factor(dA$surv_5year[sam2kp1]))
pval22_year5

## --------------------------------------------------------
## svm  
## --------------------------------------------------------
CTcomp = datC + 0.0001
CTcomp = log(CTcomp/CTcomp$`Macrophages M0`)
CTcomp = CTcomp[, -which(colnames(CTcomp) == "Macrophages M0")]

perm = 999
func1 <- function(surv, data){
  fitSVM1 = svm(as.factor(surv) ~ ., data = data)
  mean(fitSVM1$fitted == as.factor(surv))
}
 
# 1 year 
ind1  = which(is.na(dA$surv_1year))
dat1  = replicate(perm, sample(dA$surv_1year[-ind1]))
dat1  = cbind(dA$surv_1year[-ind1], dat1)
dim(dat1)

accuracy1 = apply(dat1, 2, func1, CTcomp[dA$participant[-ind1], ])
accuracy1[1]
mean(accuracy1[1] >= accuracy1[-1]) #pval

# 2 year 
ind2  = which(is.na(dA$surv_2year))
dat2  = replicate(perm, sample(dA$surv_2year[-ind2]))
dat2  = cbind(dA$surv_2year[-ind2], dat2)
dim(dat2)

accuracy2 = apply(dat2, 2, func1, CTcomp[dA$participant[-ind2], ])
accuracy2[1]
mean(accuracy2[1] >= accuracy2[-1]) #pval

# 5 year
ind5  = which(is.na(dA$surv_5year))
dat5  = replicate(perm, sample(dA$surv_5year[-ind5]))
dat5  = cbind(dA$surv_5year[-ind5], dat5)
dim(dat5)

accuracy5 = apply(dat5, 2, func1, CTcomp[dA$participant[-ind5], ])
accuracy5[1]
mean(accuracy5[1] >= accuracy5[-1]) #pval


q('no')