
library(survival)
## --------------------------------------------------------
## read in sample clinical data
## --------------------------------------------------------

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


## --------------------------------------------------------
## read in immune cell composition file
## --------------------------------------------------------

CTcomp = read.table("../data/COAD_composition.txt", header=T, 
                    sep="\t", row.names = 1, check.names=F)
rownames(CTcomp) = sub('^X','',rownames(CTcomp))
CTcomp[1:5,]

apply(CTcomp, 2, function(x) sum(x == 0))

CTcomp = CTcomp + 0.0001
CTcomp = log(CTcomp/CTcomp$`Macrophages M0`)
CTcomp = CTcomp[, -which(colnames(CTcomp) == "Macrophages M0")]

table(rownames(CTcomp) %in% clinic_dat$participant)
clinic_dat = clinic_dat[match(rownames(CTcomp),clinic_dat$participant),]
table(clinic_dat$participant == rownames(CTcomp))

table(is.na(clinic_dat$surv_time))
table(is.na(clinic_dat$death))
table(is.na(clinic_dat$stage))
table(is.na(clinic_dat$gender))
table(is.na(clinic_dat$age_at_diagnosis))

## --------------------------------------------------------
## Cox regression: age, gender
## --------------------------------------------------------
# CTcomp_qn = t(apply(CTcomp, 1, normscore)) #  quantile normalization

cox_surv_ct = matrix(NA, nrow = ncol(CTcomp), ncol = 7)
rownames(cox_surv_ct) = colnames(CTcomp)


for(i in colnames(CTcomp)){
  fit1 = coxph(Surv(surv_time, death) ~ CTcomp[,i]  + gender + 
                age_at_diagnosis, data=clinic_dat)
  print(i)
  cox_surv_ct[i,] = c(summary(fit1)$n,summary(fit1)$nevent,
                      summary(fit1)$coefficient[1,])
  print(summary(fit1))
}


colnames(cox_surv_ct)  = c('n', 'nevent', colnames(summary(fit1)$coefficient))

cox_surv_ct
signif(cox_surv_ct, 3)
write.table(signif(cox_surv_ct, 3), file = "../data/cox_surv.txt")



## --------------------------------------------------------
## Cox regression: age, gender, tumor stage
## --------------------------------------------------------

cox_surv_ct2 = matrix(NA, nrow = ncol(CTcomp), ncol = 7)
rownames(cox_surv_ct2) = colnames(CTcomp)


for(i in colnames(CTcomp)){
  fit1 = coxph(Surv(surv_time, death) ~ CTcomp[,i] + gender + stage +
                 age_at_diagnosis, data=clinic_dat)
  print(i)
  cox_surv_ct2[i,] = c(summary(fit1)$n,summary(fit1)$nevent,
                      summary(fit1)$coefficient[1,])
  print(summary(fit1))
}


colnames(cox_surv_ct2)  = c('n', 'nevent', 
                            colnames(summary(fit1)$coefficient))

cox_surv_ct2
signif(cox_surv_ct2, 3)
write.table(signif(cox_surv_ct2, 3), file = "../data/cox_surv_with_stage.txt")

## --------------------------------------------------------
## Cox regression: dichotomize
## --------------------------------------------------------
cox_surv_ct = matrix(NA, nrow = ncol(CTcomp), ncol = 8)
rownames(cox_surv_ct) = colnames(CTcomp)


for(i in colnames(CTcomp)){
  ctc = cut(CTcomp[,i], breaks = c(-0.01, median(CTcomp[,i]),1),
            labels =c("Low", "High") ) 
  fit1 = coxph(Surv(surv_time, death) ~ ctc + gender + 
                 age_at_diagnosis, data=clinic_dat)
  print(i)
  cox_surv_ct[i,] = c(summary(fit1)$n,
                      summary(fit1)$nevent,
                      sum(ctc == "Low"),
                      summary(fit1)$coefficient[1,])
  print(summary(fit1))
}


colnames(cox_surv_ct)  = c('n', 'nevent',  "NLow",
                           colnames(summary(fit1)$coefficient))

cox_surv_ct

## --------------------------------------------------------
## Cox regression: dichotomize
## --------------------------------------------------------
cox_surv_ct2 = matrix(NA, nrow = ncol(CTcomp), ncol = 8)
rownames(cox_surv_ct2) = colnames(CTcomp)


for(i in colnames(CTcomp)){
  ctc = cut(CTcomp[,i], breaks = c(-0.01, median(CTcomp[,i]),1),
            labels =c("Low", "High") ) 
  fit1 = coxph(Surv(surv_time, death) ~ ctc + gender + stage + 
                 age_at_diagnosis, data=clinic_dat)
  print(i)
  cox_surv_ct2[i,] = c(summary(fit1)$n,
                      summary(fit1)$nevent,
                       sum(ctc == "Low"),
                      summary(fit1)$coefficient[1,])
  print(summary(fit1))
}


colnames(cox_surv_ct2)  = c('n', 'nevent',  "NLow",
                           colnames(summary(fit1)$coefficient))

cox_surv_ct2


## --------------------------------------------------------
## C-index for macrophage M1
## --------------------------------------------------------

func1 <- function(nsim, ctD, clinic_dat){
  set.seed(2018 + nsim)
  n = length(ctD)
  samTrain = sample(n, size=round(n*0.7))
  samTest  = setdiff(1:n, samTrain)
  fit1 = coxph(Surv(surv_time, death) ~ ctD[samTrain], 
               data=clinic_dat[samTrain,])
  return(summary(fit1)$concordance[1])
}

Cindex = sapply(1:100, func1, CTcomp[, "Macrophages M1"], clinic_dat)
save(Cindex, file ="../data/M1_cstat.Rdata")
pdf("../figures/C-index for log ratio Macrophages M1.pdf")
boxplot(Cindex, xlab = "Macrophages M1", ylab = "C-statistics")
dev.off()

q("no")