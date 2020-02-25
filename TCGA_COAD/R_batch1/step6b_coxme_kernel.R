
library(survival)
library(ggplot2)
library(ape)
library(geiger)
library(coxme)
source("Unifrac_dist.R")

kernel_tree <- function(distance, k, sigma){
  mu = apply(distance, 2, function(x) mean(sort(x)[1:k+1]))
  epsilon = outer(mu, mu, '+') * sigma/2
  ker = 1/(epsilon * sqrt(2*pi)) * exp(- distance^2/(2 * epsilon^2))
  return(list(kernels = ker, eigenvalue = eigen(ker)$value))
}

func1 <- function(nsim, n, dA, K1){ 
  set.seed(nsim + 2018)
  library(coxme)
  library(survival)
  samTrain = sample(n, size=round(n*0.7))
  samTest  = setdiff(1:n, samTrain)
  
  KernList = list(K1[samTrain, samTrain], diag(1, length(samTrain)))
  colnames(KernList[[1]]) = rownames(KernList[[1]]) = dA$participant[samTrain]
  colnames(KernList[[2]]) = rownames(KernList[[2]]) = dA$participant[samTrain]
  Mlist = coxmeMlist(KernList, rescale=F)
  coxME_fit = try(coxme(Surv(surv_time, death) ~ (1|participant), 
                        data=dA[samTrain,], varlist = Mlist))
  
  count = 0
  while(class(coxME_fit) == 'try-error' & count < 10){
    samTrain = sample(n, size=round(n*0.7))
    samTest  = setdiff(1:n, samTrain)
    
    KernList = list(K1[samTrain, samTrain], diag(1, length(samTrain)))
    colnames(KernList[[1]]) = rownames(KernList[[1]]) = dA$participant[samTrain]
    colnames(KernList[[2]]) = rownames(KernList[[2]]) = dA$participant[samTrain]
    Mlist = coxmeMlist(KernList, rescale=F)
    coxME_fit = try(coxme(Surv(surv_time, death) ~ (1|participant), 
                          data=dA[samTrain,], varlist = Mlist))
    count = count + 1
  }
  
  if(class(coxME_fit) == 'try-error'){
    return(NA)
    # return(list(Cs = NA, sigmas = c(NA, NA)))
  }
  # ---------------------------------------------------------------------
  # prediction on test set (test.inds)
  # ---------------------------------------------------------------------
  
  sigmas = c(coxME_fit$vcoef$participant)
  t_covr = sigmas[1]*K1 + sigmas[2]*diag(1,n)
  gtrain = coxME_fit$frail$participant
  gtrain = gtrain[match(dA$participant[samTrain],names(gtrain))]
  coef.mat  = solve(t_covr[samTrain, samTrain], t_covr[samTrain, samTest])
  yPredict3 = crossprod(coef.mat, gtrain)
  sfit3 = summary(coxph(Surv(surv_time, death) ~ yPredict3, data=dA[samTest,]))
  return(sfit3$concordance[1])
  #return(list(Cs = sfit3$concordance[1], sigmas = sigmas))
}

# ----------------------------------------------------------------------
# read in CIBERSORT output and clinic data 
# ----------------------------------------------------------------------
tree1 = read.tree(file="../data/immune_cell_lineage_tree_434_gene_22ct.tre")
tree1$edge.length

ks = 40
ks
sigmas = c(0.2, 0.25, 0.5)
sigmas
alphas = 0.5

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
                                   'demographic.gender', 'participant')])

dA = dA[!is.na(dA$surv_time) &!is.na(dA$death), ]
n = nrow(dA)
n

nsim = 100
Cs = array(NA, dim = c(nsim, length(ks) * length(sigmas), length(alphas) + 1),
           dimnames = list(c(1:nsim), seq(length(ks) * length(sigmas)), 
                           c(paste0('alpha',alphas), 'Euclidean'))) 

# L2 distance
datC = datC[dA$participant, ]
table(rownames(datC) == dA$participant )
apply(datC, 2, function(x) sum(x ==0))

distance  = as.matrix(dist(datC))
KList = list()
for(k in ks){
  for(sigma in sigmas){
    KList[[paste0(k, '_', sigma)]] = kernel_tree(distance, k, sigma)
  }
}
dimnames(Cs)[[2]] = names(KList)

kernel_all = lapply(KList, '[[', 'kernels')
ev = sapply(KList, '[[', 2)
print(apply(ev,2,function(x) sum(x <0) ))

for(i in 1:length(kernel_all)){
  cat(names(kernel_all)[i], date(), '\n')
  Cs[,i, 'Euclidean'] =  sapply(1:nsim, func1, n = n, dA = dA, K1 =  kernel_all[[i]])
}

# exp(-D)
K1 = exp(-distance)
expDl2 = sapply(1:nsim, func1, n = n, dA = dA, K1 =  K1)

# GUniFrac distance
datP = NodeProp22(datC, tree1)  #for 22 cell types
D0.5 = GUniFrac(datP, tree1, alpha = 0.5)
D0.5[1:5, 1:5]

KList2 = list()
for(k in ks){
  for(sigma in sigmas){
    KList2[[paste0(k, '_', sigma)]] = kernel_tree(D0.5, k, sigma)
  }
}

kernel_all2 = lapply(KList2, '[[', 'kernels')
ev = sapply(KList2, '[[', 2)
print(apply(ev,2,function(x) sum(x <0) ))

for(i in 1:length(kernel_all2)){
  cat(names(kernel_all2)[i], date(), '\n')
  Cs[,i, 'alpha0.5'] =  sapply(1:nsim, func1, n = n, dA = dA, K1 =  kernel_all2[[i]])
}

K1 = exp(-D0.5)
expD.5 = sapply(1:nsim, func1, n = n, dA = dA, K1 =  K1)

# correlation of 
Kcor = cor((t(datC) + 1)/2)
corrC = sapply(1:nsim, func1, n = n, dA = dA, K1 =  Kcor)

# save(Cs, expDl2,  file = '../data/cstats.Rdata' )
# load("../data/cstats.Rdat")

library(reshape2)
library(ggplot2)
Cs2 = melt(Cs)[,-1]
colnames(Cs2) = c('kernels', 'distance', 'C-statistics')
Cs2[1:4,]
levels(Cs2$distance)
levels(Cs2$distance) = c("GUniFrac(0.5)", "Euclidean")
levels(Cs2$kernels) <- paste0('GK(', sub("_", ',', levels(Cs2$kernels) ), ')')

expDl2 = cbind('kernels' = "exp(-D)", 'distance' = 'Euclidean', 
              'C-statistics' = expDl2)
expD.5 = cbind('kernels' = "exp(-D)", 'distance' = 'GUniFrac(0.5)', 
               'C-statistics' = expD.5)
corrC = cbind('kernels' = "correlation", 'distance' = 'correlation', 
               'C-statistics' = corrC)

Cs3 = rbind(Cs2, expD.5, expDl2, corrC)
Cs3$`C-statistics` = as.numeric(Cs3$`C-statistics`)


pdf(paste0('../figures/coxME_survival_COAD.pdf'))
ggplot(data = Cs3, aes(x = distance, y = `C-statistics`)) + 
  geom_boxplot(aes(fill=kernels)) +
  theme(legend.position ='top') + 
  scale_fill_discrete(name = 'kernels')
dev.off()

q('no')