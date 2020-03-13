rm(list = ls())
setwd("/fh/fast/sun_w/licai/_Cancer/data_immu")
source("/fh/fast/sun_w/licai/COAD/R_batch5/Unifrac_dist.R")
library(ggplot2)
library(ape)
library(geiger)
library(vegan)

L2 = T # Euclidean distance or correlation
# distance
kernel = F
k = NULL
sigma = NULL #0.2 #0.25 #0.5 #


outdir = paste0("/fh/fast/sun_w/licai/COAD/data_lineage_dist/tcga/maf005/revision/ksigma")
system(sprintf("mkdir -p %s", outdir))

source("/fh/fast/sun_w/licai/COAD/R_batch5/microbiomeGWAS-master/microbiomeGWAS.R")

packageDir = ('/fh/fast/sun_w/licai/COAD/R_batch5/microbiomeGWAS-master')
setwd(packageDir)
system(paste0('cd ', packageDir, '; sh compile.src.sh'))

sam2kp = read.csv( "/fh/fast/sun_w/licai/COAD/data_GWAS/tcga_gwas/sample2kp_dosageData.csv", header = F, as.is = T)$V1
nSam = length(sam2kp)

# ------------------------------------------------------------
#CaucasianOnly 
# ------------------------------------------------------------
tcga_pca2 = read.table("/fh/fast/sun_w/research/TCGA/COAD/QC_PCA/final_data2_CaucasianOnly_COAD.txt",
                       sep = "\t", header = TRUE, as.is = TRUE)
dim(tcga_pca2)
tcga_pca2[1:5, ]

coad = read.table("/fh/fast/sun_w/research/TCGA/COAD/GDC/data/expression_v2_sample.txt", sep="\t",
                  as.is=TRUE, header=TRUE)
dim(coad)
coad[1:2,]
table(coad$demographic.race)

barcode = sapply(strsplit(tcga_pca2$TCGA_ID, "-" ), function(x) paste(x[1:3], collapse = '-')) 

coad = coad[match(sam2kp, coad$bcr_patient_barcode), ]
tcga_pca2 = tcga_pca2[match(sam2kp, barcode), ]
coad[1:5,47]
tcga_pca2[1:5, 1]

# ------------------------------------------------------------
# covariates 
# ------------------------------------------------------------
adjustment_data = cbind(tcga_pca2[,c(paste0("noHM_PC", 1:4))], 
                        coad[, c('diagnoses.age_at_diagnosis','bcr_patient_barcode')], 
                        genderM = coad[,"demographic.gender"] == "male")
table(adjustment_data$genderM)
table(coad[,'demographic.gender'])

rownames(adjustment_data) = coad[,"participant"]
adjustment_data[1:5,]
dim(adjustment_data)

# ------------------------------------------------------------
# Distance / corelation
# ------------------------------------------------------------
tree1 <-   read.tree(file = "/fh/fast/sun_w/licai/COAD/data_lineage_dist/immune_cell_lineage_tree_434_gene_22ct.tre")
tree1$edge.length

datC = read.table(paste0("/fh/fast/sun_w/licai/_Cancer/data_immu/", 'COAD_composition.txt'),
                  header = T, as.is = T, row.names = 1, sep='\t',
                  check.names=FALSE)
rownames(datC) = sub("^X", "", rownames(datC))
# datP = NodeProp22(datC, tree1)  #for 22 cell types
if(L2){
	dG_E = as.matrix(dist(datC))  #GUniFrac(datP, tree1, alpha = a)
}else{
	dG_E = 1 - (cor(t(datC)) + 1)/2 
}
dim(dG_E)

if(kernel){
	mu = apply(dG_E, 2, function(x) mean(sort(x)[1:k+1]))
	epsilon = outer(mu, mu, '+') * sigma/2
	dG_E = dG_E^2/(2 * epsilon^2)
}
# # renormalize
# DE = diag(1/sqrt(rowSums(dG_E)))
# dG_E = DE %*% dG_E %*% DE
# colnames(dG_E) = rownames(dG_E) = colnames(datE)

#correlation
dG_E <- dG_E[rownames(adjustment_data), rownames(adjustment_data)]
dim(dG_E)
dG_E[1:5, 1:5]

# ------------------------------------------------------------
# adjust for covariates
# ------------------------------------------------------------
cov = c(paste0('noHM_PC', 1:4), "genderM", "diagnoses.age_at_diagnosis")

dataCovariate = adjustment_data[, cov]
distRes_E = distRes(dG_E, adjustment_data[, cov])

eD_E <- funED(dG_E, soFile1 = paste0(packageDir, "/lib/dExp1.so"), 
              soFile2 = paste0(packageDir, "/lib/dExp2.so"))

library(Rmpi)
mpi.spawn.Rslaves(nslaves=mpi.universe.size() - 1)
time.1 = proc.time()

for(chr1 in 5:1){
  # load dosage data
  cat('chromosome', chr1, date(), '\n')
  
  # maf 005 6 million snp
  Maf = read.csv(paste0("/fh/fast/sun_w/licai/COAD/data_GWAS/tcga_gwas/maf/maf_chr", 
                        chr1 , '.txt'), sep='\t' )
  index = which(Maf$MAF >= 0.05 & Maf$MAF <= 1 - 0.05)
  
  
  # genotype file and snpinfo
  # genotype file
  file = paste0("/fh/fast/sun_w/research/TCGA/COAD/GDC/data_geno/step12_get_dosage/dosage/chr",chr1,".rds")
  dosage = readRDS(file)
  
  dosage_barcode <- sapply(strsplit(colnames(dosage), "-" ), function(x) paste(x[1:3], collapse = '-')) 
  snp_info <- dosage[index, 1:5]
  dosage <- dosage[index, match(sam2kp, dosage_barcode)]
  #rownames(dosage) <- snp_info$rsid
  dosage = data.matrix(dosage)
  mafs = apply(dosage, 1, sum)/(2 *ncol(dosage)) 
  
  tasks = 40
  groups = cut(1:nrow(dosage), breaks =tasks, include.lowest = T, labels = F)
  dos = lapply(1:tasks, function(x) dosage[which(groups ==x),])
  sapply(dos, dim)
  rm(dosage)
  gc()
  
  SM_E = NULL
  for( d in 1:tasks){
    SM2 = t(mpi.parApply(dos[[d]], 1,sm, distRes_E))
    SM_E = c(SM_E, SM2)
  }
  
  
  eG = funEG(mafs)
  
  p_E = funMain(SM_E, dG_E, eD_E, eG)
  
  res_E = cbind(snp_info, p_E)
  
  write.table(res_E, 
              file = paste0(outdir,'/microbiome_L2_k',k,'sigma',
                             sigma,'_GWAS_cov_chr',
                            chr1,'.txt'), quote = FALSE,
              sep = "\t", row.names = FALSE, col.names = TRUE)
}

time.2 = proc.time()
mpi.close.Rslaves()
time.2 - time.1

mpi.exit()

q("no")

# ----------------------------------------------------------------------------
# combine resutls
# ----------------------------------------------------------------------------

# distance
k = 40
sigma = 0.2 #0.25 #0.5 #

outdir = paste0("/fh/fast/sun_w/licai/COAD/data_lineage_dist/tcga/maf005/revision/k"
                ,k,"sigma",sigma)

resE_all = resC_all = NULL
for(chr1 in 1:22){
  print(chr1)
  
  res_E = read.table(paste0(outdir,'/microbiome_Gunifrac_k',k,'sigma',
                            sigma,'_GWAS_cov_chr',
                            chr1,'.txt'),
                     as.is = T, header = T)
  resE_all = rbind(res_E, resE_all)
}

write.table(resE_all, 
            file = paste0('microbiome_Gunifrac_k',k,'sigma',
                          sigma,'_GWAS_all.txt'), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


table(resE_all$pMK < 5e-8)

