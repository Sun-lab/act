# genotype data requires large memory, please run this code on cluster
# args = commandArgs(TRUE)
args = c("22")
chr1 = args[1]

a = 0.5 # alpha

source("Unifrac_dist.R")
source("microbiomeGWAS-master/microbiomeGWAS.R")

packageDir = ('microbiomeGWAS-master')
system(paste0('cd ', packageDir, '; sh compile.src.sh'))

sampleFile = "/fh/fast/sun_w/licai/COAD/data_GWAS/tcga_gwas/sample2kp_dosageData.csv"
sam2kp = read.csv(sampleFile, header = F, as.is = T)$V1
nSam = length(sam2kp)
nSam

PCfile   = "/fh/fast/sun_w/research/TCGA/COAD/QC_PCA/final_data2_CaucasianOnly_COAD.txt"
geno_dir = "/fh/fast/sun_w/research/TCGA/COAD/GDC/data_geno/step12_get_dosage"
maf_dir  = "/fh/fast/sun_w/licai/COAD/data_GWAS/tcga_gwas/maf/"
# ------------------------------------------------------------
# CaucasianOnly PCs
# ------------------------------------------------------------
tcga_pca2 = read.table(PCfile, sep = "\t", header = TRUE, as.is = TRUE)
dim(tcga_pca2)
tcga_pca2[1:5, ]
barcode = sapply(strsplit(tcga_pca2$TCGA_ID, "-" ), 
                 function(x) paste(x[1:3], collapse = '-')) 
barcode[1:5]

coad = read.table("../data/expression_v2_sample.txt", sep="\t",
                  as.is=TRUE, header=TRUE)
dim(coad)
coad[1:2,]
table(coad$demographic.race)

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
# 22 cell types
tree <-   read.tree(file = "../data/immune_cell_lineage_tree_434_gene_22ct.tre")

CTcomp = read.table("../data/COAD_composition.txt",
                        header=T,sep="\t",row.names = 1, check.names=F)
rownames(CTcomp) = gsub('^X', '', rownames(CTcomp))

# distance
datP = NodeProp22(CTcomp, tree)
dG_E <- GUniFrac(datP, tree, alpha = a)

#correlation
dG_E <- dG_E[rownames(adjustment_data), rownames(adjustment_data)]

# ------------------------------------------------------------
# adjust for covariates
# ------------------------------------------------------------
cov = c(paste0('noHM_PC', 1:4), "genderM", "diagnoses.age_at_diagnosis")

dataCovariate = adjustment_data[, cov]
distRes_E = distRes(dG_E, adjustment_data[, cov])

eD_E <- funED(dG_E, soFile1 = paste0(packageDir, "/lib/dExp1.so"), 
              soFile2 = paste0(packageDir, "/lib/dExp2.so"))


# load dosage data
cat('chromosome', chr1, date(), '\n')

# maf 005 6 million snp
Maf   = read.csv(paste0(maf_dir, "/maf_chr", chr1 , '.txt'), sep='\t' )
index = which(Maf$MAF >= 0.05 & Maf$MAF <= 1 - 0.05)

# genotype file and snpinfo
file  = paste0(geno_dir,"/dosage/chr",chr1,".rds")
dosage = readRDS(file)

dosage_barcode = sapply(strsplit(colnames(dosage), "-" ), 
                        function(x) paste(x[1:3], collapse = '-')) 
snp_info = dosage[index, 1:5]
dosage   = dosage[index, match(sam2kp, dosage_barcode)]
dosage   = data.matrix(dosage)
mafs     = apply(dosage, 1, sum)/(2 *ncol(dosage)) 

tasks  = 40
groups = cut(1:nrow(dosage), breaks =tasks, include.lowest = T, labels = F)
dos    = lapply(1:tasks, function(x) dosage[which(groups ==x),])
rm(dosage)
gc()

SM_E = NULL
for( d in 1:tasks){
  SM2  = t(apply(dos[[d]], 1,sm, distRes_E))
  SM_E = c(SM_E, SM2)
}
mafs = round(mafs, 15)
eG    = funEG(mafs)
p_E   = funMain(SM_E, dG_E, eD_E, eG)
res_E = cbind(snp_info, p_E)

write.table(res_E, file = paste0('../data/microbiome_tcga_22ct_dG',a,
                                 '_E_GWAS_cov_chr',chr1,'.txt'), quote = FALSE,
            sep = "\t", row.names = FALSE, col.names = TRUE)

q(save = 'no')


