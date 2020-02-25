args=(commandArgs(TRUE))
chr1 = as.numeric(args[1])
chr1 

date <- Sys.Date()

sampleFile = "/fh/fast/sun_w/licai/COAD/data_GWAS/tcga_gwas/sample2kp_dosageData.csv"
sam2kp = read.csv(sampleFile, header = F, as.is = T)$V1
nSam = length(sam2kp)
nSam

PCfile   = "/fh/fast/sun_w/research/TCGA/COAD/QC_PCA/final_data2_CaucasianOnly_COAD.txt"
geno_dir = "/fh/fast/sun_w/research/TCGA/COAD/GDC/data_geno/step12_get_dosage"
maf_dir  = "/fh/fast/sun_w/licai/COAD/data_GWAS/tcga_gwas/maf/"

# ------------------------------------------------------------
#CaucasianOnly 
# ------------------------------------------------------------
tcga_pca2 = read.table(PCfile,
                       sep = "\t", header = TRUE, as.is = TRUE)
dim(tcga_pca2)
tcga_pca2[1:5, ]

coad = read.table("../data/expression_v2_sample.txt", sep="\t",
                  as.is=TRUE, header=TRUE)
dim(coad)
coad[1:2,]
table(coad$demographic.race)

tcga_pca2$TCGA_ID = sapply(strsplit(tcga_pca2$TCGA_ID, "-" ), 
                           function(x) paste(x[1:3], collapse = '-')) 
coad = coad[match(sam2kp, coad$bcr_patient_barcode), ]
tcga_pca2 = tcga_pca2[match(sam2kp, tcga_pca2$TCGA_ID), ]
table(coad$bcr_patient_barcode ==tcga_pca2$TCGA_ID )

# ------------------------------------------------------------
# covariates 
# ------------------------------------------------------------
adjustment_data = cbind(tcga_pca2[,c("TCGA_ID", paste0("noHM_PC", 1:4))], 
                        coad[, c("demographic.gender", 'diagnoses.age_at_diagnosis','bcr_patient_barcode')])
rownames(adjustment_data) = coad[,"participant"]
adjustment_data[1:5,]
dim(adjustment_data)


# ------------------------------------------------------------
# immune cell composition
# ------------------------------------------------------------

CTcomp = read.table("../data/COAD_composition.txt", header=T, 
                    sep="\t", row.names = 1, check.names=F)
rownames(CTcomp) = sub('^X','',rownames(CTcomp))
CTcomp[1:5,]
CTcomp = CTcomp[rownames(adjustment_data), ]
apply(CTcomp, 2, function(x) sum(x != 0))

CTcomp = CTcomp + 0.0001
CTcomp = log(CTcomp/CTcomp$`Macrophages M0`)
CTcomp = CTcomp[, -which(colnames(CTcomp) == "Macrophages M0")]

table(rownames(adjustment_data) == rownames(CTcomp))

# ------------------------------------------------------------
# immune cell composition
# ------------------------------------------------------------

DoCiber <- function(g, cell_type, adjustment_data, datCiber){
  dat = cbind(dosage = as.numeric(g), adjustment_data)
  composition = datCiber[,cell_type ]
  if(length(unique(g)) > 1 ){
    fit1 <- try(lm(composition ~ . , data = dat ))
    res = c(summary(fit1)$coefficients[2,])
  }
  else{ 
    res = c(rep(NA, 4))
  }
  names(res) <- c('beta', 'sd', 't_value', 'p_value')
  return(res)
}


# chr1 = 22
cat('chromosome', chr1, date(), '\n')

# genotype file
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

table(rownames(CTcomp) == sapply(strsplit(colnames(dosage), "-" ), 
                                 function(x) x[3] ))
table(rownames(adjustment_data) == rownames(CTcomp))

covs = c(paste0('noHM_PC', 1:4), "demographic.gender", "diagnoses.age_at_diagnosis")

# tasks = 20 
# groups = cut(1:nrow(dosage), breaks =tasks, include.lowest = T, labels = F)
# dos = lapply(1:tasks, function(x) dosage[which(groups ==x), ])
# 
# rm(dosage)
# gc()

for(cell_type in colnames(CTcomp)){
  cat(cell_type, date(), ' \n')
  Ciber = NULL
  # for( d in 1:tasks){
  #   Ciber_m = t(apply(dos[[d]], 1, DoCiber, cell_type, 
  #                            adjustment_data = adjustment_data[, covs], 
  #                            datCiber = CTcomp))
  #   Ciber = rbind(Ciber, Ciber_m)
  # }
  #test
  Ciber = t(apply(dosage, 1,DoCiber, cell_type, adjustment_data = adjustment_data[, covs], datCiber = CTcomp))
  Ciber = cbind(Ciber, snp_info)
  write.table(Ciber, 
              file = paste0('/fh/fast/sun_w/licai/COAD/data_GWAS/tcga_gwas/log_ratio/', 
                            "tcga_gwas_maf0.05_chr",chr1, '_', cell_type, '.txt'), quote = FALSE, 
              sep = "\t", row.names = F, col.names = TRUE)
  
  
}



q(save = 'no')

for(chri in 22:1){
  com = sprintf("sbatch R CMD BATCH '--args %s' step7a_cell_type_specific_gwas.R step7a_cell_type_specific_gwas_chr%s.Rout",
                chri, chri)
  message(com)
}


# ------------------------------------------------------------
# collect results
# ------------------------------------------------------------

res_folder = '/fh/fast/sun_w/licai/COAD/data_GWAS/tcga_gwas/log_ratio/'

CTcomp = read.table("../data/COAD_composition.txt", header=T, 
                    sep="\t", row.names = 1, check.names=F)
rownames(CTcomp) = sub('^X','',rownames(CTcomp))
CTcomp[1:5,]
CTcomp = CTcomp[,-which(colnames(CTcomp) == "Macrophages M0")]

for(cell_type in colnames(CTcomp)[14:21]){
  cat(cell_type, date(), '\n')
  resA = NULL
  for(chr1 in 1:22){
  cat(chr1, date(), '\n')
    
  resi = read.table(paste0(res_folder, 
                          "tcga_gwas_maf0.05_chr",chr1, '_', cell_type, '.txt'), 
                    header = T, as.is = T, sep = '\t')
  resA = rbind(resA, resi)
  }
  write.table(resA, 
              file = paste0(res_folder, 
                            "../results_all/maf0.05_log_ratio/tcga_gwas_maf0.05_", 
                            cell_type, '.txt'), 
              quote = FALSE, 
              sep = "\t", row.names = F, col.names = TRUE)
  
}


## --------------------------------------------------------
## Manhattan plot
## --------------------------------------------------------
indir = "/fh/fast/sun_w/licai/COAD/data_GWAS/tcga_gwas/results_all/maf0.05_log_ratio/"
setwd(indir)
MakeManhattans <- function(dir, files, cols=c("gray10", "gray60"), 
                           pvalToPlot="p_value", png_name=NULL, 
                           labelSigs=TRUE, highlight=NULL, ...){
  library(qqman)
  res_all_afterQC <- read.table(files, sep='', header=T, check.names=F, 
                                colClasses = c(rep('NULL', 3), 
                                               rep(NA, 3), rep('NULL', 3)))
  res_all_afterQC$chr <- as.numeric(sub("^chr","",res_all_afterQC$chr))
  res_all_afterQC <- res_all_afterQC[!is.na(res_all_afterQC$pos),]
  columnToUse <- which(pvalToPlot == colnames(res_all_afterQC))
  mplot <- data.frame(SNP=rownames(res_all_afterQC), 
                      CHR=res_all_afterQC$chr, BP=res_all_afterQC$pos, 
                      P=as.numeric(res_all_afterQC[,columnToUse]), 
                      stringsAsFactors=FALSE)
  if(anyNA(mplot))
    mplot <- mplot[-which(apply(mplot, 1, anyNA)),]
  matchedHighlight <- NULL
  if(!is.null(highlight)){
    if(length(grep("rs", highlight)) == length(highlight))
      highlight <- find_rs_to_gigs(highlight)[,2]
    posnames <- sapply(mplot$SNP, function(x) 
      return(strsplit(x, split="_")[[1]][1]), USE.NAMES=FALSE)
    matchedHighlight <- mplot$SNP[which(posnames %in% highlight)]
  }
  if(is.null(png_name))
    png_name <- paste0("Manhattan_", gsub(".txt", "", files), ".png")
  png(file=paste0(dir, png_name), width=8, height=8, units="in", res=600, 
      pointsize=12)
  manhattan(mplot, col=cols, cex.axis=1, bty='o', highlight=highlight, ...)
  box()
  dev.off()
  #return(res_all_afterQC[oversigs,])
}
for(i in list.files(indir))
  MakeManhattans(dir="/fh/fast/sun_w/licai/cell_type_association_summary/act/TCGA_COAD/figures/", 
                 files =  i, 
                 labelSigs=FALSE, highlight=NULL)

## --------------------------------------------------------
## histgram of top findings
## --------------------------------------------------------

qqp <- function(pvals, main, confidence=.95, cutoff=1){
  
  alpha = 1-confidence
  n     = length(pvals)
  
  pvals[is.na(pvals)]=1
  pvals=sort(pvals)
  
  k=c(1:n)
  
  lower = cutoff*qbeta(alpha/2, k, n+1-k)
  upper = cutoff*qbeta((1-alpha/2), k, n+1-k)
  
  expected = cutoff*k/(n+1)
  n0 = length(which(pvals ==0))
  
  if(n0 > 0){
    warning(sprintf("there are %d p-values being 0\n", n0))
  }
  
  biggest= max(-log10(pvals[which(pvals > 0)]), -log10(expected))
  
  plot(-log10(expected), -log10(pvals), xlim=c(0,biggest),
       ylim=c(0,biggest), pch=20, xlab="-log10(expected p-value)",
       ylab="-log10(observed p-value)", cex=0.6, bty="n", main=main)
  
  lines(-log10(expected), -log10(lower), lty=2)
  lines(-log10(expected), -log10(upper), lty=2)
  abline(0, 1, col="red", lwd=1)
}

pvals1 = system("cut -f 9 microbiome_tcga_22ct_dG0.5_E_GWAS_all.txt", 
                intern = T)[-1]
pvals1 = as.numeric(pvals1)
png(sprintf('/fh/fast/sun_w/licai/cell_type_association_summary/act/TCGA_COAD/figures/qqplot_microbiome.png'
            ))
par(mfrow=c(1,2))
qqp(pvals1,  main = "", confidence=.95, cutoff=1)
pvals1 = system("cut -f 9 microbiome_tcga_22ct_kernel_alpha0.5_GWAS_all.txt",
                intern = T)[-1]

pvals1 = as.numeric(pvals1)
qqp(pvals1,  main = "", confidence=.95, cutoff=1)
dev.off()

for(i in list.files()){
  res = read.csv(i, header=T, sep="\t", stringsAsFactors = F)
  # hist(res$p_value, main = gsub(".txt", "", i))
  title1 = gsub("tcga_gwas_maf0.05_", "" ,gsub(".txt", "", i)) 
  png(sprintf('/fh/fast/sun_w/licai/cell_type_association_summary/act/TCGA_COAD/figures/hist_%s.png'
              , title1))
  # res = res[which(res$p_value < 0.001), ]
  qqp(res$p_value,  main = title1, confidence=.95, cutoff=1)
  dev.off()
}

## --------------------------------------------------------
## top 20k findings command line
## --------------------------------------------------------
# remove all the spaces in file names
# for f in *\ *; do mv "$f" "${f// /_}"; done

# top 20k res
# for file in *; 
# do sort -gk4 "$file"  | head -20001  > top20k_"$file";
# done;

# sig res
# for file in *;
# do awk '$4 < 0.00000005' "$file"  | head -20001  > sig_"$file";
# done;


load("/shared/cs_researcher/Newcomb_P/Molecular Correlates_ISACC/Molecular Correlates/Prognostic Modeling/licai_working/GWAS/N0147_C08_GWAS/N0147_C08_survival_res_all_MAF0.01.Rdata")

for (i in list.files(pattern = '^top')){
  res = read.csv(i, header=T, sep="\t", stringsAsFactors = F)
  # rs_id = result_rs$rsid[match(rownames(res), result_rs$rs)]
  title1 = gsub("top20k_tcga_gwas_maf0.05_", "", i)
  title1 = gsub('.txt', '', title1)
  title1 = gsub("_", " ", title1)
  pval = res_all_afterQC$p.fixed[match(res$rsid, res_all_afterQC$rs)]
  pdf(sprintf('/fh/fast/sun_w/licai/cell_type_association_summary/act/TCGA_COAD/figures/survival_gwas_pvalue_%s.pdf'
              , title1))
  hist(pval, main = title1, xlab = 'p-value in survival GWAS')
  dev.off()
  res$p.fixed_survival_gwas = pval
  res <- res[order(res$p.fixed_survival_gwas),]
  write.table(res, file = i,  quote = FALSE, 
              sep = "\t", row.names = T, col.names = TRUE)
}

