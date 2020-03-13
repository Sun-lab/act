## --------------------------------------------------------
## plot
## --------------------------------------------------------
setwd("/fh/fast/sun_w/licai/COAD/data_lineage_dist/tcga/maf005/revision/")
list.dirs()
# combine data
# head -1 microbiome_L2_ksigma_GWAS_cov_chr1.txt > ../microbiome_L2_ksigma_GWAS_cov.txt; tail -n +2 -q microbiome_L2_ksigma_GWAS_cov_chr*.txt >> ../microbiome_L2_ksigma_GWAS_cov.txt
# head -1 microbiome_cor_GWAS_cov_chr1.txt > ../microbiome_cor_GWAS_cov.txt ; tail -n +2 -q microbiome_cor_GWAS_cov_chr*.txt >> ../microbiome_cor_GWAS_cov.txt
# head -1 microbiome_Gunifrac_k40sigma0.2_GWAS_cov_chr1.txt > ../microbiome_Gunifrac_k40sigma0.2_GWAS_cov.txt ; tail -n +2 -q microbiome_Gunifrac_k40sigma0.2_GWAS_cov_chr*.txt >> ../microbiome_Gunifrac_k40sigma0.2_GWAS_cov.txt
# head -1 microbiome_L2_k40sigma0.2_GWAS_cov_chr1.txt > ../microbiome_L2_k40sigma0.2_GWAS_cov.txt ; tail -n +2 -q microbiome_L2_k40sigma0.2_GWAS_cov_chr*.txt >> ../microbiome_L2_k40sigma0.2_GWAS_cov.txt
# head -1 microbiome_L2_k40sigma0.25_GWAS_cov_chr1.txt > ../microbiome_L2_k40sigma0.25_GWAS_cov.txt ; tail -n +2 -q microbiome_L2_k40sigma0.25_GWAS_cov_chr*.txt >> ../microbiome_L2_k40sigma0.25_GWAS_cov.txt
# head -1 microbiome_L2_k40sigma0.5_GWAS_cov_chr1.txt > ../microbiome_L2_k40sigma0.5_GWAS_cov.txt ; tail -n +2 -q microbiome_L2_k40sigma0.5_GWAS_cov_chr*.txt >> ../microbiome_L2_k40sigma0.5_GWAS_cov.txt
# head -1 microbiome_Gunifrac_k40sigma0.5_GWAS_cov_chr1.txt > ../microbiome_Gunifrac_k40sigma0.5_GWAS_cov.txt ; tail -n +2 -q microbiome_Gunifrac_k40sigma0.5_GWAS_cov_chr*.txt >> ../microbiome_Gunifrac_k40sigma0.5_GWAS_cov.txt

## --------------------------------------------------------
## top 20k findings 
## --------------------------------------------------------

# top 20k res
# for file in microbiome*;
# do sort -gk9 "$file"  | head -20001  > top_sig/top20k_"$file";
# done;

# sig res
# for file in microbiome*;
# do awk '$9 < 0.00000005' "$file"  | head -20001  > top_sig/sig_"$file";
# done;

setwd("top_sig/")
load("/shared/cs_researcher/Newcomb_P/Molecular Correlates_ISACC/Molecular Correlates/Prognostic Modeling/licai_working/GWAS/N0147_C08_GWAS/N0147_C08_survival_res_all_MAF0.01.Rdata")

res_all = NULL
for (fli in list.files(pattern = 'top20k')){
  res = read.table(fli, header=T, sep="\t", stringsAsFactors = F)
  pval = res_all_afterQC$p.fixed[match(res$rsid, res_all_afterQC$rs)]
  title1 = gsub("microbiome_", "", fli)
  title1 = gsub('_GWAS_cov.txt', '', title1)
  title1 = gsub('L2', 'Euclidean', title1)
  title1 = gsub('ksigma', '', title1)
  pdf(sprintf('/fh/fast/sun_w/licai/cell_type_association_summary/act/TCGA_COAD/figures/survival_gwas_pvalue_%s.pdf'
              , title1),  pointsize= 20)
  par(mar=c())
  title1 = gsub('top20k', '', title1)
  title1 = gsub("_", " ", title1)
  title1 = gsub('k', ' GK ', title1)
  title1 = gsub('sigma', ',', title1)
  title1 = gsub("cor", "correlation", title1)
  
  hist(pval, main = title1, xlab = 'p-value in survival GWAS')
  title1 = gsub("_", " ", title1)

  dev.off()
  
  res$p.fixed_survival_gwas = pval
  res <- res[order(res$p.fixed_survival_gwas),]
  tot = sum(!is.na(res$p.fixed_survival_gwas))
  
  obs = sum(res$p.fixed_survival_gwas< 0.05, na.rm = T)
  expe = tot * 0.05
  btest =   binom.test(obs, tot, 0.05, "greater")$p.value
  res_all = rbind(res_all, c(title1, tot, obs, expe, 
                             formatC(btest, format= "e")))
  # write.table(res, file = i,  quote = FALSE,
  #             sep = "\t", row.names = T, col.names = TRUE)
}

colnames(res_all) = c("Distance", "Total SNPs", "Nsignificant SNPs", "expected", 
                      "binomial p.value")
write.csv(res_all, file = "/fh/fast/sun_w/licai/cell_type_association_summary/act/TCGA_COAD/revision/survival_pvalue_enrichment_test.csv", 
          quote = F, 
          col.names = T, row.names = F)
## --------------------------------------------------------
## qq plot 
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

setwd("../")

for(fli in list.files(pattern = '^microbiome_')){
  res = system(sprintf("cut -f 9 %s", fli),
               intern = T)[-1]
  res = as.numeric(res)
  title1 = gsub("microbiome_", "", fli)
  title1 = gsub("cor", "correlation", title1)
  title1 = gsub('_GWAS_cov.txt', '', title1)
  title1 = gsub('L2_', 'Euclidean', title1)
  title1 = gsub('ksigma', '', title1)
  cat(title1,"\n")
  png(sprintf('/fh/fast/sun_w/licai/cell_type_association_summary/act/TCGA_COAD/figures/qqplot_%s.png'
              , title1))
  title1 = gsub('_', ' ', title1)
  title1 = gsub('k', ' GK ', title1)
  title1 = gsub('sigma', ',', title1)
  cat(title1,"\n")
  qqp(res, main = title1)
  dev.off()
}

setwd("/fh/fast/sun_w/licai/cell_type_association_summary/act/TCGA_COAD/figures/")

q('no')
