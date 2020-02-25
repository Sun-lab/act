# Workflow

**Note that the files (LM22-ref-sample.txt, LM22-classes.txt, LM22.txt ) and source code of CIBERSORT can be downloaded from https://cibersort.stanford.edu/ after register. The license of CIBERSORT does not allow us to share these files in this public repository.**

step1: check the gene expression data of purified sample 
	Input: LM22-ref-sample.txt & LM22-classes.txt & LM22.txt 
	
step2: collect gene expression data htseq.counts files
	Input: meta information & clinical information & htseq.counts
	
step3: prepare gene expression data - filter out low expressed genes
	Input: expression_data_v1.txt & gencode.v22.genes.txt & metadata_clinic.txt
	Based on: step2
	
step4: Normalize gene expression by read-depth
	Input: expression_v2_counts.txt & expression_v2_info.txt & expression_v2_sample.txt & LM22-ref-sample_logged.txt & LM22-ref-sample_info.txt
	Based on: step3 & step1
	
step5: Run CIBERSORT
	Input: LM22.txt & expresssion_COAD_data.txt
	Based on: CIBERSORT.R & step4 
	
step6: survival analysis (a: cell type specific, b: cell type composition)
	Input: expression_v2_sample.txt & COAD_composition.txt
	Based on: step5 & Unifrac_dist.R(b) & immune_cell_lineage_tree_434_gene_22ct.tre(b)
	
step7: GWAS
	Input: PCA files & sample files & genotype files & expression_v2_sample.txt
	Based on: step5 & microbiomeGWAS.R & immune_cell_lineage_tree_434_gene_22ct.tre 


