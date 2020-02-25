rm(list = ls())

count_file = "../data/expression_v2_counts.txt"
info_file  = "../data/expression_v2_info.txt"
samp_file  = "../data/expression_v2_sample.txt"

# ------------------------------------------------------------
# read in data
# ------------------------------------------------------------

datE = read.table(file = count_file, sep = "\t",
header = TRUE, as.is = TRUE)
dim(datE)
datE[1:2, 1:5]

infoE = read.table(file = info_file, sep = "\t",
                    header = TRUE, as.is = TRUE,  quote="")
dim(infoE)
infoE[1:2, 1:5]
table(infoE$chr)

table(rownames(datE) == infoE$gene)

datE = data.matrix(datE)

dim(datE)
datE[1:2, 1:5]

coad = read.table(samp_file, sep="\t",
            as.is=TRUE, header=TRUE)
dim(coad)
coad[1:2,]

# ----------------------------------------------------------------------------
# Normalize gene expression by read-depth
# ----------------------------------------------------------------------------

load('../data/Gene_Lengths.RData')
class(GeneLengths.Mat)
head(GeneLengths.Mat)


table(GeneLengths.Mat$Gencode.ID %in% rownames(datE)) 
length(rownames(datE))

GeneLength = GeneLengths.Mat[match(rownames(datE), GeneLengths.Mat$Gencode.ID), ]
dim(GeneLength)

#log(T_{ij} / read-depth of sample i / gene length of gene j)
tot = colSums(datE)
s75 = apply(datE, 2, quantile, prob = 0.75)

cor(tot, s75)

datE = log(t(t(datE + 1) / s75) / GeneLength$Exonic) 
dim(datE)


# ------------------------------------------------------------
# read in reference data
# ------------------------------------------------------------

lm22 = read.table("../data/LM22-ref-sample_logged.txt",  sep = "\t",
  header = TRUE, as.is = TRUE)
dim(lm22)
lm22[1:2,1:2]
lm22 = data.matrix(lm22)

sam = read.table("../data/LM22-ref-sample_info.txt" , sep = "\t",
  header = TRUE, as.is = TRUE)
dim(sam)
sam[1:2,]

table(sam$label)

# ------------------------------------------------------------
# take intersection
# ------------------------------------------------------------

genes = intersect(infoE$hgnc_symbol, rownames(lm22))
length(genes)
length(unique(genes))
genes[c(1:5, length(genes))]

mat1 = match(genes, rownames(lm22))
mat2 = match(genes, infoE$hgnc_symbol)

lm22 = lm22[mat1,]
dim(lm22)
lm22[1:2,1:2]

infoE = infoE[mat2,]
dim(infoE)
infoE[1:2, 1:5]

datE = datE[mat2,]
dim(datE)
rownames(datE) = genes
datE[1:2, 1:5]

table(infoE$hgnc_symbol == rownames(lm22))

# ------------------------------------------------------------
# write out data
# ------------------------------------------------------------

write.table(datE, file = "../data/expresssion_COAD_data.txt", append = FALSE,
quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

q(save = "no")