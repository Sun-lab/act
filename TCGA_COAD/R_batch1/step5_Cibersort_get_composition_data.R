library(nnls)
library(quantreg)

# ------------------------------------------------------------
# run CIBERSORT
# ------------------------------------------------------------

source('CIBERSORT.R')

CiberRun <- function(X, Y, QN = T){
# read in data

  X = data.matrix(X)
  Y = data.matrix(Y)

# anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y = exp(Y)}

# quantile normalization of mixture file
  if(QN == TRUE){
    tmpc = colnames(Y)
    tmpr = rownames(Y)
    Y = normalize.quantiles(Y)
    colnames(Y) = tmpc
    rownames(Y) = tmpr
  }

# intersect genes
  Xgns = row.names(X)
  Ygns = row.names(Y)
  genes = sort(intersect(Xgns, Ygns))
  length(genes)

  Y = Y[match(genes,Ygns),]
  X = X[match(genes,Xgns),]

  dim(X)
  dim(Y)

# standardize sig matrix
  Xsd = (X - mean(X)) / sd(as.vector(X))

# iterate through mixtures

  w.svm = matrix(nrow=ncol(Y), ncol=ncol(X))
  colnames(w.svm) = colnames(X)
  rownames(w.svm) = colnames(Y)

  for(itor in 1:ncol(Y)){
  
    y = Y[,itor]
    ysd = (y - mean(y)) / sd(y)
  
  # run SVR core algorithm
    result = CoreAlg(Xsd, ysd)

  # get results
    w.svm[itor,] = result$w
  }
  return(w.svm)
}

# ------------------------------------------------------------
# run CIBERSORT using LM22.txt and original gene expression data
# ------------------------------------------------------------

sig_matrix   = '../data/LM22.txt'
mixture_file = '../data/expresssion_COAD_data.txt'

X = read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
Y = read.table(mixture_file, header=T, sep="\t", row.names=1,check.names=F)

lm22.ori = CiberRun(X, Y, QN = F)
lm22.ori[1:5,1:5]
write.table(lm22.ori, file = "../data/COAD_composition.txt", 
            append = FALSE,quote = FALSE, 
            sep = "\t", row.names = TRUE, col.names = TRUE)

q(save = 'no')