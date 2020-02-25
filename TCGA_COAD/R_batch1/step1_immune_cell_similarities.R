
# -------------------------------------------------------------------
# read in data
# LM22-ref-sample.txt: gene expression data of purified samples
# LM22-classes.txt: matching between gene expression reference sample 
#                   vs. cell types. some samples do not match to any 
#                   one of the 22 cell types, we will remove them.
# -------------------------------------------------------------------

rDat = read.table("../data/LM22-ref-sample.txt", as.is=TRUE, sep="\t",
                  header=TRUE, quote="")
dim(rDat)
rDat[1:2,1:5]

cDat = read.table("../data/LM22-classes.txt", as.is=TRUE, sep="\t",
                  header=FALSE, quote="")
dim(cDat)
cDat[1:2,1:5]

rownames(rDat) = rDat$Relabel
rDat = data.matrix(rDat[,-1])

rownames(cDat) = cDat$V1
cDat = data.matrix(cDat[,-1])

table(colSums(cDat))
w2kp = which(colSums(cDat) > 0)

rDat = rDat[,w2kp]
cDat = cDat[,w2kp]

dim(rDat)
rDat[1:2,1:5]

dim(cDat)
cDat[1:2,1:5]

apply(cDat, 1, table)

# -----------------------------------------------------------------
# simplify sample name
# -----------------------------------------------------------------

colnames(rDat)

sms = gsub("..HG.U133A...", ".", colnames(rDat), fixed=TRUE)
sms = gsub("_U133A.", "", sms)
sms = gsub("_U133A_", "_", sms)
sms = gsub(".HG.U133A.", ".", sms, fixed=TRUE)
sms = gsub("^X\\.", "", sms, perl=TRUE)
sms = gsub("..", ".", sms, fixed=TRUE)
sms = gsub(".$", "", sms, perl=TRUE)
sms = gsub("^CXCR5hiICOShi\\d\\.", "", sms, perl=TRUE)
sms = gsub(".biological.", ".", sms, fixed=TRUE)

sms

colnames(rDat) = sms

# -----------------------------------------------------------------
# get cell type information
# -----------------------------------------------------------------

info = data.frame(sampleID=colnames(rDat), cellType=rep(NA, ncol(rDat)))
dim(info)
info[1:2,]

for(i in 1:nrow(info)){
  wwi = which(cDat[,i] == 1)
  if(length(wwi) != 1){
    stop("number of cell type is not 1\n")
  }
  info$cellType[i] = rownames(cDat)[wwi]
}

info
table(info$cellType)

ct = info$cellType
ct[which(info$cellType=="Mast cells resting")]   = "Mast resting"
ct[which(info$cellType=="Mast cells activated")] = "Mast activated"
ct[which(info$cellType=="T cells gamma delta")]  = "T gd"
ct[which(info$cellType=="T cells CD8")]  = "CD8+ T"
ct[which(info$cellType=="T cells CD4 memory resting")] = "CD4+ T memory resting"
ct[which(info$cellType=="T cells CD4 memory activated")]  = "CD4+ T memory activated"
ct[which(info$cellType=="NK cells resting")]    = "NK resting"
ct[which(info$cellType=="NK cells activated")]  = "NK activated"
ct[which(info$cellType=="B cells naive")]   = "B naive"
ct[which(info$cellType=="B cells memory")]  = "B memory"
ct[which(info$cellType=="Plasma cells")]    = "Plasma"
ct[which(info$cellType=="Dendritic cells resting")]   = "Dendritic resting"
ct[which(info$cellType=="Dendritic cells activated")] = "Dendritic activated"
ct[which(info$cellType=="T cells follicular helper")] = "Tfh"
ct[which(info$cellType=="T cells CD4 naive")] = "CD4+ T"
ct[which(info$cellType=="T cells regulatory (Tregs)")] = "Tregs"

info$label = ct

uct = unique(ct)

panT = c("CD4+ T", "CD8+ T", "CD4+ T memory activated", "CD4+ T memory resting")
panT = c(panT, "T gd", "Tfh", "Tregs")
panB = c("B memory", "B naive", "Plasma")
panNK  = c("NK resting", "NK activated")
panMono  = c("Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2")
panMono  = c(panMono, "Dendritic resting",  "Dendritic activated")
panGranu = c("Mast resting", "Mast activated", "Neutrophils", "Eosinophils")

setdiff(uct, c(panT, panB, panNK, panMono, panGranu))

newCt = c(panGranu, panMono, panNK, panB, panT)
newCt

# -----------------------------------------------------------------
# order the data based on new order of cell types, to faciliate
# some plots to separate cell types
# -----------------------------------------------------------------

newRDat = NULL

for(ct1 in newCt){
  ww1 = which(info$label == ct1)
  newRDat = cbind(newRDat, rDat[,ww1])
}

dim(newRDat)
newRDat[1:2,1:5]

newInfo = info[match(colnames(newRDat), info$sampleID),]

dim(newInfo)
newInfo[1:2,]

info = newInfo
rDat = newRDat

# ------------------------------------------------------------
# add study information, based table S2 of CIBERSORT paper
# ------------------------------------------------------------

table(info$label)

studies = rep("", nrow(info))

wAbbas    = grep("GSE22886", info$sampleID, fixed=TRUE)
studies[wAbbas] = "Abbas"

wChtanova = grep("^A_", info$sampleID, perl=TRUE)
studies[wChtanova] = "Chtanova"

wCho      = grep("GSE7138", info$sampleID, fixed=TRUE)
studies[wCho] = "Cho"

wOcklenburg   = grep("GSE4527", info$sampleID, fixed=TRUE)
studies[wOcklenburg] = "Ocklenburg"

info[which(studies==""),]
studies[which(studies=="")] = "Rasheed"

table(studies)

info$study = studies
table(info$study, info$label)

# ------------------------------------------------------------
# check distribution for genes selected for LM22
# ------------------------------------------------------------

lDat = read.table("../data/LM22.txt", as.is=TRUE, sep="\t",
                  header=TRUE, quote="")
dim(lDat)
lDat[1:2,1:5]

rownames(lDat) = lDat$Gene.symbol
lDat = data.matrix(lDat[,-1])
dim(lDat)
lDat[1:2,1:5]

names(lDat)
lDatLog = log10(lDat)

uct  = colnames(lDat)
cols = rainbow(length(uct), start=0.2, end=1)

cMean = apply(lDat, 2, mean)
cSD   = apply(lDat, 2, sd)

pdf("../figures/LM22_547_genes_summary.pdf", width=8, height=4)
par(mfrow=c(1,2), mar=c(5,4,1,1), bty="n", cex=0.9)
plot(cMean, cSD, type="n", bty="n", xlab="mean", ylab="sd", bty="n")

for(i in 1:length(uct)){
  ct1 = uct[i]
  points(cMean[i], cSD[i], pch=i, col=cols[i])
}

legend("topleft", legend=uct[1:11], pch=1:11, col=cols[1:11], 
       bty="n", cex=0.6)
legend("bottomright", legend=uct[12:22], pch=12:22, col=cols[12:22], 
       bty="n", cex=0.6)

plot(density(lDatLog[,1]), col=cols[1], lwd=2, main="", ylim=c(0,0.6), 
     xlab="log10(expression)")
for(i in 2:22){lines(density(lDatLog[,i]),lwd=2,col=cols[i])}

dev.off()

# ------------------------------------------------------------
# check distribution
# ------------------------------------------------------------

rDatLog = log10(rDat)

uct  = unique(info$label)

pdf("../figures/LM22_ref_summary.pdf", width=8, height=4)

par(mfrow=c(1,2), mar=c(5,4,1,1), bty="n")
cMean = apply(rDat, 2, mean)
cSD   = apply(rDat, 2, sd)

plot(cMean, cSD, type="n", bty="n", xlab="mean", ylab="sd")

cols = rainbow(length(uct), start=0.2, end=1)

for(i in 1:length(uct)){
  ct1 = uct[i]
  ww1 = which(info$label==ct1)
  points(cMean[ww1], cSD[ww1], pch=i, col=cols[i])
}

legend("topleft", legend=uct[1:15], pch=1:15, col=cols[1:15], 
       bty="n", cex=0.6)
legend("bottomright", legend=uct[16:22], pch=16:22, col=cols[16:22], 
       bty="n", cex=0.6)

summary(cMean)
summary(cSD)

plot(density(rDatLog[,1]),col=cols[1],lwd=2, ylim=c(0,0.9), main="",
     xlab="log10(gene expression")
for(i in 2:ncol(rDatLog)){lines(density(rDatLog[,i]),lwd=2,col=cols[i])}

dev.off()

# ------------------------------------------------------------
# select probes that are differentially expressed across 
# cell types
# ------------------------------------------------------------

labels    = uct
rDatLog   = log(rDat)


write.table(rDatLog, file = "../data/LM22-ref-sample_logged.txt",
            append = FALSE, quote = FALSE, sep = "\t", row.names = T, 
            col.names = T)


ct  = as.factor(info$label)
pvs = rep(NA, nrow(rDatLog))

for(i in 1:nrow(rDatLog)){
  
  if(i %% 5000 == 1){cat (i, date(), "\n")}
  
  yi = rDatLog[i,]
  li = lm(yi ~ ct)
  ai = anova(li)
  
  pvs[i] = ai$Pr[1]
}

table(pvs < 1e-10)
table(pvs < 1e-20)

rDatLog   = rDatLog[which(pvs < 1e-10),]
dim(rDatLog)

# ------------------------------------------------------------
# do PCAs, original data
# ------------------------------------------------------------

dat4Pr = rDatLog - rowMeans(rDatLog, na.rm=TRUE)

dat4Pr[is.na(dat4Pr)] = 0
covdat = t(dat4Pr) %*% dat4Pr / nrow(dat4Pr)
dim(covdat)
prdat  = eigen(covdat)

prdat$values[1:20]/sum(prdat$values)

PC1 =  prdat$vectors[,1]
PC2 =  prdat$vectors[,2]
PC3 =  prdat$vectors[,3]
PC4 =  prdat$vectors[,4]

pdf("../figures/expression_PCs_logged_filtered_data.pdf", width=10.5, height=7)
par(mar=c(5,4,1,1), mfrow=c(2,3), bty="n", cex=0.9)
barplot(prdat$values[1:20], main="", xlab="Index", ylab="Eigen-value")

legend("topright", legend=labels[1:10], col=cols[1:10], bty="n",
       pch=1:10)

plot(0:1, 0:1, xaxt="n", xlab="", yaxt="n", ylab="", type="n")
legend("topleft", legend=labels[11:22], col=cols[11:22], bty="n",
       pch=11:22)

plot(PC1, PC2,  bty="n", cex=0.9, type="n")
for(j in 1:length(labels)){
  wj = which(info$label == labels[j])
  points(PC1[wj], PC2[wj], col=cols[j], pch=j, cex=0.9)
}

plot(PC1, PC3,  bty="n", cex=0.9, type="n")
for(j in 1:length(labels)){
  wj = which(info$label == labels[j])
  points(PC1[wj], PC3[wj], col=cols[j], pch=j, cex=0.9)
}

plot(PC1, PC4,  bty="n", cex=0.9, type="n")
for(j in 1:length(labels)){
  wj = which(info$label == labels[j])
  points(PC1[wj], PC4[wj], col=cols[j], pch=j, cex=0.9)
}

plot(PC2, PC3,  bty="n", cex=0.9, type="n")
for(j in 1:length(labels)){
  wj = which(info$label == labels[j])
  points(PC2[wj], PC3[wj], col=cols[j], pch=j, cex=0.9)
}


dev.off()

# ------------------------------------------------------------
# write out data
# ------------------------------------------------------------

write.table(rDatLog, file = "../data/LM22-ref-sample_logged_filtered.txt", 
            append = FALSE, quote = FALSE, sep = "\t", row.names = T,
            col.names = T)

write.table(info, file = "../data/LM22-ref-sample_info.txt", 
            append = FALSE, quote = FALSE, sep = "\t", row.names = F, 
            col.names = T)

sessionInfo()

quit(save="no")


