rm(list=ls())

library(gplots)
library(RColorBrewer)
library(edgeR)
library(limma)
library(openxlsx)

### read metadata
metadata <- read.csv("~/metadata.csv", header = T)

### read counts table (Pathway database, Virulence factor database)
countdata <- read.csv("~/countdata.csv", header = T, row.names = "Product") #row.names = dataset rownames


## Create DGList object
x <- DGEList(counts=countdata)


##---------------------------------------------------
## Building the experimental design information  ----
##---------------------------------------------------
# (metadata[,10]) <- grouping variable
group <- as.factor(as.character(metadata[,10]))
x <- DGEList(counts=countdata, group = group)
x$samples

x$samples$group <- group
id <- metadata[,2]
x$samples$id <- id

### Some boxplot
nomcol <- x$samples$id
par(mar=c(10,2,1,1))
boxplot(cpm(x, log=T), main="cpms in log scale (not normalized)", names= nomcol, las=2, colours = x$samples$group)


###-----------------------------------
### Normalization
### ----------------------------------

x <- edgeR::calcNormFactors(x, method = "TMM")
boxplot(cpm(x, log = T), main="cpms in log scale (TMM normalization)", names= nomcol, las=2)

#### PCA's analysis

col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "PuRd")
col.group <- as.character(col.group)

lcpm <- cpm(x, log = TRUE) 
par(mar=c(2,2,2,2))

prcTMM <- prcomp(lcpm)

### graphics

### TMM normalization
plot(prcTMM$rotation[,1], prcTMM$rotation[,2], col=col.group, main = "PCA (TMM normalization)")
text(prcTMM$rotation[,1],prcTMM$rotation[,2], col= col.group, labels = nomcol)

### PCA NW
prcNW <- prcomp(lcpm[,metadata$percentil_merged=="NW"])
plot(prcNW$rotation[,1], prcNW$rotation[,2], col=col.group, main="PCA NW community")
text(prcNW$rotation[,1],prcNW$rotation[,2], col= col.group, labels = nomcol[c(metadata$percentil_merged =="NW")])

### PCA OWOB
prcOB <- prcomp(lcpm[,metadata$percentil_merged=="OB"])
plot(prcOB$rotation[,1], prcOB$rotation[,2], col=col.group, main="PCA OW-OB community")
text(prcOB$rotation[,1],prcOB$rotation[,2], col= col.group, labels = nomcol[c(metadata$percentil_merged =="OB")])


###-------------------
### building the design
###--------------------

design <- model.matrix(~0+group, data=x$samples)
colnames(design) <- levels(x$samples$group)

###----------------------------------------
### estimating dispersions
###----------------------------------------

x <- estimateDisp(x, design)

x <- estimateGLMTagwiseDisp(x, design)


###-----------------------------------
### Testing for DE genes
###-----------------------------------

Qfit <- glmQLFit(x, design)
## likelihood ratio test
fit <- glmFit(x, design)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

### defining contrast

matrix.contr <- makeContrasts(NW.vs.OB = NW-OB,
                              OB.vs.NW = OB-NW,
                              levels = design)

### comparisons

qlf.NWvsOB <- glmQLFTest(Qfit, contrast = matrix.contr[,"NW.vs.OB"])
qlf.OBvsNW <- glmQLFTest(Qfit, contrast = matrix.contr[,"OB.vs.NW"])
lf.NWvsOB <- glmLRT(fit, contrast = matrix.contr[, "NW.vs.OB"])
lf.OBvsNW <- glmLRT(fit, contrast = matrix.contr[,"OB.vs.NW"])


### MD Plots
plotMD(qlf.NWvsOB, main = "NW vs OW-OB")
abline(h=c(-1, 1), col="blue")
plotMD(qlf.OBvsNW, main = "OW-OB vs NW")
abline(h=c(-1, 1), col="blue")
plotMD(lf.NWvsOB, main = "NW vs OW-OB")
abline(h=c(-1, 1), col="blue")
plotMD(lf.OBvsNW, main = "OW-OB vs NW")
abline(h=c(-1, 1), col="blue")

### DE genes summary
summary(decideTests(qlf.NWvsOB))
summary(decideTests(qlf.OBvsNW))
summary(decideTests(lf.NWvsOB))
summary(decideTests(lf.OBvsNW))

### Write tables

## Contrast 
qlf.NWvsOB <- topTags(qlf.NWvsOB, n = Inf)
qlf.OBvsNW <- topTags(qlf.OBvsNW, n = Inf)
lf.NWvsOB <- topTags(lf.NWvsOB, n = Inf)
lf.OBvsNW <- topTags(lf.OBvsNW, n = Inf)

#write.csv(qlf.NWvsOB, "~/qlf.NWvsOB.csv")
#write.csv(qlf.OBvsNW, "~/qlf.OBvsNW.csv")
#write.csv(lf.NWvsOB, "~/lf.NWvsOB.csv")
#write.csv(lf.OBvsNW, "~/lf.OBvsNW.csv")

