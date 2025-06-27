source("utils.R")

# Import Data ----

cts <- read.delim("./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3counts_transcript.tsv")
geneLookup <- data.frame(TXNAME = cts$TXNAME, GENEID = cts$GENEID)
rm(cts)

niLFC <- read.csv("./NODvsIRT/NODvsIRTlfc.csv")
niVar7 <- read.csv("./NODvsIRT/NODvsIRT_variants2.csv")
niSplice <- read.csv("./NODvsIRT/NODvsIRT_splice.csv")
niSpliceDE <- read.csv("./NODvsIRT/NODvsIRT_splice_DE.csv")

nmLFC <- read.csv("./NODvsMRT/NODvsMRTlfc.csv")
nmVar7 <- read.csv("./NODvsMRT/NODvsMRT_variants2.csv")

imLFC <- read.csv("./IRTvsMRT/IRTvsMRTlfc.csv")
imVar7 <- read.csv("./IRTvsMRT/IRTvsMRT_variants2.csv")

# Subsetting ----
## NODvsIRT----
# for NODvsIRT_splice, the subsetting was done manually. We can take this subset and apply it to all other samples

spliced_genes <- niSplice$GeneID

## NODvsMRT ----
nmSplice <- data.frame(matrix(nrow = length(spliced_genes), ncol = ncol(nmVar7)))
colnames(nmSplice) <- colnames(nmVar7)
count <- 0
for (i in spliced_genes) {
    count <- count + 1
    for (j in 1:nrow(nmVar7)) {
        if (i == nmVar7$GeneID[j]) {
            nmSplice[count,] <- nmVar7[j,]
        }
    }
}

write.csv(nmSplice, "./NODvsMRT/NODvsMRT_splice.csv")

## IRTvsMRT ----
imSplice <- data.frame(matrix(nrow = length(spliced_genes), ncol = ncol(imVar7)))
colnames(imSplice) <- colnames(imVar7)
count <- 0
for (i in spliced_genes) {
    count <- count + 1
    for (j in 1:nrow(imVar7)) {
        if (i == imVar7$GeneID[j]) {
            imSplice[count,] <- imVar7[j,]
        }
    }
}

write.csv(imSplice, "./IRTvsMRT/IRTvsMRT_splice.csv")

# De Correlation ----

## NODvsIRT ----
colnames(niLFC) <- c("BambuTx", colnames(niLFC[,2:ncol(niLFC)]))

count <- 0
niSpliceDE <- data.frame(matrix(ncol = ncol(niLFC)))
colnames(niSpliceDE) <- colnames(niLFC)
for (i in nmSplice$GeneID) {
    for (j in 1:nrow(niLFC)) {
        if (niLFC$GENEID[j] == i) {
            count = count + 1
            niSpliceDE[count,] <- niLFC[j,]
        }
    }
}

rownames(niSpliceDE) <- niSpliceDE$BambuTx
niSpliceDE <- niSpliceDE[,2:ncol(niSpliceDE)]
head(niSpliceDE)


## NODvsMRT ----
colnames(nmLFC) <- c("BambuTx", colnames(nmLFC[,2:ncol(nmLFC)]))

count <- 0
nmSpliceDE <- data.frame(matrix(ncol = ncol(nmLFC)))
colnames(nmSpliceDE) <- colnames(nmLFC)
for (i in nmSplice$GeneID) {
    for (j in 1:nrow(nmLFC)) {
        if (nmLFC$GENEID[j] == i) {
            count = count + 1
            nmSpliceDE[count,] <- nmLFC[j,]
        }
    }
}

rownames(nmSpliceDE) <- nmSpliceDE$BambuTx
nmSpliceDE <- nmSpliceDE[,2:ncol(nmSpliceDE)]
head(nmSpliceDE)

write.csv(niSpliceDE, "./NODvsMRT/NODvsMRT_splice_DE.csv")

## IRTvsMRT ----
# Isoform Props ----

library(foreach)
library(doParallel)

isoProps <- isoformProp2(cts)
write.csv(isoProps, "./isoform_proportions.csv")

isoProp10 <- isoProps[isoProps$prop > 0.10,]
isoProp05 <- isoProps[isoProps$prop > 0.05,]

NOD <- list(
    subset = "NOD",
    df = data.frame(TXNAME = cts$TXNAME, GENEID = cts$GENEID, CINOD = cts$CINOD_aln, DINOD = cts$DINOD_aln, EINOD = cts$EINOD_aln)
)
write.csv(NOD$df, "./subsets/NOD.csv", row.names = F)
IRT <- list(
    subset = "IRT",
    df = data.frame(TXNAME = cts$TXNAME, GENEID = cts$GENEID, CIRTJ = cts$CIRTJ_aln, DIRTJ = cts$DIRTJ_aln, EIRTJ = cts$EIRTJ_aln)
)
write.csv(IRT$df, "./subsets/IRT.csv", row.names = F)
MRT <- list(
    subset = "MRT",
    df = data.frame(TXNAME = cts$TXNAME, GENEID = cts$GENEID, CMRTJ = cts$CMRTJ_aln, DMRTJ = cts$DMRTJ_aln, EMRTJ = cts$EMRTJ_aln)
)
write.csv(MRT$df, "./subsets/MRT.csv", row.names = F)

allsets <- list(
    NOD = NOD, 
    IRT = IRT,
    MRT = MRT)



colnames(imLFC) <- c("BambuTx", colnames(imLFC[,2:ncol(imLFC)]))

count <- 0
imSpliceDE <- data.frame(matrix(ncol = ncol(imLFC)))
colnames(imSpliceDE) <- colnames(imLFC)
for (i in nmSplice$GeneID) {
    for (j in 1:nrow(imLFC)) {
        if (imLFC$GENEID[j] == i) {
            count = count + 1
            imSpliceDE[count,] <- imLFC[j,]
        }
    }
}

head(imSpliceDE)
rownames(imSpliceDE) <- imSpliceDE$BambuTx
imSpliceDE <- imSpliceDE[,2:ncol(imSpliceDE)]

write.csv(imSpliceDE, "./IRTvsMRT/IRTvsMRT_splice_DE.csv")


# Plotting ----

rownames(niLFC) <- niLFC[,1]

## NODvsIRT ----

layout(matrix(c(1,2), nrow = 2, ncol = 1))
for (i in niSplice$GeneID) {
    plotIsoform(gene = strsplit(i, ";")[[1]][2], annotation = "./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf", exon_marker = F, prop = isoProps)
    plotSpliceReg(niLFC, "NODvsIRT", i)
}

dev.off()
layout(matrix(c(1,2), nrow = 2, ncol = 1))
plotIsoform(gene = "MtrunA17_Chr8g0344721", annotation = "./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf", exon_marker = F, prop = isoProps)
plotSpliceReg(niLFC, "NODvsIRT", "gene_biotype mRNA; MtrunA17_Chr8g0344721")

dev.off()
layout(matrix(c(1,2), nrow = 1, ncol = 2))
plotIsoform(gene = "MtrunA17_Chr8g0344721", annotation = "./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf", exon_marker = T, prop = isoProps)
plotHeatmapIso(propTable, "gene_biotype mRNA; MtrunA17_Chr8g0344721", c("blue", "yellow"))

## NODvsMRT ----

dev.off()

rownames(nmLFC) <- nmLFC$BambuTx

for (i in nmSplice$GeneID) {
    plotSpliceReg(nmLFC, "NODvsMRT", i)
    plotIsoform(gene = strsplit(i, ";")[[1]][2], annotation = "./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf")
}

layout(matrix(c(1,2), nrow = 2, ncol = 1))
plotIsoform(gene = "MtrunA17_Chr2g0322801", annotation = "./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf", exon_marker = T)
plotSpliceReg(nmLFC, "NODvsMRT", "gene_biotype mRNA; MtrunA17_Chr2g0322801")

## IRTvsMRT ----

dev.off()

rownames(imLFC) <- imLFC$BambuTx

for (i in imSplice$GeneID) {
    plotSpliceReg(imLFC, "IRTvsMRT", i)
    plotIsoform(gene = strsplit(i, ";")[[1]][2], annotation = "./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf")
}

layout(matrix(c(1,2), nrow = 2, ncol = 1))
plotIsoform(gene = "MtrunA17_Chr2g0322801", annotation = "./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf", exon_marker = T)
plotSpliceReg(imLFC, "IRTvsMRT", "gene_biotype mRNA; MtrunA17_Chr2g0322801")

# Comparison ----

# Common Upregulated

library(edgeR)

fit <- readRDS("./edgeR_glm_fit.rds")
    
## NODvsIRT ----

NODvsIRT <- makeContrasts(nod-irt, levels = fit$design)
qlfNODvsIRT <- glmQLFTest(fit, contrast = NODvsIRT)
niDEX <- decideTests(qlfNODvsIRT, p = 0.05)
summary(niDEX)

niReg <- regTx(niDEX)
head(niReg)

## NODvsMRT ----

NODvsMRT <- makeContrasts(nod-mrt, levels = fit$design)
qlfNODvsMRT <- glmQLFTest(fit, contrast = NODvsMRT)
nmDEX <- decideTests(qlfNODvsMRT, p = 0.05)
summary(nmDEX)

nmReg <- regTx(nmDEX)
head(nmReg)

## IRTvsMRT ----

IRTvsMRT <- makeContrasts(irt-mrt, levels = fit$design)
qlfIRTvsMRT <- glmQLFTest(fit, contrast = IRTvsMRT)
imDEX <- decideTests(qlfIRTvsMRT, p = 0.05)
summary(imDEX)

imReg <- regTx(imDEX)
head(imReg)

# Comp ----

niimCompReg <- compareReg(niReg, imReg)
ninmCompReg <- compareReg(niReg, nmReg)
nmimCompReg <- compareReg(nmReg, imReg)

masterReg <- list(
    NODvsIRT = niReg,
    NODvsMRT = nmReg,
    IRTvsMRT = imReg,
    NIxIM = niimCompReg,
    NIxNM = ninmCompReg,
    NMxIM = nmimCompReg
)

saveRDS(masterReg, "./masterReg.rds")

## NODvsIRT X IRTvsMRT ----

png("./comp/NODvsIRTxIRTvsMRTup.png")
plotVen(length(niReg$upReg) - length(niimCompReg$up), length(niimCompReg$up), length(imReg$upReg) - length(niimCompReg$up),
        "Upregulated Genes Shared by NODvsIRT and IRTvsMRT",
        labl = "NODvsIRT",
        labc = "Both",
        labr = "IRTvsMRT"
)
dev.off()

niimRegUnique <- list(
    niRegUpUnique = data.frame(TXNAME = NA, GENEID = NA),
    imRegUpUnique = data.frame(TXNAME = NA, GENEID = NA),
    niimRegUpUnique = data.frame(TXNAME = NA, GENEID = NA),
    niRegDownUnique = data.frame(TXNAME = NA, GENEID = NA),
    imRegDownUnique = data.frame(TXNAME = NA, GENEID = NA),
    niimRegDownUnique = data.frame(TXNAME = NA, GENEID = NA)
    )

niimRegUpUnique <- list(
    niRegUpUnique = data.frame(TXNAME = NA, GENEID = NA),
    imRegUpUnique = data.frame(TXNAME = NA, GENEID = NA),
    niimRegUpUnique = data.frame(TXNAME = NA, GENEID = NA)
)

testTable <- regUnique(imReg$upReg, niReg$upReg, niimCompReg$up, niimRegUpUnique, geneLookup, 
                       setNames = c("imRegUp", "niRegUp", "niimCompReg")
)

imRegUpUnique <- c()
count <- 1
for (i in imReg$upReg) {
    test <- sum(niimCompReg$up == i)
    if (test == 1) {
        next
    } else {
        imRegUpUnique[count] <- i
        count <- count + 1
    }
}

count <- 1
for(i in imRegUpUnique) {
    for(j in 1:nrow(geneLookup)) {
        if (i == geneLookup$TXNAME[j]) {
            niimRegUnique$imRegUpUnique[count,] <- geneLookup[j,]
            count <- count + 1
        } else {
            next
        }
    }
}

niRegUpUnique <- c()
count <- 1
for (i in niReg$upReg) {
    test <- sum(niimCompReg$up == i)
    if (test == 1) {
        next
    } else {
        niRegUpUnique[count] <- i
        count <- count + 1
    }
}



png("./comp/NODvsIRTxIRTvsMRTdown.png")
plotVen(length(niReg$downReg) - length(niimCompReg$down), length(niimCompReg$down), length(imReg$downReg) - length(niimCompReg$down),
        "Downregulated Genes Shared by NODvsIRT and IRTvsMRT",
        labl = "NODvsIRT",
        labc = "Both",
        labr = "IRTvsMRT"
)
dev.off()

for (i in niimCompReg$up) {
    plotSpliceReg(niLFC, "IRTvsMRT", paste("gene_biotype mRNA;", i))
    plotIsoform(gene = i, annotation = "./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf")
}

## NODvsIRT X NODvsMRT ----

png("./comp/NODvsIRTxNODvsMRTup.png")
plotVen(length(niReg$upReg) - length(ninmCompReg$up), length(ninmCompReg$up), length(nmReg$upReg) - length(ninmCompReg$up),
        "Upregulated Genes Shared by NODvsIRT and NODvsMRT",
        labl = "NODvsIRT",
        labc = "Both",
        labr = "NODvsMRT"
)
dev.off()

png("./comp/NODvsIRTxNODvsMRTdown.png")
plotVen(length(niReg$downReg) - length(ninmCompReg$down), length(ninmCompReg$down), length(nmReg$downReg) - length(ninmCompReg$down),
        "Downregulated Genes Shared by NODvsIRT and NODvsMRT",
        labl = "NODvsIRT",
        labc = "Both",
        labr = "NODvsMRT"
)
dev.off()

for (i in ninmCompReg$up) {
    plotSpliceReg(niLFC, "IRTvsMRT", paste("gene_biotype mRNA;", i))
    plotIsoform(gene = i, annotation = "./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf")
}

plotIsoform(gene = "MtrunA17_Chr8g0344721", annotation = "./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf")

## NODvsMRT X IRTvsMRT ----

png("./comp/NODvsMRTxIRTvsMRTup.png")
plotVen(length(nmReg$upReg) - length(nmimCompReg$up), length(nmimCompReg$up), length(imReg$upReg) - length(nmimCompReg$up),
        "Upregulated Genes Shared by NODvsMRT and IRTvsMRT",
        labl = "NODvsMRT",
        labc = "Both",
        labr = "IRTvsMRT"
)
dev.off()

png("./comp/NODvsMRTxIRTvsMRTdown.png")
plotVen(length(nmReg$downReg) - length(nmimCompReg$down), length(nmimCompReg$down), length(imReg$downReg) - length(nmimCompReg$down),
        "Downregulated Genes Shared by NODvsMRT and IRTvsMRT",
        labl = "NODvsMRT",
        labc = "Both",
        labr = "IRTvsMRT"
)
dev.off()

for (i in ninmCompReg$up) {
    plotSpliceReg(niLFC, "IRTvsMRT", paste("gene_biotype mRNA;", i))
    plotIsoform(gene = i, annotation = "./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf")
}


# Comp Props ----

NODprop <- read.csv("./isoform_proportions NOD")
IRTprop <- read.csv("./isoform_proportions IRT")
MRTprop <- read.csv("./isoform_proportions MRT")
allprop <- read.csv("./isoform_proportions.csv")

NODprop <- NODprop[,2:ncol(NODprop)]
IRTprop <- IRTprop[,2:ncol(IRTprop)]
MRTprop <- MRTprop[,2:ncol(MRTprop)]
allprop <- allprop[,2:ncol(allprop)]

NODprop$prop <- nan.to.zero(NODprop$prop)
IRTprop$prop <- nan.to.zero(IRTprop$prop)
MRTprop$prop <- nan.to.zero(MRTprop$prop)
allprop$prop <- nan.to.zero(allprop$prop)

propTable <- data.frame(TXNAME = NODprop$TXNAME, GENEID = NODprop$GENEID,
                        NOD = NODprop$prop,
                        IRT = IRTprop$prop,
                        MRT = MRTprop$prop,
                        ALL = allprop$prop
)

gene = "gene_biotype mRNA; MtrunA17_Chr2g0322801"

plotPropComp(propTable, gene)

plotPropCompSep(propTable, gene)

plotHeatmapIso(propTable[,], gene, c("blue", "pink"))
plotHeatmapIso(propTable, "gene_biotype mRNA; MtrunA17_Chr8g0344721", c("blue", "yellow"))
