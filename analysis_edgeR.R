library(Rsubread)
library(openxlsx)
library(edgeR)
library(bambu)
library(Repitools)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

genome <- file.path("./MtrunA17r5.0-20161119-ANR.genome.fasta")
annotation <- file.path("./MtrunA17r5.0-ANR-EGN-r1.9.gtf")
prepAnno <- prepareAnnotations(annotation)
samples <- list.files(path = "./aln", recursive = T, full.names = T)

#CHANGE NDR
analysis <- bambu(reads = samples, annotations = prepAnno, genome = genome, NDR = 0.7)
saveRDS(analysis, file = "./bambu_out_NDR_7/bambu_analysis.rds")

writeBambuOutput(analysis, path = "./bambu_out_NDR_7/", prefix = "Harris_JR_RNA_004_NDR_7")
writeToGTF(rowRanges(analysis), file = "./bambu_out_NDR_7/Harris_JR_RNA_004_NDR_7_extended_anntation.gtf")

rm(prepAnno)

png(filename = "./plots/bambu_heatmap.png")
plotBambu(analysis, type = "heatmap")
dev.off()

png(filename = "./plots/MtrunA17_Chr8g0390331_isoform_count.png")
plotBambu(analysis, type = "annotation", gene_id = "gene_biotype mRNA; MtrunA17_Chr8g0390331")
dev.off()

analysis_noquant <- bambu(reads = samples, annotations = prepAnno, genome = genome, quant = F)
writeBambuOutput(analysis_noquant, path = "./bambu_out", prefix = "Harris_JR_RNA_004_noquant")
writeToGTF(rowRanges(analysis_noquant), file = "Harris_JR_RNA_004_extended_anntation_noquant.gtf")

png(filename = "./plots/bambu_heatmap.png")
plotBambu(analysis, type = "heatmap")
dev.off()

png(filename = "./plots/MtrunA17_Chr8g0390331_isoform_count.png")
plotBambu(analysis, type = "annotation", gene_id = "gene_biotype mRNA; MtrunA17_Chr8g0390331")
dev.off()

#following example -------------------------------------------

analysis <- readRDS("./bambu_out_NDR_3/bambu_analysis.rds")
plotBambu(analysis, type = "annotation", gene_id = "gene_biotype mRNA; MtrunA17_Chr1g0200031")
cts <- assay(analysis)
head(cts)

cts <- read.delim("./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3counts_transcript.tsv")
rownames(cts) <- cts$TXNAME
cts <- cts[,2:ncol(cts)]
head(cts)

#anno2 <- read.delim("./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf", sep = "\t", header = F)
#anno2 <- anno2[anno2$V3 == "transcript", ]

group <- rep(c("nod" ,"irt" ,"mrt"), 3)

# dge <- DGEList(counts = cts, group = group, genes = anno2)
dge <- DGEList(counts = cts, group = group)
head(dge)

expDesign <- model.matrix(~0+group, data = dge$samples)
colnames(expDesign) <- levels(dge$samples$group)
expDesign

filterLowcts <- filterByExpr(dge, design = expDesign)
dge <- dge[filterLowcts, ]

dge <- normLibSizes(dge)

png("./plots/gene_expression_clustering.png")
par(xpd = NA, font = 2)
plotMDS(dge)
title("Clustering of Gene Expression Profiles")
dev.off()

dge <- estimateDisp(dge)
fit <- glmQLFit(dge, design = expDesign, robust = TRUE)
png("./plots/QLDispersion.png")
plotQLDisp(fit)
title("Quasi-Likelihood Dispersion Estimate")
dev.off()

NODvsIRT <- makeContrasts(nod-irt, levels = expDesign)
qlfNODvsIRT <- glmQLFTest(fit, contrast = NODvsIRT)
lfcNODvsIRT <- data.frame(qlfNODvsIRT$genes, qlfNODvsIRT$table)
write.csv(lfcNODvsIRT, "NODvsIRTlfc.csv")
topTags(qlfNODvsIRT)
dexp1 <- decideTests(qlfNODvsIRT, p.value = 0.05)
summary(dexp1)
alt_splice_1 <- diffSpliceDGE(fit, contrast = NODvsIRT, geneid = "GENEID")
top_splice <- topSpliceDGE(alt_splice_1)
top_splice
variants <- spliceVariants(dge, dge$genes, dge$common.dispersion)
plot(variants$table$PValue, variants$table$logCPM)
sig_var <- variants[variants$table$PValue < 0.01,]
sig_var2 <- variants[variants$table$PValue < 0.0000001,]
sig_var2

write.csv(top_splice, "./NODvsIRT_topSplice.csv")
write.csv(sig_var, "./NODvsIRT_variants1.csv")
write.csv(sig_var2, "./NODvsIRT_variants2.csv")

splice <- read.csv("NODvsIRT_splice.csv")
lfcNODvsIRT$BambuTx <- rownames(lfcNODvsIRT)
rownames(lfcNODvsIRT) <- 1:nrow(lfcNODvsIRT)

count <- 0
spliceLFC <- data.frame(matrix(ncol = ncol(lfcNODvsIRT)))
colnames(spliceLFC) <- colnames(lfcNODvsIRT)
for (i in splice$GeneID) {
    for (j in 1:nrow(lfcNODvsIRT)) {
        if (lfcNODvsIRT$GENEID[j] == i) {
            count = count + 1
            spliceLFC[count,] <- lfcNODvsIRT[j,]
        }
    }
}

spliceLFC
rownames(spliceLFC) <- spliceLFC$BambuTx
spliceLFC <- spliceLFC[,1:ncol(spliceLFC)-1]

write.csv(spliceLFC, "NODvsIRT_splice_DE.csv")

source("utils.R")

splice_1 <- splice[splice$Alt...0n..1y..2m. == 1,]

for (i in splice_1$GeneID) {
    plotSpliceReg(spliceLFC, "NODvsIRT", i)
    plotIsoform(gene = strsplit(i, ";")[[1]][2], annotation = "./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf")
}

png("./plots/MtrunA17_Chr1g0200031_LFC_NODvsIRT.png")
plotSpliceReg(spliceLFC, "NODvsIRT", "gene_biotype mRNA; MtrunA17_Chr1g0200031")
dev.off()
png("./plots/MtrunA17_Chr1g0200031_isoforms.png")
plotIsoform(gene = "MtrunA17_Chr1g0200031", annotation = "./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf", exon_marker = T)
dev.off()

dev.off()
jpeg("./plots/MtrunA17_Chr1g0200031_isoforms_LFC_NODvsIRT.jpeg", width = 800, height = 500, quality = 100)
layout(matrix(c(1,2), nrow = 2, ncol = 1))
iso <- plotIsoform(gene = "MtrunA17_Chr8g0340121", annotation = "./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf", exon_marker = T)
par(xpd = NA)
text(x = min(iso$xrange), y = 2, labels = "A", adj = c(13,-5))
plotSpliceReg(spliceLFC, "NODvsIRT", "gene_biotype mRNA; MtrunA17_Chr8g0340121")
par(xpd = NA)
text(x = 0, y = 6, labels = "B", adj = c(8, -5))

dev.off()

png("./plots/MtrunA17_Chr2g0322801_iso_reg.png",
    height = 600,
    width = 1000)
layout(matrix(c(1,2), nrow = 2, ncol = 1))
plotIsoform(gene = "MtrunA17_Chr2g0322801", annotation = "./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf", exon_marker = T)
plotSpliceReg(spliceLFC, "NODvsIRT", "gene_biotype mRNA; MtrunA17_Chr2g0322801")
dev.off()

single_cell_transcriptome <- read.csv("./Pereira_2024_single_cell_transcriptome-mmc6_CN1_toCN7.csv")

count <- 0
single_cell_xref <- c()
for (i in single_cell_transcriptome$gene_id) {
    for (j in splice$GeneID) {
        spGene <- trimWhiteSpace(gsub("_", "", strsplit(j, ";")[[1]][2]))
        if (i == spGene) {
            single_cell_xref <- append(single_cell_xref, i)
            count <- count + 1
        }
    }
}

total <- nrow(splice)
prop <- count/total
print(prop)
print(single_cell_xref)


IRTvsNOD <- makeContrasts(irt-nod, levels = expDesign)
qlfIRTvsNOD <- glmQLFTest(fit, contrast = IRTvsNOD)
lfcIRTvsNOD <- data.frame(qlfIRTvsNOD$genes, qlfIRTvsNOD$table)
write.csv(lfcIRTvsNOD, "IRTvsNODlfc.csv")
topTags(qlfIRTvsNOD)
dexp1 <- decideTests(qlfIRTvsNOD, p.value = 0.05)
summary(dexp1)



IRTvsMRT <- makeContrasts(mrt-irt, levels = expDesign)
qlfIRTvsMRT <- glmQLFTest(fit, contrast = IRTvsMRT)
alt_splice_2 <- diffSpliceDGE(fit, contrast = IRTvsMRT, geneid = "GENEID")
topSpliceDGE(alt_splice_2)

qlfIRTvsNOD$table$padj <- p.adjust(qlfIRTvsNOD$table$PValue, method = "BH")
tol10b.enriched <- qlfIRTvsNOD$table$logFC > 0 & qlfIRTvsNOD$table$padj < 0.01
tol10b.depleted <- qlfIRTvsNOD$table$logFC < 0 & qlfIRTvsNOD$table$padj < 0.01
plot(qlfIRTvsNOD$table$logFC, -10*log10(qlfIRTvsNOD$table$PValue),
     main = "Volcano Plot Inoculated Root Tip vs Nodule",
     xlab = "logFC",
     ylab = "-10*log10(P-val)"
)
points(qlfIRTvsNOD$table$logFC[tol10b.enriched], -10*log10(qlfIRTvsNOD$table$PValue[tol10b.enriched]),
       col = "red"
)
points(qlfIRTvsNOD$table$logFC[tol10b.depleted], -10*log10(qlfIRTvsNOD$table$PValue[tol10b.depleted]),
       col = "green"
)

DEXdata <- cts[tol10b.enriched | tol10b.depleted, ]
heatmap(as.matrix(DEXdata[,2:ncol(DEXdata)]))

