library(bambu)
library(DEXSeq)
library(ggplot2)

source("./utils.R")

genome <- file.path("./MtrunA17r5.0-20161119-ANR.genome.fasta")
annotation <- file.path("./MtrunA17r5.0-ANR-EGN-r1.9.gtf")
prepAnno <- prepareAnnotations(annotation)
samples <- list.files(path = "./aln", recursive = T, full.names = T)

analysis <- bambu(reads = samples, annotations = prepAnno, genome = genome)
writeBambuOutput(analysis, path = "./bambu_out", prefix = "Harris_JR_RNA_004_")
writeToGTF(rowRanges(analysis), file = "Harris_JR_RNA_004_extended_anntation.gtf")

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


# by tissue 
# running bambu
cts <- round(assay(analysis))

dxd <- DEXSeqDataSet(countData = cts,
                     sampleData = coldata,
                     design = ~sample + exon + tissue:exon,
                     featureID = rowData(analysis)$TXNAME,
                     groupID = rowData(analysis)$GENEID)
dxd = estimateSizeFactors(dxd)
dxd = estimateDispersions(dxd)
plotDispEsts(dxd)
dxd = testForDEU(dxd)
dxd = estimateExonFoldChanges(dxd)
dxr = DEXSeqResults(dxd)

dxds <- DEXSeq(dxd)
dxds_df <- as.data.frame(dxds)
write.csv(apply(dxds_df, 2, as.character), "./regulation/regulation_by_treatment.csv")

png("./plots/MA_plot_by_treatment_1.png")
y = 10
plotMA(dxds, ylim = c(-y, y))
title("./plots/Gene Regulation by Treament")
dev.off()
png("./plots/MA_plot_by_treatment_2.png")
y = 20
plotMA(dxds, ylim = c(-y, y))
title("./plots/Gene Regulation by Treament")
dev.off()

filter <- "log2fold_.mock_.inoc"
upreg_log2_2 <- filter_results(dxds, 2, "log2fold_.mock_.inoc")
upreg_log2_5 <- filter_results(dxds, 5, "log2fold_.mock_.inoc")
downreg_log2_2 <- filter_results(dxds, -2, "log2fold_.mock_.inoc")
downreg_log2_5 <- filter_results(dxds, -5, "log2fold_.mock_.inoc")

write.csv(apply(upreg_log2_5, 2, as.character), file = "./regulation/upregulated_5x_by_treatment.csv", row.names = F)
write.csv(apply(upreg_log2_2, 2, as.character), file = "./regulation/upregulated_2x_by_treatment.csv", row.names = F)
write.csv(apply(downreg_log2_2, 2, as.character), file = "./regulation/downregulated_2x_by_treatment.csv", row.names = F)
write.csv(apply(downreg_log2_5, 2, as.character), file = "./regulation/downregulated_5x_by_treatment.csv", row.names = F)

rm(dxd, dxds, dxds_df, downreg_log2_2, downreg_log2_5, upreg_log2_2, upreg_log2_5, cts)

# by treatment and tissue 
cts <- round(assay(analysis))
coldata <- read.csv("./sample_data.csv")
coldata

dxd2 <- DEXSeqDataSet(countData = cts,
                     sampleData = coldata,
                     design = ~sample + exon + tissue:exon + treatment:exon,
                     featureID = rowData(analysis)$TXNAME,
                     groupID = rowData(analysis)$GENEID)
sampleAnnotation(dxd2)
# normalization
dxd2 = estimateSizeFactors(dxd2)
# estiamte dispersion with formula
formulaFullModel = ~ sample + exon + tissue:exon + treatment:exon
formulaReducedModel = ~ sample + exon + tissue:exon
dxd2 = estimateDispersions(dxd2, formula = formulaFullModel)
plotDispEsts(dxd2)
dxd2 = testForDEU(dxd2,
                 reducedModel = formulaReducedModel,
                 fullModel = formulaFullModel)
dxd2 = estimateExonFoldChanges(dxd2, fitExpToVar = "tissue")
dxr2 <- DEXSeqResults(dxd2)
head(dxr2)

table( before = dxr$padj < 0.1, now = dxr2$padj < 0.1 )

y <- 20
plotMA(dxr2, ylim = c(-y,y))

dxr2$groupID[dxr2$groupID == "gene_biotypemRNA;MtrunA17_Chr1g0147631"]

plotDEXSeq(dxr2, "MtrunA17_Chr8g0390331", fitExpToVar = "tissue")

write.csv(apply(as.data.frame(DEXSeq(dxd2, fitExpToVar = "tissue")), 2, as.character), "./regulation/regulation_by_tissue_and_treatment.csv")


png("./plots/MA_plot_by_tissue_1.png")
y = 10
plotMA(dxds, ylim = c(-y, y))
title("./plots/Gene Regulation by Tissue")
dev.off()
png("./plots/MA_plot_by_tissue_2.png")
y = 20
plotMA(dxds, ylim = c(-y, y))
title("./plots/Gene Regulation by Tissue")
dev.off()

upreg_log2_2 <- filter_results(dxds, 2, "log2fold_.rt_.nod")
upreg_log2_5 <- filter_results(dxds, 5, "log2fold_.rt_.nod")
downreg_log2_2 <- filter_results(dxds, -2, "log2fold_.rt_.nod")
downreg_log2_5 <- filter_results(dxds, -5, "log2fold_.rt_.nod")

write.csv(apply(upreg_log2_5, 2, as.character), file = "./regulation/upregulated_5x_by_tissue.csv", row.names = F)
write.csv(apply(upreg_log2_2, 2, as.character), file = "./regulation/upregulated_2x_by_tissue.csv", row.names = F)
write.csv(apply(downreg_log2_2, 2, as.character), file = "./regulation/downregulated_2x_by_tissue.csv", row.names = F)
write.csv(apply(downreg_log2_5, 2, as.character), file = "./regulation/downregulated_5x_by_tissue.csv", row.names = F)

rm(dxd, dxds, dxds_df, downreg_log2_2, downreg_log2_5, upreg_log2_2, upreg_log2_5, cts)

# mean cts
nod <- data.frame(cts[,1], cts[,4], cts[,7])
irtj <- data.frame(cts[,2], cts[,5], cts[,8])
mrtj <- data.frame(cts[,3], cts[,6], cts[,9])

nod$avg <- apply(nod, 1, mean)
irtj$avg <- apply(irtj, 1, mean)
mrtj$avg <- apply(mrtj, 1, mean)
mean_cts <- data.frame(NOD = nod$avg, IRTJ = irtj$avg, MRTJ = mrtj$avg)
rownames(mean_cts) <- rownames(cts)

mean_sample_data <- read.csv("./mean_sample_data.csv")
mean_sample_data
mean_dxd <- DEXSeqDataSet(countData = round(mean_cts),
                          sampleData = mean_sample_data,
                          design = ~sample + exon + treatment + tissue:exon,
                          featureID = rowData(analysis)$TXNAME,
                          groupID = rowData(analysis)$GENEID)


sampleAnnotation(mean_dxd)
# normalization
mean_dxd = estimateSizeFactors(mean_dxd)
# estiamte dispersion with formula
formulaFullModel = ~ sample + exon + treatment + tissue:exon
formulaReducedModel = ~ sample + exon 
mean_dxd = estimateDispersions(mean_dxd, formula = formulaFullModel)
plotDispEsts(mean_dxd)
mean_dxd = testForDEU(mean_dxd,
                      reducedModel = formulaReducedModel,
                      fullModel = formulaFullModel)
mean_dxd = estimateExonFoldChanges(mean_dxd, fitExpToVar = "tissue")
dxr2 <- DEXSeqResults(mean_dxd)
head(mean_dxr)

table( before = dxr$padj < 0.1, now = dxr2$padj < 0.1 )

y <- 20
plotMA(dxr2, ylim = c(-y,y))



# no NOD samples
samples_noNOD <- list.files(path = "./aln_noNOD", recursive = T, full.names = T)

analysis_noNOD <- bambu(reads = samples_noNOD, annotations = prepAnno, genome = genome)
writeBambuOutput(analysis_noNOD, path = "./bambu_out", prefix = "Harris_JR_RNA_004_noNOD_")
writeToGTF(rowRanges(analysis_noNOD), file = "Harris_JR_RNA_004_extended_anntation_noNOD.gtf")

png(filename = "./plots/bambu_heatmap_noNOD.png")
plotBambu(analysis_noNOD, type = "heatmap")
dev.off()

png(filename = "./plots/MtrunA17_Chr8g0390331_isoform_count_noNOD.png")
plotBambu(analysis_noNOD, type = "annotation", gene_id = "gene_biotype mRNA; MtrunA17_Chr8g0390331")
dev.off()



# importing counts
cts <- read.delim("./bambu_out/Harris_JR_RNA_004_counts_gene.tsv")
head(cts)
rownames(cts) <- cts$GENEID
cts <- cts[,2:ncol(cts)]

t_cts <- read.delim("./bambu_out/Harris_JR_RNA_004_counts_transcript.tsv")
head(t_cts)
rownames(t_cts) <- t_cts$TXNAME
t_cts <- t_cts[,2:ncol(t_cts)]

coldata <- read.csv("./sample_data.csv")
rownames(coldata) <- coldata$sample
coldata <- coldata[,2:ncol(coldata)]
coldata

boxplot(cts)
boxplot(t_cts)
plot(NA, xlim = c(0, nrow(cts)), ylim = c(0, max(cts)))
points(x = 1:nrow(cts), y = cts$CINOD_aln, col = "blue")
points(x = 1:nrow(cts), y = cts$CIRTJ_aln, col = "red")

npf1b_cts <- cts[rownames(cts) == "gene_biotype mRNA; MtrunA17_Chr8g0390331", ]
npf1b_t_cts <- t_cts[t_cts$GENEID == "gene_biotype mRNA; MtrunA17_Chr8g0390331", ]

dxd <- DEXSeqDataSet(round(t_cts[,2:ncol(t_cts)]),
                     coldata,
                     design = ~treatment + exon + tissue:exon,
                     featureID = rownames(t_cts),
                     groupID = t_cts$GENEID)

dxd = estimateSizeFactors(dxd)
dxd = estimateDispersions(dxd)
plotDispEsts(dxd)
dxd = testForDEU(dxd, fullModel = ~treatment + exon + tissue:exon, reducedModel = ~treatment + exon)
dxd = estimateExonFoldChanges(dxd, fitExpToVar = "tissue")
dxr = DEXSeqResults(dxd)
nrow_padj <- table(dxr$padj < 0.1)[2]
dxr_df <- as.data.frame(dxr)
padj01_temp <- dxr_df[dxr_df$padj < 0.1, ]
padj01 <- data.frame(matrix(NA, nrow = nrow_padj, ncol = ncol(padj01_temp)))
count <- 0
for (i in 1:nrow(padj01_temp)) {
    if (!is.na(padj01_temp$padj[i])) {
        count <- count + 1
        padj01[count,] <- padj01_temp[i,]
    }
}
colnames(padj01) <- colnames(dxr_df)
padj01

plotMA(dxr, ylim = c(-10, 10))
write.csv(apply(dxr_df, 2, as.character), file = "./dxr.csv")

plotDEXSeq(dxr, "MtrunA17_Chr5g0427811", legend = T, displayTranscripts = T, fitExpToVar = "tissue")
plotDEXSeq(dxr, "MtrunA17_Chr8g0390331", legend = T, displayTranscripts = T, fitExpToVar = "tissue")
plotDEXSeq(dxr, "BambuTx1", legend = T, displayTranscripts = T, fitExpToVar = "tissue")


dxd <- DEXSeqDataSet(round(t_cts[,2:ncol(t_cts)]),
                     coldata,
                     design = ~rep + exon + tissue:exon,
                     featureID = rownames(t_cts),
                     groupID = t_cts$GENEID)

dxd = estimateSizeFactors(dxd)
dxd = estimateDispersions(dxd)
plotDispEsts(dxd)
dxd = testForDEU(dxd, fullModel = ~rep + exon + tissue:exon, reducedModel = ~rep + exon)
dxd = estimateExonFoldChanges(dxd, fitExpToVar = "tissue")
dxr = DEXSeqResults(dxd)

plotDEXSeq(dxr, "MtrunA17_Chr5g0427811", legend = T, displayTranscripts = T, fitExpToVar = "tissue")
plotDEXSeq(dxr, "MtrunA17_Chr8g0390331", legend = T, displayTranscripts = T, fitExpToVar = "tissue")
plotDEXSeq(dxr, "BambuTx1", legend = T, displayTranscripts = T, fitExpToVar = "tissue")

dxr[196998,]




library(pasilla)
data(pasillaDEXSeqDataSet, package="pasilla")
dxd <- estimateSizeFactors( dxd )
dxd <- estimateDispersions( dxd )
dxd <- testForDEU( dxd )
dxr <- DEXSeqResults( dxd )
