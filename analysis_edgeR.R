library(Rsubread)
library(openxlsx)
library(edgeR)
library(bambu)
library(Repitools)
library(ggplot2)
library(RColorBrewer)

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


group <- rep(c("nod" ,"irt" ,"mrt"), 3)

anno_path <- "./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf"
anno <- readFromGTF(anno_path)

#feature counts from Rsubread`
features <- featureCounts(samples,
                          annot.ext = anno_path,
                          isGTFAnnotationFile = T,
                          isLongRead = T,
                          strandSpecific = rep(1, length(samples)),
                          genome = genome
)

dge <- DGEList(counts = features$counts, group = group, genes = features$annotation)

anno_df <- as.data.frame(anno)
head(anno_df)

dge <- DGEList(counts = cts, group = group)
head(dge)

expDesign <- model.matrix(~0+group, data = dge$samples)
colnames(expDesign) <- levels(dge$samples$group)
expDesign

filterLowcts <- filterByExpr(dge, design = expDesign)
dge <- dge[filterLowcts, ]

dge <- normLibSizes(dge)

plotMDS(dge)

dge <- estimateDisp(dge)
fit <- glmQLFit(dge, design = expDesign, robust = TRUE)
plotQLDisp(fit)

IRTvsNOD <- makeContrasts(nod-irt, levels = expDesign)
qlfIRTvsNOD <- glmQLFTest(fit, contrast = IRTvsNOD)
topTags(qlfIRTvsNOD)
dexp1 <- decideTests(qlfIRTvsNOD, p.value = 0.05)
summary(dexp1)
alt_splice_1 <- diffSpliceDGE(fit, contrast = IRTvsNOD, geneid = "GENEID")
top_splice <- topSpliceDGE(alt_splice_1)
top_splice
variants <- spliceVariants(dge, dge$genes, dge$common.dispersion)
plot(variants$table$PValue, variants$table$logCPM)
sig_var <- variants[variants$table$PValue < 0.01,]
sig_var2 <- variants[variants$table$PValue < 0.0000001,]
sig_var2

write.csv(top_splice, "./IRTvsNOD_topSplice.csv")
write.csv(sig_var, "./IRTvsNOD_variants1.csv")
write.csv(sig_var2, "./IRTvsNOD_variants2.csv")


xlsxName <- "IRTvsNOD.xlsx"
wb <- createWorkbook()

addWorksheet(wb, "Regulation")
writeData(wb, sheet = 1, summary(dexp1))
addWorksheet(wb, "Top Splice")
writeData(wb, sheet = 2, top_splice)
addWorksheet(wb, "Variants")
writeData(wb, sheet = 3, variants)
addWorksheet(wb, "Significant Variants 0.01")
writeData(wb, sheet = 4, sig_var)
addWorksheet(wb, "Significant Variants 10^-7")
writeData(wb, sheet = 5, sig_var2)
saveWorkbook(wb, xlsxName)


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


IRTvsMRT <- makeContrasts(mrt-irt, levels = expDesign)
qlfIRTvsMRT <- glmQLFTest(fit, contrast = IRTvsMRT)
topTags(qlfIRTvsMRT)


MRTvsNOD <- makeContrasts(nod-mrt, levels = expDesign)
qlfMRTvsNOD <- glmQLFTest(fit, contrast = MRTvsNOD)
topTags(qlfMRTvsNOD)





#from saved counts -----------------------------------------------
source("./utils.R")
cts <- read.delim("./bambu_out/Harris_JR_RNA_004_counts_gene.tsv")
cts_transcript <- read.delim("./bambu_out/Harris_JR_RNA_004_counts_transcript.tsv")
cts_CPM <- read.delim("./bambu_out/Harris_JR_RNA_004_CPM_transcript.tsv")
cts_full_read <- read.delim("./bambu_out/Harris_JR_RNA_004_fullLengthCounts_transcript.tsv")
cts_unique <- read.delim("./bambu_out/Harris_JR_RNA_004_uniqueCounts_transcript.tsv")

rownames(cts) <- cts$TXNAME
cts <- cts[,2:ncol(cts)]
head(cts)
group <- rep(c("nod" ,"irt" ,"mrt"), 3)
dge <- DGEList(counts = cts, group = group)

expDesign <- model.matrix(~0+group, data = dge$samples)
colnames(expDesign) <- levels(dge$samples$group)
expDesign

filterLowcts <- filterByExpr(dge, design = expDesign)
dge <- dge[filterLowcts, ]

dge <- normLibSizes(dge)
dge <- estimateDisp(dge)
fit <- glmQLFit(dge, design = expDesign)

IRTvsNOD <- makeContrasts(nod-irt, levels = expDesign)
qlfIRTvsNOD <- glmQLFTest(fit, contrast = IRTvsNOD)
topTags(qlfIRTvsNOD)
IRTvsMRT <- makeContrasts(mrt-irt, levels = expDesign)
qlfIRTvsMRT <- glmQLFTest(fit, contrast = IRTvsMRT)
topTags(qlfIRTvsMRT)
MRTvsNOD <- makeContrasts(nod-mrt, levels = expDesign)
qlfMRTvsNOD <- glmQLFTest(fit, contrast = MRTvsNOD)
topTags(qlfMRTvsNOD)

diff <- diffSpliceDGE(fit, contrast = IRTvsNOD, geneid = dge$genes$GENEID)
topDiff <- topSpliceDGE(diff)
write.csv(as.data.frame(topSpliceDGE(diff)), file = paste("./edgeR_out/diff_splice_", name, "IRT_NOD.csv"))
plotSpliceDGE(diff)
diff2 <- diffSpliceDGE(fit, contrast = IRTvsMRT, geneid = dge$genes$GENEID)
topDiff2 <- topSpliceDGE(diff2)
write.csv(as.data.frame(topSpliceDGE(diff2)), file = paste("./edgeR_out/diff_splice_", name, "IRT_MRT.csv"))
plotSpliceDGE(diff2)
diff3 <- diffSpliceDGE(fit, contrast = MRTvsNOD, geneid = dge$genes$GENEID)
topDiff3 <- topSpliceDGE(diff3)
write.csv(as.data.frame(topSpliceDGE(diff3)), file = paste("./edgeR_out/diff_splice_", name, "MRT_NOD.csv"))
plotSpliceDGE(diff3)

#---------------------------------------

dge <- prep_dge(cts_gene)

#read chapter 3 on grouping
expDesign <- model.matrix(~0+group, data = dge$samples)
colnames(expDesign) <- levels(dge$samples$group)
expDesign

fit <- get_fit(dge, expDesign)

IRTvsNOD <- makeContrasts(nod-irt, levels = expDesign)
qlfIRTvsNOD <- glmQLFTest(fit, contrast = IRTvsNOD)
topTags(qlfIRTvsNOD)
IRTvsMRT <- makeContrasts(mrt-irt, levels = expDesign)
qlfIRTvsMRT <- glmQLFTest(fit, contrast = IRTvsMRT)
topTags(qlfIRTvsMRT)
MRTvsNOD <- makeContrasts(nod-mrt, levels = expDesign)
qlfMRTvsNOD <- glmQLFTest(fit, contrast = MRTvsNOD)
topTags(qlfMRTvsNOD)

#clustering
logcpm <- cpm(dge, log = T)
plotMDS(logcpm)

#alt exp
alt_exp_analysis <- analyze_alt_expr(fit, dge, "counts_gene")
alt_exp_analysis <- analyze_alt_expr(fit, dge, "counts_trancript")
alt_exp_analysis <- analyze_alt_expr(fit, dge, "CPM_trancript")
alt_exp_analysis <- analyze_alt_expr(fit, dge, "fullLengthCounts_trancript")
alt_exp_analysis <- analyze_alt_expr(fit, dge, "uniqueCounts_trancript")


topDiff <- alt_exp_analysis[[1]]
topDiff2 <- alt_exp_analysis[[2]]
topDiff3 <- alt_exp_analysis[[3]]

diff_df <- make_diff_df(topDiff, topDiff2, topDiff3)

diff_hm(diff_df)
