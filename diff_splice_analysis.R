library(edgeR)

# Import Counts ----
cts <- read.delim("./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3counts_transcript.tsv")
rownames(cts) <- cts$TXNAME
cts <- cts[,2:ncol(cts)]
head(cts)

#anno2 <- read.delim("./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf", sep = "\t", header = F)
#anno2 <- anno2[anno2$V3 == "transcript", ]

# Model ----
group <- rep(c("nod" ,"irt" ,"mrt"), 3)
dge <- DGEList(counts = cts, group = group)
head(dge)
expDesign <- model.matrix(~0+group, data = dge$samples)
colnames(expDesign) <- levels(dge$samples$group)
expDesign
filterLowcts <- filterByExpr(dge, design = expDesign)
dge <- dge[filterLowcts, ]
dge <- normLibSizes(dge)
dge <- estimateDisp(dge)
fit <- glmQLFit(dge, design = expDesign, robust = TRUE)
saveRDS(fit, "./edgeR_glm_fit.rds")

# Differnetial Analysis ----
## NODvsIRT ----
basePath <- "./NODvsIRT/"
NODvsIRT <- makeContrasts(nod-irt, levels = expDesign)
qlfNODvsIRT <- glmQLFTest(fit, contrast = NODvsIRT)
lfcNODvsIRT <- data.frame(qlfNODvsIRT$genes, qlfNODvsIRT$table)
write.csv(lfcNODvsIRT, paste(basePath, "NODvsIRTlfc.csv", sep = ""))
topTags(qlfNODvsIRT)
dexp1 <- decideTests(qlfNODvsIRT, p.value = 0.05)
dexp1Sum <- summary(dexp1)

alt_splice_1 <- diffSpliceDGE(fit, contrast = NODvsIRT, geneid = "GENEID")
top_splice <- topSpliceDGE(alt_splice_1)
top_splice
variants <- spliceVariants(dge, dge$genes, dge$common.dispersion)
plot(variants$table$PValue, variants$table$logCPM)
sig_var <- variants[variants$table$PValue < 0.01,]
sig_var2 <- variants[variants$table$PValue < 0.0000001,]
sig_var2


png(paste(basePath, "gene_expression_clustering.png", sep = ""))
par(xpd = NA, font = 2)
plotMDS(dge)
title("Clustering of Gene Expression Profiles")
dev.off()

png(paste(basePath, "QLDispersion.png", sep = ""))
plotQLDisp(fit)
title("Quasi-Likelihood Dispersion Estimate")
dev.off()

write.csv(top_splice, paste(basePath, "NODvsIRT_topSplice.csv", sep = ""))
write.csv(sig_var, paste(basePath, "NODvsIRT_variants1.csv", sep = ""))
write.csv(sig_var2, paste(basePath, "NODvsIRT_variants2.csv", sep = ""))

## NODvsMRT ----
basePath <- "./NODvsMRT/"
NODvsMRT <- makeContrasts(nod-mrt, levels = expDesign)
qlfNODvsMRT <- glmQLFTest(fit, contrast = NODvsMRT)
lfcNODvsMRT <- data.frame(qlfNODvsMRT$genes, qlfNODvsMRT$table)
write.csv(lfcNODvsMRT, paste(basePath, "NODvsMRTlfc.csv", sep = ""))
topTags(qlfNODvsMRT)
dexp1 <- decideTests(qlfNODvsMRT, p.value = 0.05)
dexp1Sum <- summary(dexp1)

alt_splice_1 <- diffSpliceDGE(fit, contrast = NODvsMRT, geneid = "GENEID")
top_splice <- topSpliceDGE(alt_splice_1)
top_splice
variants <- spliceVariants(dge, dge$genes, dge$common.dispersion)
plot(variants$table$PValue, variants$table$logCPM)
sig_var <- variants[variants$table$PValue < 0.01,]
sig_var2 <- variants[variants$table$PValue < 0.0000001,]
sig_var2


png(paste(basePath, "gene_expression_clustering.png", sep = ""))
par(xpd = NA, font = 2)
plotMDS(dge)
title("Clustering of Gene Expression Profiles")
dev.off()

png(paste(basePath, "QLDispersion.png", sep = ""))
plotQLDisp(fit)
title("Quasi-Likelihood Dispersion Estimate")
dev.off()

write.csv(top_splice, paste(basePath, "NODvsMRT_topSplice.csv", sep = ""))
write.csv(sig_var, paste(basePath, "NODvsMRT_variants1.csv", sep = ""))
write.csv(sig_var2, paste(basePath, "NODvsMRT_variants2.csv", sep = ""))

## IRTvsMRT ----
basePath <- "./IRTvsMRT/"
IRTvsMRT <- makeContrasts(irt-mrt, levels = expDesign)
qlfIRTvsMRT <- glmQLFTest(fit, contrast = IRTvsMRT)
lfcIRTvsMRT <- data.frame(qlfIRTvsMRT$genes, qlfIRTvsMRT$table)
write.csv(lfcIRTvsMRT, paste(basePath, "IRTvsMRTlfc.csv", sep = ""))
topTags(qlfIRTvsMRT)
dexp1 <- decideTests(qlfIRTvsMRT, p.value = 0.05)
dexp1Sum <- summary(dexp1)

alt_splice_1 <- diffSpliceDGE(fit, contrast = IRTvsMRT, geneid = "GENEID")
top_splice <- topSpliceDGE(alt_splice_1)
top_splice
variants <- spliceVariants(dge, dge$genes, dge$common.dispersion)
plot(variants$table$PValue, variants$table$logCPM)
sig_var <- variants[variants$table$PValue < 0.01,]
sig_var2 <- variants[variants$table$PValue < 0.0000001,]
sig_var2


png(paste(basePath, "gene_expression_clustering.png", sep = ""))
par(xpd = NA, font = 2)
plotMDS(dge)
title("Clustering of Gene Expression Profiles")
dev.off()

png(paste(basePath, "QLDispersion.png", sep = ""))
plotQLDisp(fit)
title("Quasi-Likelihood Dispersion Estimate")
dev.off()

write.csv(top_splice, paste(basePath, "IRTvsMRT_topSplice.csv", sep = ""))
write.csv(sig_var, paste(basePath, "IRTvsMRT_variants1.csv", sep = ""))
write.csv(sig_var2, paste(basePath, "IRTvsMRT_variants2.csv", sep = ""))