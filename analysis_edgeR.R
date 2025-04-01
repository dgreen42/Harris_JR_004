library(edgeR)

cts <- read.delim("./bambu_out/Harris_JR_RNA_004_counts_transcript.tsv")
rownames(cts) <- cts$TXNAME
cts <- cts[,2:ncol(cts)]
head(cts)

group <- rep(c("nod" ,"irt" ,"mrt"), 3)
dge <- DGEList(counts = cts, group = group)

#read chapter 3 on grouping
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

#clustering
logcpm <- cpm(dge, log = T)
plotMDS(logcpm)

#alt exp
diff <- diffSpliceDGE(qlfIRTvsNOD, contrast = IRTvsNOD, geneid = dge$genes$GENEID)
