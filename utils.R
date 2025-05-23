library(bambu)

plotSpliceReg <- function(data, set, GeneID) {
    isoforms <- data[data$GENEID == GeneID,]
    barplot(isoforms$logFC, main = paste0("LFC of ", GeneID, "\n", set))
    if (nrow(isoforms) > 3) {
        font.size <- 0.6
    } else {
        font.size <- 0.8
    }
    for (i in 1:nrow(isoforms)) {
        axis(1, i, rownames(isoforms)[i], las = 1, cex.axis = font.size)
    }
}

plotIsoform <- function(gene, isoforms = NULL, annotation) {}

grepd <- system2("grep", args = "'MtrunA17_Chr8g0344721' ./bambu_out_NDR_3/Harris_JR_RNA_004_NDR_3_extended_anntation.gtf", stdout = T)
rawFeatures <- strsplit(grepd, split = "\t")
featureFrame <- data.frame(matrix(NA, ncol = length(rawFeatures[[1]]), nrow = length(rawFeatures)))
colnames(featureFrame) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
for (i in 1:length(rawFeatures)) {
    featureFrame[i,] <- rawFeatures[[i]]
}
featureFrame$attribute <- strsplit(featureFrame$attribute, split = ";")
for (i in 1:nrow(featureFrame)) {
    geneid <- paste(strsplit(featureFrame[i,]$attribute[[1]][1], split = "\"")[[1]][2],
                gsub("\"", "", featureFrame[i,]$attribute[[1]][2]))
    transcriptid <- featureFrame[i,]$attribute[[1]][3]
}


plot(NA, xlim = c(7870000,7879000), ylim = c(0,3))
lines(x = c(7870294, 7870803), y = c(1,1), lwd = 10, lend = 2)
