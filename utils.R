
plotSpliceReg <- function(data, GeneID) {
    isoforms <- data[data$GENEID == GeneID,]
    barplot(isoforms$logFC, main = paste0("LFC of ", GeneID))
    for (i in 1:nrow(isoforms)) {
        axis(1, i, rownames(isoforms)[i], las = 1, cex.axis = 0.8)
    }
}