# These are utility functions that are not used anymore

isoformProp <- function(counts, subset = NULL) {
    current <- NULL
    after <- NULL
    gene <- NULL
    geneCount <- 1 
    set <- data.frame(TXNAME = NA, GENEID = NA, total = NA) 
    geneSum <- NULL
    for(i in 1:nrow(counts)) {
        if (is.null(gene)) {
            gene <- counts$GENEID[i]
        } 
        print(gene)
        if (i == nrow(counts)) {
            break
        } else {
            current <- counts[i,]
            after <- counts[i+1,]
            if (after$GENEID == current$GENEID) {
                # hard coded for first tow columns as transcript and gene info, the rest are samples
                set[geneCount,] <- c(current$TXNAME, current$GENEID, sum(current[3:ncol(current)]))
            } else {
                set[geneCount,] <- c(current$TXNAME, current$GENEID, sum(current[3:ncol(current)]))
                gene <- NULL
                geneCount <- geneCount + 1
            }
        }
    }
    return(set)
}