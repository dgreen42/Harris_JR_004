# These are utility functions that are not used anymore

regTx <- function(DEX) {
    # upReg <- sum(DEX == 1)
    # downReg <- sum(DEX == -1)
    # noSig <- sum(DEX == 0)
    
    upRegTx <- c()
    downRegTx <- c()
    noSigTx <- c()
    for (i in 1:nrow(DEX)) {
        if (DEX[[i]] == 1) {
            upRegTx[i] <- rownames(DEX)[i]
        } else if (DEX[[i]] == -1) {
            downRegTx[i] <- rownames(DEX)[i]
        } else if (DEX[[i]] == 0) {
            noSigTx[i] <- rownames(DEX)[i]
        }
    }
    
    list(
        upReg = ramna(upRegTx),
        downReg = ramna(downRegTx),
        noSig = ramna(noSigTx)
    )
}

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