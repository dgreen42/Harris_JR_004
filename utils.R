# functions for downstream data analysis

#' Title
#'
#' @param dxds 
#' @param threshold 
#' @param filter 
#'
#' @returns
#' @export
#'
#' @examples
filter_results <- function(dxds = NULL, threshold = 2, filter) {
    if (threshold > 0 ) {
        dxds_df <- as.data.frame(dxds)
        upreg_log2 <- data.frame(matrix(ncol  = length(colnames(dxds_df))))
        colnames(upreg_log2) <- colnames(dxds_df)
        count = 0
        for (i in 1:297716) {
            if (!is.na(eval(parse(text = paste0("dxds_df[i,]$", filter))))) {
                if (eval(parse(text = paste0("dxds_df[i,]$", filter))) > threshold) {
                    count <- count + 1
                    upreg_log2[count,] <- dxds_df[i,]
                }
            }
        }
        return(upreg_log2)
    } else if (threshold < 0) {
        dxds_df <- as.data.frame(dxds)
        downreg_log2 <- data.frame(matrix(ncol  = length(colnames(dxds_df))))
        colnames(downreg_log2) <- colnames(dxds_df)
        count = 0
        for (i in 1:297716) {
            if (!is.na(eval(parse(text = paste0("dxds_df[i,]$", filter))))) {
                if (eval(parse(text = paste0("dxds_df[i,]$", filter))) < threshold) {
                    count <- count + 1
                    downreg_log2[count,] <- dxds_df[i,]
                }
            }
        }
        return(downreg_log2)
    } else if (threshold == 0) {
        print("Please enter a threshold that is greater or less than zero")
        return(NULL)
    }
}

#' Title
#'
#' @param topDiff 
#' @param topDiff2 
#' @param topDiff3 
#'
#' @returns
#' @export
#'
#' @examples
make_diff_df <- function(topDiff, topDiff2, topDiff3) {
    gene_names <- unique(c(topDiff$GENEID, topDiff2$GENEID, topDiff3$GENEID))
    diff_df <- data.frame(matrix(0, nrow = length(gene_names), ncol = 3))
    rownames(diff_df) <- gene_names
    colnames(diff_df) <- c("IRTvsNOD", "IRTvsMRT", "MRTvsNOD")
    
    for (i in 1:nrow(diff_df)) {
        for (j in topDiff$GENEID) {
            if (rownames(diff_df[i,]) == j) {
                diff_df$IRTvsNOD[i] = 1
            }
        }
    }
    
    for (i in 1:nrow(diff_df)) {
        for (j in topDiff2$GENEID) {
            if (rownames(diff_df[i,]) == j) {
                diff_df$IRTvsMRT[i] = 1
            }
        }
    }
    
    for (i in 1:nrow(diff_df)) {
        for (j in topDiff3$GENEID) {
            if (rownames(diff_df[i,]) == j) {
                diff_df$MRTvsNOD[i] = 1
            }
        }
    }
    return(diff_df)
}

diff_hm <- function(diff_df) {
    heatmap(data.matrix(diff_df),
            margins = c(5,5),
            cexCol = 0.8,
            cexRow = 0.8,
            col = colorRampPalette(brewer.pal(3, "Blues"))(2),
    )
    
    legend("topleft",
           legend = c("Absent", "Present"),
           fill = colorRampPalette(brewer.pal(3, "Blues"))(2),
    )
}

analyze_alt_expr <- function(fit, dge, name = "analysis") {
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
    return(list(topDiff, topDiff2, topDiff3))
}

prep_dge <- function(cts) {
    rownames(cts) <- cts$TXNAME
    cts <- cts[,2:ncol(cts)]
    head(cts)
    group <- rep(c("nod" ,"irt" ,"mrt"), 3)
    dge <- DGEList(counts = cts, group = group)
    return(dge)
}

get_fit <- function(dge, expDesign) {
    filterLowcts <- filterByExpr(dge, design = expDesign)
    dge <- dge[filterLowcts, ]
    
    dge <- normLibSizes(dge)
    dge <- estimateDisp(dge)
    fit <- glmQLFit(dge, design = expDesign)
    return(fit)
}