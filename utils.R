# functions for downstream data analysis

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