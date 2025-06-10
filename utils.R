
plotSpliceReg <- function(data, set, GeneID) {
    par(mai = c(1.02,0.82,0.82,0.42))
    isoforms <- data[data$GENEID == GeneID,]
    barplot(isoforms$logFC,
            main = paste0("LFC of ", GeneID, "\n", set),
            xlab = "Transcript",
            ylab = "log fold change (LFC)"
            )
    if (nrow(isoforms) > 3) {
        font.size <- 0.6
    } else {
        font.size <- 0.8
    }
    for (i in 1:nrow(isoforms)) {
        axis(1, i, rownames(isoforms)[i], las = 1, cex.axis = font.size)
    }
}

plotIsoform <- function(gene, isoforms = NULL, annotation, exon_marker = F) {
    grepd <- system2("grep", args = paste(gene, annotation), stdout = T)
    rawFeatures <- strsplit(grepd, split = "\t")
    featureFrame <- data.frame(matrix(NA, ncol = length(rawFeatures[[1]]), nrow = length(rawFeatures)))
    colnames(featureFrame) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
    
    for (i in 1:length(rawFeatures)) {
        featureFrame[i,] <- rawFeatures[[i]]
    }
    featureFrame$attribute <- strsplit(featureFrame$attribute, split = ";")
    
    for (i in 1:nrow(featureFrame)) {
        if (featureFrame[i,]$feature == "transcript") {
            geneid <- paste(strsplit(featureFrame[i,]$attribute[[1]][1], split = "\"")[[1]][2],
                            gsub("\"", "", featureFrame[i,]$attribute[[1]][2]))
            transcriptid <- strsplit(featureFrame[i,]$attribute[[1]][3], split = "\"")[[1]][2]
            featureFrame$geneid[i] = geneid
            featureFrame$transcriptid[i] = transcriptid 
            featureFrame$exonid[i] = NA
        } else if (featureFrame[i,]$feature == "exon") {
            geneid <- paste(strsplit(featureFrame[i,]$attribute[[1]][1], split = "\"")[[1]][2],
                            gsub("\"", "", featureFrame[i,]$attribute[[1]][2]))
            transcriptid <- strsplit(featureFrame[i,]$attribute[[1]][3], split = "\"")[[1]][2]
            exonid <- strsplit(featureFrame[i,]$attribute[[1]][4], split = "\"")[[1]][2]
            featureFrame$geneid[i] = geneid
            featureFrame$transcriptid[i] = transcriptid 
            featureFrame$exonid[i] = exonid 
        }
    }
    
    transcripts <- unique(featureFrame$transcriptid)
    par(mai = c(1.02,2,0.82,0.42))
    xlimit =  c(as.integer(min(featureFrame$start)), as.integer(max(featureFrame$end)))
    plot(NA,
         xlim = xlimit,
         ylim = c(0,length(transcripts)),
         main = gene,
         xlab = "Location (bp)",
         ylab = NA,
         yaxt = "n",
         bty = "n",
    )
    axis(side = 2,
         at = 1:length(transcripts),
         labels = transcripts,
         las = 1,
         cex.axis = 0.8)
    count <- 0
    for (i in transcripts) {
        count <- count + 1
        transFrame <- featureFrame[featureFrame$transcriptid == i, ]
        for(j in 1:nrow(transFrame)) {
            row <- transFrame[j,]
            if (row$feature == "transcript") {
                lines(x = c(row$start, row$end), y = c(count,count), lwd = 1, lend = 1)
            } else if (row$feature == "exon") {
                lines(x = c(row$start, row$end), y = c(count,count), lwd = 10, lend = 1)
                if (exon_marker == T) {
                    abline(v = row$start, lty = 2)
                    abline(v = row$end, lty = 2)
                }
            }
        }
    }
    return(list(xrange = xlimit, transcripts = transcripts)) 
}

plotVen <- function(left, center, right, title, labl = NULL, labc = NULL, labr = NULL) {
    par(xpd = NA)
    plot(NA, xlim = c(0,100), ylim = c(0,100),
         xlab = NA, ylab = NA,
         xaxt = "n", yaxt = "n",
         bty = "n")
    points(35, 40, cex = 40, col = rgb(1,0,0,0.5), pch = 16)
    points(65, 40, cex = 40, col = rgb(0,0,1,0.5), pch = 16)
    text(20, 40, paste(left))
    text(50, 40, paste(center))
    text(80, 40, paste(right))
    title(paste(title))
    text(20, -30, paste(labl))
    text(50, -30, paste(labc))
    text(80, -30, paste(labr))
}

regTx <- function(DEX) {
    upReg <- sum(DEX == 1)
    downReg <- sum(DEX == -1)
    noSig <- sum(DEX == 0)
    
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

ramna <- function(x) {
    y <- NULL
    count <- 1
    for (i in x) {
        if (is.na(i)) {
            next
        } else {
            y[count] <- i
            count <- count + 1
        }
    }
    return(y)
}

compareReg <- function(set1, set2) {
    count <- 1
    reg <- list(
        up = c(),
        down = c(),
        summary = c()
    )
    for (i in set1$upReg) {
        for(j in set2$upReg) {
            if (i == j) {
                reg$up[count] <- i
                count <- count + 1
            }
        }
    }
    reg$summary[1] <- count
    
    for (i in set1$downReg) {
        for(j in set2$downReg) {
            if (i == j) {
                reg$down[count] <- i
                count <- count + 1
            }
        }
    }
    reg$summary[2] <- count
    
    for (i in set1$noSig) {
        for(j in set2$noSig) {
            if (i == j) {
                count <- count + 1
            }
        }
    }
    reg$summary[3] <- count
    
    return(reg)
}
