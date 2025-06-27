# These are utility functions
library(RColorBrewer)

plotSpliceReg <- function(data, set, GeneID) {
    par(mai = c(1.02,0.82,0.82,0.42), xpd = F)
    isoforms <- data[data$GENEID == GeneID,]
    barplot(isoforms$logFC,
            main = paste0("LFC of ", GeneID, "\n", set),
            xlab = "Transcript",
            ylab = "log fold change (LFC)"
            )
    abline(h = 0, col = "red", lwd = 3)
    if (nrow(isoforms) > 3) {
        font.size <- 0.6
    } else {
        font.size <- 0.8
    }
    for (i in 1:nrow(isoforms)) {
        axis(1, i, rownames(isoforms)[i], las = 1, cex.axis = font.size)
    }
}

plotIsoform <- function(gene, isoforms = NULL, annotation, exon_marker = F, prop = NULL) {
    par(xpd = F)
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
    
    if (is.null(prop)) {
        transcripts <- unique(featureFrame$transcriptid)
        par(mai = c(1.02,2,1,0.42))
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
        mtext("Isoform",
              side = 3,
              padj = -1.2, 
              adj = -0.3,
              
        )
        return(list(xrange = xlimit, transcripts = transcripts)) 
        
    } else {
        transcripts <- unique(featureFrame$transcriptid)
        par(mai = c(1.02,2,1,1))
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
        propidx <- grep(gene, prop$GENEID)
        props <- c()
        plist <- data.frame(matrix(NA, ncol = ncol(prop)))
        colnames(plist) <- colnames(prop)
        count <- 1
        for(i in propidx) {
            plist[count,] <- prop[i,]
            count <- count + 1
        }
        print(plist)
        print(transcripts)
        count <- 1
        #fix ordering of props
        for (i in transcripts) {
           for (j in 1:nrow(plist)) {
               if (i == plist$TXNAME[j]) {
                   props[count] <- round(as.double(plist$prop[j]), digit = 2)
                   print(props[count])
                   count <- count + 1
               } else {
                   next
               }
           }
        }
        if (length(props) > 4) {
            if (length(props) > 6) {
                font.size = 0.6
            } else {
                font.size = 0.8
            }
        } else {
            font.size = 1
        }
        axis(side = 4,
             at = 1:length(propidx),
             labels = props,
             cex.axis = font.size
        )
        mtext("Isoform Proportion",
              side = 3,
              adj = 1.45,
              padj = -1.2
        )
        mtext("Isoform",
              side = 3,
              padj = -1.2, 
              adj = -0.3,
              
        )
        
        return(list(xrange = xlimit, transcripts = transcripts, props = plist)) 
    }
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

nan.to.zero <- function(x) {
    y <- NULL
    count <- 1
    for (i in x) {
        if (is.nan(i)) {
            y[count] <- 0
            count <- count + 1
        } else {
            y[count] <- i
            count <- count + 1
        }
    }
    return(y)

}

compareReg <- function(set1, set2) {
    summary <- list(
        upReg = NA,
        downReg = NA,
        noSig = NA
    )
    countup <- 1
    countdown <- 1
    countnosig <- 1
    reg <- list(
        up = c(),
        down = c(),
        summary = list()
    )
    for (i in set1$upReg) {
        for(j in set2$upReg) {
            if (i == j) {
                reg$up[countup] <- i
                countup <- countup + 1
            }
        }
    }
    summary$upReg <- countup
    
    for (i in set1$downReg) {
        for(j in set2$downReg) {
            if (i == j) {
                reg$down[countdown] <- i
                countdown <- countdown + 1
            }
        }
    }
    summary$downReg <- countdown
    
    for (i in set1$noSig) {
        for(j in set2$noSig) {
            if (i == j) {
                countnosig <- countnosig + 1
            }
        }
    }
    summary$noSig <- countnosig
    reg$summary <- summary
    
    return(reg)
}

getUniqueReg <- function(reg, compReg) {
    regUpUnique <- c()
    count <- 1
    for (i in reg) {
        test <- sum(compReg == i)
        if (test == 1) {
            next
        } else {
            regUpUnique[count] <- i
            count <- count + 1
        }
    }
    return(regUpUnique)
}

getGeneIDsRegUnique <- function(regUnique, master, geneLookup, n) {
    count <- 1
    for(i in regUnique) {
        for(j in 1:nrow(geneLookup)) {
            if (i == geneLookup$TXNAME[j]) {
                master[[n]][count,] <- geneLookup[j,]
                count <- count + 1
            } else {
                next
            }
        }
    }
    return(master)
}

regUnique <- function(reg1, reg2, compReg, master, geneLookup, setNames = NULL) {
    print("Starting Reg1")
    regUpUnique1 <- getUniqueReg(reg1, compReg)
    print("Reg1 complete")
    print("Starting Reg2")
    regUpUnique2 <- getUniqueReg(reg2, compReg)
    print("Reg2 complete")
    
    print("Starting Gene Lookup")
    master <- getGeneIDsRegUnique(regUpUnique1, master, geneLookup, 1)
    print("Reg1 complete")
    print("Starting Reg2")
    master <- getGeneIDsRegUnique(regUpUnique2, master, geneLookup, 2)
    print("Reg2 complete")
    print("Starting Comp")
    master <- getGeneIDsRegUnique(compReg, master, geneLookup, 3)
    print("Comp Complete")
    
    if (!is.null(setNames)) {
        names(master) <- setNames
    }
    
    return(master)
}

isoformProp2 <- function(counts) {
    len <- nrow(counts)
    props <- data.frame(TXNAME = NA, GENEID = NA, genetotal = NA, isototal = NA, prop = NA)
    genes <- unique(counts$GENEID)
    count <- 1
    for(gene in genes) {
        set <- counts[counts$GENEID == gene,]
        genetotal <- sum(set[,3:ncol(set)])
        for(i in 1:nrow(set)) {
            props[count,] <- c(set$TXNAME[i],
                               set$GENEID[i],
                               genetotal,
                               sum(set[i,3:ncol(set)]),
                               sum(set[i,3:ncol(set)])/genetotal
            )
            count <- count + 1
            if (count > len*0.25) {
                print("--- 25% complete ---")
            } else if (count > len*0.5) {
                print("--- 50% complete ---")
            } else if (count > len*0.75) {
                print("--- 75% complete ---")
            }
        }
    }
    return(props)
}

plotPropComp <- function(propTable, gene) {
    subset <- propTable[propTable$GENEID == gene,]
    complist <- c()
    count <- 1
    for(i in 1:nrow(subset)) {
        for (j in 3:ncol(subset)) {
            complist[count] <-subset[i,j]
            count <- count + 1
        }
    }
    div <- length(complist) / nrow(subset)
    cols <- c()
    for(i in 1:nrow(subset)) {
        cols <- append(cols, rep(colors()[i*4], div))
    }
    barplot(complist,
            width = 1,
            cex.names = 0.6,
            col = cols
    )
    namelocs <- c()
    for(i in 1:nrow(subset)) {
        namelocs[i] <- div * i 
    }
    axis(1,
         at = namelocs,
         labels = subset$TXNAME,
         cex.axis = 0.5
    )
}

plotPropCompSep <- function(propTable, gene) {
    subset <- propTable[propTable$GENEID == gene,]
    layout(matrix(1:4, ncol = 4, nrow = 1))
    par(mai = c(1.3, 0.5, 0.8, 0.5))
    title("Transcript Proportions Across Tissues")
    barplot(subset$NOD, main = "NOD", ylim = c(0,0.6),
            names.arg = subset$TXNAME,
            cex.names = 0.8,
            las = 2,
    )
    barplot(subset$IRT, main = "IRT", ylim = c(0,0.6),
            names.arg = subset$TXNAME,
            cex.names = 0.8,
            las = 2,
    )
    barplot(subset$MRT, main = "MRT", ylim = c(0,0.6),
            names.arg = subset$TXNAME,
            cex.names = 0.8,
            las = 2,
    )
    barplot(subset$all, main = "ALL", ylim = c(0,0.6),
            names.arg = subset$TXNAME,
            cex.names = 0.8,
            las = 2,
    )
}

plotPropHeatamp <- function(propTable, gene) {
    subset <- propTable[propTable$GENEID == gene,]
    propMat <- matrix(0, nrow = nrow(subset), ncol = ncol(subset) - 2)
    for (i in 1:nrow(subset)) {
        for (j in 3:ncol(subset)) {
            propMat[i,j - 2] <- as.numeric(subset[i,j])
        }
    }
    rownames(propMat) <- subset$TXNAME
    colnames(propMat) <- colnames(subset[,3:ncol(subset)])
    n <- 10
    nums <- round(seq(0, 1, length.out = n), digit = 1)
    cols <- hcl.colors(n, palette = "Reds")
    hm <- heatmap(
        propMat, Rowv = NA, Colv = NA, col = cols,
        cexRow = 0.8,
        cexCol = 0.8,
    )
    legend("bottomleft",
           legend = nums,
           col = c(cols),
           pch = 15,
           pt.cex = 3
    )
    return(hm)
}

plotHeatmapIso <- function(propTable, gene, color) {
    if (!is.data.frame(propTable)) {
        stop("must provide a data.frame")
    } else {
        subset <- propTable[propTable$GENEID == gene,]
        propMat <- matrix(0, nrow = nrow(subset), ncol = ncol(subset) - 2)
        for (i in 1:nrow(subset)) {
            for (j in 3:ncol(subset)) {
                propMat[i,j - 2] <- as.numeric(subset[i,j])
            }
        }
        rownames(propMat) <- subset$TXNAME
        colnames(propMat) <- colnames(subset[,3:ncol(subset)])
        
        divx <- seq(1, ncol(propMat))
        divy <- seq(1, nrow(propMat))
        xlimit = c(0, ncol(propMat))
        ylimit = c(0, nrow(propMat))
        par(xpd = T, mai = c(0.5, 2, 0.3, 1))
        plot(NA, xlim = xlimit, ylim = ylimit,  bty = "n", xaxt = "n", yaxt = "n",
             xlab = NA, ylab = NA)
        xaxisat <- 1:ncol(propMat) - 0.5
        yaxisat <- 1:nrow(propMat) - 0.5
        axis(1, at = xaxisat, labels = colnames(propMat), padj = -2, tick = F)
        axis(2, at = yaxisat, labels = rownames(propMat), tick = F, las = 2)
        propMat100 <- round(propMat * 100, digit = 0) + 1
        colors <- colorRampPalette(color)(max(propMat100))
        for (i in divx) {
            for (j in divy) {
                polygon(
                    x = c(i - 1, i, i, i - 1),
                    y = c(j - 1, j - 1, j, j),
                    col = colors[propMat100[j,i]]
                )
            }
        }
        nums <- round(seq(0, max(propMat), length.out = 10), digit = 2)
        lcolnums <- round(seq(1, length(colors), by = length(colors)/length(nums)), digit = 0)
        lcol <- c()
        count <- 1
        for(i in lcolnums) {
            lcol[count] <- colors[i]
            count <- count + 1
        }
        xmin <- max(xlimit) + 0.25
        xmax <- max(xlimit) + 1.25
        ymin <- max(ylimit) * 0.2
        ymax <- max(ylimit) * 0.92
        xleg <- c(xmin, xmax, xmax, xmin)
        yleg <- c(ymin, ymin, ymax, ymax)
        polygon(x = xleg,
                y = yleg
        )
        xadd <- 0.25
        yadd <- 0.1
        xsizer <- 0.2
        ysizer <- max(ylimit) *0.0667
        for(i in 1:length(nums)) {
            xcoord <- xmin + xadd
            ycoord <- ymin + yadd
            polygon(x = c(xcoord, xcoord + xsizer, xcoord + xsizer, xcoord),
                    y = c(ycoord, ycoord, ycoord + ysizer, ycoord + ysizer),
                    col = lcol[i]
            )
            text(x = xcoord + xsizer + 0.25,
                 y = ycoord + ysizer - 0.1,
                 labels = nums[i]
            )
            yadd <- yadd + ysizer
        }
        # legend(x = xleg,
        #        y = yleg,
        #        inset = c(-0.2,0),
        #        legend = nums,
        #        col = lcol,
        #        pch = 15,
        #        pt.cex = 3
        # )
    }
    
    returnl <- list(
        propMat = propMat,
        xlim = xlimit,
        ylim = ylimit,
        xleglim = xleg,
        yleglim = yleg
    )
    
    return(returnl)
}

