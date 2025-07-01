library(igraph)
library(limma)

data <- read.csv("./NODvsIRT/NODvsIRT_splice.csv")
for(i in 1:nrow(data)) {
    if(data$name[i] == "") {
        name <- trimWhiteSpace(strsplit(data$GeneID[i], ";")[[1]][2])
        data$name[i] <- name
    } else {
        next
    }
}

presub <- rbind(data[1:26,], data[28:nrow(data),])
presub <- data.frame(presub$GeneID, presub$name, presub[,13:ncol(presub)])
for(i in 1:nrow(presub)) {
    for(j in 1:ncol(presub)) {
        if(is.na(presub[i,j])) {
            presub[i,j] <- 0
        }
    }
}

subset <- data.frame(matrix(NA, ncol = ncol(presub)))
colnames(subset) <- colnames(presub)
count <- 1
for (i in 1:nrow(presub)) {
    if (sum(presub[i,4:ncol(presub)]) != 0)  {
        subset[count,] <- presub[i,]
        count <- count + 1
    } else {
        next
    }
}
latd <- data.frame(presub.GeneID = "gene_biotype mRNA; MtrunA17_Chr1g0147631",
                   presub.name = "LATD/NIP",
                   description..legoo. = "NRT1/ PTR FAMILY; MtNIP/LATD; nitrate transporter / peptide transporter;nitrate/peptide transporter family 1.7",
                   NAD1 = -1, 
                   NIN = 1,
                   NF.YA1 = 1,
                   DASH = 0,
                   ERN1 = 0,
                   DMI3 = 0,
                   DME = 0,
                   NSP1 = 0,
                   CRE1 = 0,
                   EFD = 0
)
subset[nrow(subset) + 1,] <- c(latd, rep(0, ncol(subset) - ncol(latd)))


vert <- data.frame(n = 1:(length(subset[,4:ncol(subset)]) + nrow(subset)),
                   name = c(colnames(subset[,4:ncol(subset)]), subset$presub.GeneID),
                   shortname = c(colnames(subset[,4:ncol(subset)]), subset$presub.name)
)

test <- 0
for(i in presub$presub.name) {
    for (j in colnames(presub)) {
        if (i == j) {
            print(i)
            print(j)
            test <- test + 1
        }
    }
}
test

edges <- data.frame(from = NA, to = NA)
ecolors <- c()
genes <- colnames(subset[,4:ncol(subset)])
count <- 1
for(i in 1:nrow(subset)) {
    row <- subset[i,]
    geneid <- row$presub.GeneID
    to <- vert[vert$name == geneid,]
    for (j in 1:length(genes)) {
        con <- row[,colnames(row) == genes[j]]
        con_id <- genes[j]
        if (con == 1) {
            from <- vert[vert$name == con_id,]
            edges[count,] <- c(from$n, to$n)
            ecolors[count] <- "#ff0019"
            count <- count + 1
        } else if (con == -1) {
            from <- vert[vert$name == con_id,]
            edges[count,] <- c(from$n, to$n)
            ecolors[count] <- "#192bcf"
            count <- count + 1
        }
    }
}

g <- graph_from_edgelist(as.matrix(edges))

reallyshortnames <- c()
count <- 1
for (i in vert$shortname) {
    print(i)
    new <- strsplit(i, "_")[[1]][2]
    if (is.na(new)) {
        reallyshortnames[count] <- i
    } else {
        reallyshortnames[count] <- new
    }
    count <- count + 1
}
vert$reallyshortname <- reallyshortnames

lay <- layout_(g, nicely())
colors <- c(rep("red", length(genes)), rep("#4287f5", nrow(vert) - length(genes)))
vert$color <- colors
deg <- degree(g, v = V(g), mode = "in")
vertex_attr(g) <- list(
    names = vert$reallyshortname,
    color = vert$color,
    degree = deg
)
lab.sizes <- c(rep(1, 10), rep(0.7, nrow(vert) - 6))
deg.add <- c(rep(10, length(genes)), rep(6, nrow(vert) - length(genes)))

par(xpd = T)
plot(g,
     layout = lay,
     edge.width = 2.5,
     edge.arrow.size = 0.3,
     edge.color = ecolors,
     vertex.label = vert$reallyshortname,
     vertex.label.cex = lab.sizes,
     vertex.label.font = 2,
     vertex.size = deg.add,
     asp = -5
     )

legend("topright",
       title = "Edges",
       legend = c("Induced/Activates", "Respresses/Inhibits"),
       col = c("#068033", "#192bcf"),
       lty = 1,
       lwd = 2,
       cex = 0.7,
       inset = c(0.2,-0.1),
       box.lwd = 0
)

legend("bottomleft",
       title = "Nodes",
       legend = c("Nod Genes", "Spliced Genes"),
       col = c("red", "#4287f5"),
       pch = 16,
       pt.cex = 1.5,
       cex = 0.7,
       inset = c(0,-0.1),
       box.lwd = 0
)
