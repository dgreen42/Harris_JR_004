library(bambu)

genome <- file.path("./MtrunA17r5.0-20161119-ANR.genome.fasta")
annotation <- file.path("./MtrunA17r5.0-ANR-EGN-r1.9.gtf")
prepAnno <- prepareAnnotations(annotation)
samples <- list.files(path = "./aln", recursive = T, full.names = T)

analysis <- bambu(reads = samples, annotations = prepAnno, genome = genome, NDR = 0.3)
saveRDS(analysis, file = "./bambu_out_NDR_7/bambu_analysis.rds")

writeBambuOutput(analysis, path = "./bambu_out_NDR_7/", prefix = "Harris_JR_RNA_004_NDR_7")
writeToGTF(rowRanges(analysis), file = "./bambu_out_NDR_7/Harris_JR_RNA_004_NDR_7_extended_anntation.gtf")

