source("utils.R")

args <- commandArgs(trailingOnly = T)
file <- args[1]
name <- args[2]
print(file)
cts <- read.csv(file)
print(head(cts))
analysis <- isoformProp2(cts)
cts$prop <- nan.to.zero(cts$prop)


write.csv(analysis, paste("./isoform_proportions_", name, ".csv", sep = ""), row.names = F)

q()

NODprop <- read.csv("./isoform_proportions NOD")
IRTprop <- read.csv("./isoform_proportions IRT")
MRTprop <- read.csv("./isoform_proportions MRT")

NODprop <- NODprop[,2:ncol(NODprop)]
IRTprop <- IRTprop[,2:ncol(IRTprop)]
MRTprop <- MRTprop[,2:ncol(MRTprop)]

NODprop$prop <- nan.to.zero(NODprop$prop)
IRTprop$prop <- nan.to.zero(IRTprop$prop)
MRTprop$prop <- nan.to.zero(MRTprop$prop)
