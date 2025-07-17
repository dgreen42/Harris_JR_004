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