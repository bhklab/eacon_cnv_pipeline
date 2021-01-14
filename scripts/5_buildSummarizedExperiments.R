# 0.1 -- Load Dependencies
renv::activate()

library(EaCoN)
library(Biobase)
library(SummarizedExperiment)
library(qs)

# 0.2 -- Parse snakemake parameters
input <- snakemake@input
params <- snakemake@params
output <- snakemake@output

# 3 -- Convert ExpressionSet objects to SummarizedExperiments

esetFiles <- list.files(params$results, paste0(params$analysis_name, '.*ESet..RData'), full.names=TRUE)

esetNames <- list()
for (file in esetFiles) esetNames[file] <- load(file)

esets <- lapply(esetNames, get)
names(esets) <- gsub('.*/|_ESet..RData', '', esetFiles)

SEs <- lapply(esets, FUN=as, 'SummarizedExperiment')

# 4 -- Save results to disk
for (i in seq_along(SEs)) {
    qsave(SEs[[i]], file=file.path(params$results, paste0(names(SEs)[i], '_SumExp.qs')))
}