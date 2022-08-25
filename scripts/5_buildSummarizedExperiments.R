# 0.1 -- Load Dependencies
renv::activate()
library(EaCoN)
library(Biobase)
library(SummarizedExperiment)
library(data.table)
library(qs)

# 0.2 -- Parse snakemake parameters
input <- snakemake@input
params <- snakemake@params
output <- snakemake@output

# 1 -- Load input data
gr.cnv <- qread(input$gr_cnv, nthread=params$nthreads)

# 2 -- Build ExpressionSets
pairs_file <- fread(input$pairs_file)
setDF(pairs_file)
EaCoN:::buildPSetOut(gr.cnv, params$analysis_name, file.path(params$results),
    meta=pairs_file[, -c(1)])

# 3 -- Convert ExpressionSet objects to SummarizedExperiments

esetFiles <- list.files(params$results,
    paste0(params$analysis_name, '.*ESet..RData'), full.names=TRUE)

esetNames <- list()
for (file in esetFiles) esetNames[file] <- load(file)

esets <- lapply(esetNames, get)
names(esets) <- gsub('.*/|_ESet..RData', '', esetFiles)

SEs <- lapply(esets, FUN=as, 'SummarizedExperiment')

# 4 -- Save results to disk
for (i in seq_along(SEs)) {
    qsave(SEs[[i]],
        file=file.path(params$results, paste0(names(SEs)[i], '_SumExp.qs')))
}