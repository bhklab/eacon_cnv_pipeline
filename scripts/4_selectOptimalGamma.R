# 0.1 Load dependencies
renv::activate()
library(EaCoN)
library(data.table)
library(qs)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(S4Vectors)
library(org.Hs.eg.db)

# -- 0.2 Parse snakemake arguments
input <- snakemake@input
params <- snakemake@params
output <- snakemake@output

# -- 1. Read in gamma files

gammaFiles <- list.files(params$out_dir, '.*gammaEval.*txt', recursive=TRUE,
    full.names=TRUE)
print(gammaFiles)

gamma_df <- rbindlist(lapply(gammaFiles, FUN=fread))
