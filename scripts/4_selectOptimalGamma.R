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

gamma_df <- rbindlist(
    setNames(
        lapply(gammaFiles, FUN=fread),
        gsub(".gammaEval.txt", "", basename(gammaFiles))
    ),
    idcol="sample_name"
)

best_fits <- gamma_df[, .SD[which.max(GoF)], by="sample_name"]

best_fit_files <- Map(
    function(x, y)
        grep(pattern=y, list.files(x, recursive=TRUE, full.names=TRUE), value=TRUE),
    x=file.path(params$out_dir, best_fits$sample_name),
    y=paste0(".*gamma", sprintf("%.2f", best_fits$gamma), "/.*RDS$")
)

best_fit_data <- BiocParallel::bplapply(best_fit_files, FUN=readRDS)
