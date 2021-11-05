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

# 1 -- Collect and parse QC data
