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
qc_files <- list.files(input$procdata, "STT.*qc.txt", full.names=TRUE, 
    recursive=TRUE)
dt_list <- lapply(qc_files, FUN=fread)
qc_dt <- rbindlist(dt_list)

# 2 -- Output the QC results
fwrite(qc_dt, output$qc_csv)