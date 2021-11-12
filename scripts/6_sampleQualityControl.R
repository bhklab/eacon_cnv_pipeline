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
qc_dt$sample_name <- gsub("_.*$", "", basename(qc_files))

# 2 -- Output the QC results
fwrite(qc_dt, output$qc_csv)

# 3 -- Load the SummarizedExperiments and filter on QC criteria
SE_list <- lapply(input$summarized_experiments, FUN=qread)

# 4 -- Identify samples passing QC
qc_dt <- qc_dt[
    MAPD < params$mapd & `nd-waviness-sd` < params$ndwavinesssd &
        SNPQC > params$snpqc,
]
if (params$cellurity > 0) {
    qc_dt <- qc_dt[`TUSCAN-cellularity` > params$cellularity]
}
keep_samples <- qc_dt$sample_name

# 5 -- Subset SummarizedExperiments
SE_list <- Map(f=subset, x=SE_list, subset=list(keep_samples))