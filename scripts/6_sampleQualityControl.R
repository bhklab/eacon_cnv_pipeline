# 0.1 -- Load Dependencies
renv::activate()
library(EaCoN)
library(Biobase)
library(RaggedExperiment)
library(data.table)
library(qs)

# 0.2 -- Parse snakemake parameters
input <- snakemake@input
params <- snakemake@params
output <- snakemake@output

# 1 -- Collect and parse QC data
qc_files <- list.files(params$procdata, ".*qc.txt", full.names=TRUE,
    recursive=TRUE)
dt_list <- lapply(qc_files, FUN=fread)
qc_dt <- rbindlist(dt_list)
qc_dt$sample_name <- gsub("_.*$", "", basename(qc_files))

# 2 -- Output the QC results
fwrite(qc_dt, output$qc_csv)

# 3 -- Load the RaggedExperiment and filter on QC criteria
ragged_exp <- qread(input$ragged_exp)

# 4 -- Identify samples passing QC
qc_dt <- qc_dt[
    MAPD <= params$mapd & `waviness-sd` <= params$ndwavinesssd &
        SNPQC >= params$snpqc,
]
if (params$cellularity > 0) {
    if ("TUSCAN-cellularity" %in% colnames(qcdt))
        qc_dt <- qc_dt[`TUSCAN-cellularity` >= params$cellularity]
}
keep_samples <- intersect(colnames(ragged_exp), qc_dt$sample_name)
print(keep_samples)

# 5 -- Subset RaggedExperiment
ragged_exp <- ragged_exp[, keep_samples]
metadata(ragged_exp)$qc_parameters <- params[-c(1, 2)]

# 6 -- Write RaggedExperiment to disk
qsave(ragged_exp, file=output$ragged_exp, nthreads=params$nthreads)
