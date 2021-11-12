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
qc_files <- list.files(params$procdata, ".*qc.txt", full.names=TRUE, 
    recursive=TRUE)
dt_list <- lapply(qc_files, FUN=fread)
qc_dt <- rbindlist(dt_list)
qc_dt$sample_name <- gsub("_.*$", "", basename(qc_files))

# 2 -- Output the QC results
fwrite(qc_dt, output$qc_csv)

# 3 -- Load the SummarizedExperiments and filter on QC criteria
cnv_list <- lapply(input$cnv_objects, FUN=qread)

# 4 -- Identify samples passing QC
qc_dt <- qc_dt[
    MAPD < params$mapd & `nd-waviness-sd` < params$ndwavinesssd &
        SNPQC > params$snpqc,
]
if (params$cellularity > 0) {
    qc_dt <- qc_dt[`TUSCAN-cellularity` > params$cellularity]
}
# Some samples may pass qc but failed to convege when selecting optimal gamma
# The results of optimal gamma selection are non-deterministic and depend on
#   the pipeline configuration
# As a result we need to interesect the actual samples with those passing qc
existing_samples <- unique(Reduce(c, lapply(cnv_list, colnames)))
keep_samples <- intersect(existing_samples, qc_dt$sample_name)

# 5 -- Subset SummarizedExperiments
.subset_by_class <- function(x, keep) { 
    if (is(x, "SummarizedExperiment")) {
        x[, keep]
    } else if(is(x, "GRangesList")) {
        x[keep] 
    } else { 
        stop("Unsupported object class!") 
    }
}
cnv_list <- Map(f=.subset_by_class, x=cnv_list, keep=list(keep_samples))
for (i in seq_along(cnv_list)) {
    metadata(cnv_list[[i]])$qc_parameters <- params[-1]
}

# 6 -- Write SummarizedExperment to disk
for (i in seq_along(cnv_list)) {
    qsave(cnv_list[[i]], output$cnv_objects[i])
}