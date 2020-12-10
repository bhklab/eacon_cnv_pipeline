# 0 -- Load dependencies
renv::activate()

library(EaCoN)

# 0 -- Parse Snakemake arguments
input <- snakemake@input
params <- snakemake@params
nthreads <- snakemake@threads
output <- snakemake@output

# 2 -- Generate .html qc reports for

# 1 -- Make copy number calls based on L2R segmentation results
results <- ASCN.ff.Batch(input, nthread=nthreads, write.data=FALSE, return.data=TRUE)


# 2 -- Save results to disk
for (i in seq_along(results)) {
    saveRDS(results[[i]], file=file.path(params$procdata, params$sample_name[i],
        paste0('2.', params$analysis_name, '.', params$sample_name[i], '.',
            params$array_family, '.segmented.rds'))
}