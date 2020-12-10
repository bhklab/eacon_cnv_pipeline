# 0 -- Load dependencies
renv::activate()

library(EaCoN)

# 0 -- Parse Snakemake arguments
input <- snakemake@input
params <- snakemake@params
nthreads <- snakemake@threads
output <- snakemake@output

# 1 -- L2R and BAF joint segment the preprocessed data from the previous rule
results <- Segment.ff.Batch(input, nthread=nthreads, write.data=FALSE, return.data=TRUE)

# 2 -- Save results to disk
for (i in seq_along(results)) {
    saveRDS(results[[i]], file=file.path(params$procdata, params$sample_name[i],
        paste0('2.', params$analysis_name, '.', params$sample_name[i], '.',
            params$array_family, '.segmented.rds'))
}