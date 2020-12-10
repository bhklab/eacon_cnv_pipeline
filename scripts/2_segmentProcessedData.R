# 0 -- Load dependencies
renv::activate()

library(EaCoN)

# 0 -- Parse Snakemake arguments
input <- snakemake@input
params <- snakemake@params
nthreads <- snakemake@threads
output <- snakemake@output

# 1 -- L2R and BAF joint segment the preprocessed data from the previous rule
Segment.ff.Batch(unlist(input), nthread=nthreads, force=TRUE)