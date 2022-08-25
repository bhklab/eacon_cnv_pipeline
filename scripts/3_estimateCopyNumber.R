# 0 -- Load dependencies
library(EaCoN)

# 0 -- Parse Snakemake arguments
input <- snakemake@input
params <- snakemake@params
nthreads <- snakemake@threads
output <- snakemake@output

# 2 -- Generate .html qc reports for

# 1 -- Make copy number calls based on L2R segmentation results
EaCoN:::ASCN.ff.Batch(unlist(input), nthread=nthreads, force=TRUE,
    gammaRange=unlist(params$gamma_range), cluster.type="FORK")
