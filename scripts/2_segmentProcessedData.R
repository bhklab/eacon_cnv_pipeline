# 0 -- Load dependencies
library(EaCoN)

# 0 -- Parse Snakemake arguments
input <- snakemake@input
params <- snakemake@params
nthreads <- snakemake@threads
output <- snakemake@output

# 1 -- L2R and BAF joint segment the preprocessed data from the previous rule
EaCoN:::Segment.ff.Batch(RDS.files=unlist(input), nthread=nthreads,
    segmenter=params$segmenter, smooth.k=params$smoothk,
    BAF.filter=params$BAF_filter, nrf=params$nrf, SER.pen=params$SER_pen,
    force=TRUE)