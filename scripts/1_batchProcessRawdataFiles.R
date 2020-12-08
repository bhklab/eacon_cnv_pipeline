# 0 -- Load dependencies
renv::activate()

library(EaCoN)

# 0 -- Parse Snakemake arguments
input <- snakemake@input
params <- snakemake@params
nthreads <- snakemake@threads
output <- snakemake@output

# 1 -- 
CS.Process.Batch(input$pairs_file, 
    nthread=nthreads, out.dir='procdata', force=TRUE)