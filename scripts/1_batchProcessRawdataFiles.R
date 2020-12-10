# 0 -- Load dependencies
renv::activate()

library(EaCoN)

# 0 -- Parse Snakemake arguments
input <- snakemake@input
params <- snakemake@params
nthreads <- snakemake@threads
output <- snakemake@output

# 1 -- Preprocess and normalize the raw data; does
if (grepl('cytoscan', params$array_type, ignore.case=TRUE)) {
    CS.Process.Batch(input$pairs_file, 
        nthread=nthreads, out.dir=params$procdata, force=TRUE)
} else if (grepl('oncoscan', params$array_type, ignore.case=TRUE)) {
    OS.Process.Batch(input$pairs_file, 
        nthread=nthreads, out.dir=params$procdata, force=TRUE)
} else if (grepl('snf6', params$array_type, ignore.care=TRUE)) {
    SNF6.Process.Batch(input$pairs_file, out.dir=params$procdata)
} else if (grepl('wes', params$array_type)) {
    stop("WES has not been implemented in this pipeline yet, please see 
        https://github.com/gustaveroussy/EaCoN for information on
        setting up your own analysis script.")
} else {
    stop("Supported array families are wes, cytoscan, oncoscan and snp6")
}
