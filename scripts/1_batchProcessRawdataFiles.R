# 0 -- Load dependencies
renv::activate()

library(EaCoN)

# 0 -- Parse Snakemake arguments
input <- snakemake@input
params <- snakemake@params
nthreads <- snakemake@threads
output <- snakemake@output

# 1 -- 
if (grepl('cytoscan', params$array_family, ignore.case=TRUE)) {
    CS.Process.Batch(input$pairs_file, 
        nthread=nthreads, out.dir='procdata', force=TRUE)
} else if (grepl('oncoscan', params$array_family, ignore.case=TRUE)) {
    OS.Process.Batch(input$pairs_file, 
        nthread=nthreads, out.dir=params$procdata, force=TRUE)
} else if (grepl('snf6', params$array_family, ignore.care=TRUE)) {
    SNF6.Process.Batch(input$paris_file, out.dir=params$procdata)
} else {
    stop("Supported array families are cytoscan, oncoscan and snp6")
}
