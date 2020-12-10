# 0 -- Load dependencies
renv::activate()

library(EaCoN)
library(qs)

# 0 -- Parse Snakemake arguments
input <- snakemake@input
params <- snakemake@params
nthreads <- snakemake@threads
output <- snakemake@output

# 1 -- Preprocess and  normalize the .CEL files
if (grepl('cytoscan', params$array_family, ignore.case=TRUE)) {
    results <- CS.Process.Batch(input$pairs_file, 
        nthread=nthreads, out.dir='procdata', force=TRUE, write.data=FALSE,
        return.data=TRUE)
} else if (grepl('oncoscan', params$array_family, ignore.case=TRUE)) {
    results <- OS.Process.Batch(input$pairs_file, 
        nthread=nthreads, out.dir=params$procdata, force=TRUE, write.data=FALSE,
        return.data=TRUE)
} else if (grepl('snf6', params$array_family, ignore.case=TRUE)) {
    results <- SNF6.Process.Batch(input$paris_file, out.dir=params$procdata,
        write.data=FALSE, return.data=TRUE)
} else {
    stop("Supported array families are cytoscan, oncoscan and snp6")
}

# 2 -- Save the file to disk with the proper names
for (i in seq_along(results)) {
    qsave(results[[i]], file=file.path(params$procdata, 
        paste0('1.', params$analysis_name, '_', params$sample_name[i], '_',
            params$array_family, '_preprocessed.qs')))
}