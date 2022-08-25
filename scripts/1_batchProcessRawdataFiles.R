# 0.1 -- Load dependencies
library(EaCoN)


# 0.2 -- Parse Snakemake arguments
input <- snakemake@input
params <- snakemake@params
nthreads <- snakemake@threads
output <- snakemake@output


# 1 -- Format the paths in the pairs_file to match this project directory
#   structure from config.yaml
pairs_df <- read.table(input$pairs_file, sep="\t")
# remove existing path, if there is one
pairs_df$ATChannelCel <- gsub("^.*\\/", "", pairs_df$ATChannelCel)
pairs_df$GCChannelCel <- gsub("^.*\\/", "", pairs_df$GCChannelCel)
# create a new path relative to specified rawdata directory
pairs_df$ATChannelCel <- file.path(getwd(), params$rawdata,
    pairs_df$ATChannelCel)
pairs_df$GCChannelCel <- file.path(getwd(), params$rawdata,
    pairs_df$GCChannelCel)
# output the file to temporary storage so it can be read by EaCoN
pairs_file <- file.path(tempdir(), "CEL_pairs_file.csv")
write.table(pairs_df, file=pairs_file, sep="\t")


# 2 -- Preprocess and normalize the raw data; does
if (grepl('cytoscan', params$array_type, ignore.case=TRUE)) {
    CS.Process.Batch(pairs_file,
        nthread=nthreads, out.dir=params$procdata, force=TRUE,
        cluter.type="FORK")
} else if (grepl('oncoscan', params$array_type, ignore.case=TRUE)) {
    OS.Process.Batch(pairs_file,
        nthread=nthreads, out.dir=params$procdata, force=TRUE,
        cluster.type="FORK")
} else if (grepl('snf6', params$array_type, ignore.care=TRUE)) {
    SNF6.Process.Batch(pairs_file, out.dir=params$procdata, cluster.type="FORK")
} else if (grepl('wes', params$array_type)) {
    stop("WES has not been implemented in this pipeline yet, please see
        https://github.com/gustaveroussy/EaCoN for information on
        setting up your own analysis script.")
} else {
    stop("Supported assay families are wes, cytoscan, oncoscan and snp6")
}
