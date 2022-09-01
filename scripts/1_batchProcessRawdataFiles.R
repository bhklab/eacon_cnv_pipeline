# 0.1 -- Load dependencies
renv::activate()
library(EaCoN)


# 0.2 -- Parse Snakemake arguments
input <- snakemake@input
params <- snakemake@params
nthreads <- snakemake@threads
output <- snakemake@output


# 0.3 -- Load platform specific dependencies
if (grepl("snp|cytoscan|oncoscan", params$array_type, ignore.case=TRUE))
    library(affy.CN.norm.data)

if (grepl('snp', params$array_type)) {
        library(apt.snp6.1.20.0)
        library(rcnorm)
        library(GenomeWideSNP.6.na35.r1)
}
## FIXME:: Add conditional dependencies for other platforms

switch(params$reference,
    'BSgenome.Hsapiens.UCSC.hg19'=library(BSgenome.Hsapiens.UCSC.hg19),
    'BSgenome.Hsapiens.UCSC.hg38'={
        if (grepl("snp", params$array_type, ignore.case=TRUE))
            stop("Must use BSgenome.Hsapiens.UCSC.hg19 for GenomeWide SNP6 arrays!")
        library(BSgenome.Hsapiens.UCSC.hg38)
    }
)


# 1 -- Load or create the metadata file specifying CEL paths
if (file.exists(input$pairs_file)) {
    pairs_df <- read.table(input$pairs_file, sep="\t", header=TRUE,
        stringsAsFactors=FALSE)
} else {
    # find all CEL files in the raw data
    cel_file_paths <- list.files(params$rawdata, pattern="*.CEL$",
        recursive=TRUE, full.names=TRUE)
    pairs_df <- data.frame(
        cel_files=cel_file_paths,
        # assumes the second element in path is the sample name
        SampleName=vapply(cel_file_paths, FUN=function(x) strsplit(x)[[1]][2],
            FUN.VALUE=character(1))
    )
    if (!is.null(input$pairs_file) || input$paris_file == "") {
        # create path if it doesn't exist
        pairs_path <- dirname(input$pairs_file)
        if (!file.exists(pairs_path)) dir.create(pairs_path, recursive=TRUE)
        # write out a pairs file
        write.table(pairs_df, input$pairs_file)
    }
}


# 2 -- Format the paths in the pairs_file to match this project directory
#   structure from config.yaml
if (grepl('cytocscan|oncoscan', params$array_type, ignore.case=TRUE)) {
    stopifnot(c("ATChannelCel", "GCChannelCel", "SampleName") %in% colnames(pairs_df))
    # remove existing path, if there is one
    pairs_df$ATChannelCel <- gsub("^.*\\/", "", pairs_df$ATChannelCel)
    pairs_df$GCChannelCel <- gsub("^.*\\/", "", pairs_df$GCChannelCel)
    # create a new path relative to specified rawdata directory
    pairs_df$GCChannelCel <- file.path(getwd(), params$rawdata,
        pairs_df$GCChannelCel)
} else if (grepl('snp6', params$array_type, ignore.case=TRUE)) {
    stopifnot(all(c("cel_files", "SampleName") %in% colnames(pairs_df)))
    # pairs_df$cel_files <- gsub("^.*\\/", "", pairs_df$cel_files)
    # pairs_df$cel_files <- file.path(getwd(), params$rawdata, pairs_df$cel_files)
}
# output the file to temporary storage so it can be read by EaCoN
pairs_file <- file.path(tempdir(), "CEL_pairs_file.csv")
write.table(pairs_df, file=pairs_file, sep="\t")


# 3 -- Preprocess and normalize the raw data; does
if (grepl('cytoscan', params$array_type, ignore.case=TRUE)) {
    EaCoN:::CS.Process.Batch(pairs_file,
        nthread=nthreads, out.dir=params$procdata, force=TRUE,
        cluter.type=params$cluster_type)
} else if (grepl('oncoscan', params$array_type, ignore.case=TRUE)) {
    EaCoN:::OS.Process.Batch(pairs_file,
        nthread=nthreads, out.dir=params$procdata, force=TRUE,
        cluster.type=params$cluster_type)
} else if (grepl('snp', params$array_type, ignore.case=TRUE)) {
    EaCoN:::SNP6.Process.Batch(pairs_file, out.dir=params$procdata, force=TRUE,
        nthread=nthreads, cluster.type=params$cluster_type)
} else if (grepl('wes', params$array_type, ignore.case=TRUE)) {
    stop("WES has not been implemented in this pipeline yet, please see
        https://github.com/gustaveroussy/EaCoN for information on
        setting up your own analysis script.")
} else {
    stop("Supported assay families are wes, cytoscan, oncoscan and snp6")
}
