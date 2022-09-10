# 0.1 -- Load dependencies
stop("This is a work in progress! Please use 1_batchProcessRawdataFiles.R")

renv::activate()
library(EaCoN)


# 0.2 -- Parse Snakemake arguments
input <- snakemake@input
params <- snakemake@params
nthreads <- snakemake@threads
output <- snakemake@output
wildcards <- snakemake@wildcards


# 0.3 -- Load platform specific dependencies
if (grepl("snp|cytoscan|oncoscan", params$array_type, ignore.case=TRUE))
    library(affy.CN.norm.data)

if (grepl('snp', params$array_type)) {
        library(apt.snp6.1.20.0)
        library(rcnorm)
        library(GenomeWideSNP.6.na35.r1)
}
## FIXME:: Add conditional dependencies for other platforms!

switch(params$reference,
    'BSgenome.Hsapiens.UCSC.hg19'=library(BSgenome.Hsapiens.UCSC.hg19),
    'BSgenome.Hsapiens.UCSC.hg38'={
        if (grepl("snp", params$array_type, ignore.case=TRUE))
            stop("Must use BSgenome.Hsapiens.UCSC.hg19 for GenomeWide SNP6 arrays!")
        library(BSgenome.Hsapiens.UCSC.hg38)
    }
)

pairs_df <- read.table(params$pairs_file, sep="\t", header=TRUE,
    stringAsFactors=FALSE)
which_row <- grepl(wildcards$sample_name, pairs_df$cel_file)
sample_name <- pairs_df[which_row, "SampleName", drop=TRUE]
file_path <- paris_df[which_row, "cel_files"]

# 1 -- Preprocess and normalize the raw data
if (grepl('cytoscan', params$array_type, ignore.case=TRUE)) {
    result <- EaCoN:::CS.Process(file_path=file_path, sample_name=sample_name,
        out.dir=params$procdata, force=TRUE, cluter.type=params$cluster_type,
        return.data=TRUE)
} else if (grepl('oncoscan', params$array_type, ignore.case=TRUE)) {
    result <- EaCoN:::OS.Process(file_path=file_path, sample_name=sample_name,
        out.dir=params$procdata, force=TRUE, cluter.type=params$cluster_type,
        return.data=TRUE)
} else if (grepl('snp', params$array_type, ignore.case=TRUE)) {
    result <- EaCoN:::SNP6.Process(file_path=file_path, sample_name=sample_name,
        out.dir=params$procdata, force=TRUE, cluter.type=params$cluster_type,
        return.data=TRUE)
} else if (grepl('wes', params$array_type, ignore.case=TRUE)) {
    stop("WES has not been implemented in this pipeline yet, please see
        https://github.com/gustaveroussy/EaCoN for information on
        setting up your own analysis script.")
} else {
    stop("Supported assay families are wes, cytoscan, oncoscan and snp6")
}

# -- 2. Save the results to disk
writeRDS(result, file=output)
