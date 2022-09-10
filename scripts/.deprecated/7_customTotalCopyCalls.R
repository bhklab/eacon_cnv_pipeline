# 0.1 -- Load dependencies
renv::activate()
library(GenomicRanges)
library(RaggedExperiment)
library(SummarizedExperiment)
library(matrixStats)
library(data.table)
library(qs)


# 0.2 -- Parse snakemake parameters
input <- snakemake@input
params <- snakemake@params
output <- snakemake@output


# 1 -- Read in Bioconductor objects
ragged_exp <- qread(input$ragged_exp)
bins <- metadata(ragged_exp)$annotated_genome_bins
rfun <- function(scores, query, ranges) {
    mean(scores, na.rm=TRUE)
}

tcn <- qreduceSummarizedExperiment(ragged_exp, bins, i="TCN", simplifyReduce=rfun)
rowData(tcn) <- mcols(bins)

tcn_mat <- assay(tcn, "TCN")

# 2 -- Assign custom TCN ranges
tcn_cutoffs <- params$tcn_cutoffs
tcn_cutoffs <- lapply(tcn_cutoffs, function(x) x[order(names(x))])



# 3 -- Write new objects to disk
