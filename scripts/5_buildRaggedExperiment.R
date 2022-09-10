# 0.1 -- Load Dependencies
renv::activate()
library(EaCoN)
library(data.table)
library(qs)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(BiocParallel)
library(RaggedExperiment)

# -- 0.2 Parse snakemake parameters
input <- snakemake@input
params <- snakemake@params
output <- snakemake@output

# -- 1. Load the optimal gamma for each sample
best_fits <- fread(input)

# -- 2. Find the .RDS files associated with the best fit
best_fit_files <- Map(
    function(x, y)
        grep(pattern=y, list.files(x, recursive=TRUE, full.names=TRUE), value=TRUE),
    x=file.path(params$out_dir, best_fits$sample_name),
    y=paste0(".*gamma", sprintf("%.2f", best_fits$gamma), "/.*RDS$")
)
l2r_files <- Map(
    function(x, y)
        grep(pattern=y, list.files(x, recursive=TRUE, full.names=TRUE), value=TRUE),
    x=file.path(params$out_dir, best_fits$sample_name),
    y=".*L2R/.*RDS$"
)

# -- 3. Load the best fit ASCN and L2R data and build GRanges objects
ascn_data_list <- BiocParallel::bplapply(best_fit_files, FUN=readRDS)
l2r_data_list <- BiocParallel::bplapply(l2r_files, FUN=readRDS)

gr_list <- GRangesList(
    Map(f=buildGRangesFromASCNAndL2R, ascn_data_list, l2r_data_list)
)
# removing directory paths
names(gr_list) <- basename(names(gr_list))

# -- 4. Construct a RaggedExperiment object
ragged_exp <- as(gr_list, "RaggedExperiment")

# include annotated bins to summarize the RaggedExperiment with
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genome_bins <- binReferenceGenome()
annotated_bins <- annotateGRangesWithTxDB(genome_bins, txdb=txdb)

# include annotated genes to summarized RaggedExperiment with
gene_granges <- genes(txdb)
annotated_genes <- annotateGRangesWithTxDB(gene_granges, txdb=txdb)

metadata(ragged_exp) <- list(
    annotated_genome_bins=annotated_bins,
    annotated_genes=annotated_genes,
    simpifyReduce=function(scores, ranges, qranges) {
        if (is.numeric(scores)) {
            x <- mean(scores, na.rm=TRUE)
        } else {
            count_list <- as.list(table(scores))
            x <- paste0(
                paste0(names(count_list), ":", unlist(count_list)),
                collapse=","
            )
        }
        return(x)
    }
)

# -- Save files to disk
qsave(ragged_exp, file=output[1], nthreads=params$nthread)
