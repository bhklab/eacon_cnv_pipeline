# 0.1 -- Load dependencies
renv::activate()

library(GenomicRanges)
library(SummarizedExperiment)
library(matrixStats)
library(data.table)
library(qs)


# 0.2 -- Parse snakemake parameters
input <- snakemake@input
params <- snakemake@params
output <- snakemake@output


# 0.3 -- Define helper functions

# Function for Finding Segment Overlaps Between a List of -----------------

#' cnTools: Aggregate Genomic Ranges
#'
#' Function credit goes to Rene Quevedo from the CNTools package maintained by
#'   the Pugh Lab at Princess Margaret Cancer Center
#'   
#' @description Because different samples have different segmentations, this
#'   function attempts to create a unified matrix with the fewest number of
#'   segments needed to represent all copy-number segments from all samples.
#'   This allows for the elementMetadata() to store the associated copy-number
#'   values for those samples all in one matrix.
#' 
#' @param list.gr [GRangesList]: a list of granges objects or GRangesList
#' 
#' @return A single GRanges object with each row of the elementMetadata()
#'  containing all the samples copy-number values for that segment
#'
#' @import GenomicRanges
#' @export
#' 
#' @examples aggregateGr(list(gr1, gr2))
aggregateGr <- function(list.gr){
    list.gr <- as.list(list.gr)
    # Loop to combine the first two GRanges object of the list, save it to the 
    # first element of the list, pop out the second element
    while(length(list.gr) >= 2) { 
        x <- list.gr[[1]]
        y <- list.gr[[2]]
        # Make a composite GRanges object that contains all possible segments
        int.gr <- disjoin(sort(c(granges(x), granges(y))))
        int.gr <- sort(c(int.gr, gaps(int.gr)))

        # Function to create an nrow(int.gr) matrix containing all the mapped 
        # elements from x and y
        xy.list <- lapply(list(x, y), function(z, int.gr){
            mat.blank <- matrix(nrow=length(int.gr),
                                ncol=1)
            olaps <- findOverlaps(int.gr, z, type = 'within')

            mat.fill <- apply(elementMetadata(z), 2, function(meta.z){
                mat.blank[queryHits(olaps), ] <- meta.z[subjectHits(olaps)]
                mat.blank
            })
            as.data.frame(mat.fill, check.names=FALSE)
        }, int.gr=int.gr)

        # Combining the x-y elements back into the intersected GRanges object
        elementMetadata(int.gr) <- cbind(elementMetadata(int.gr), 
            do.call("cbind", xy.list))

        list.gr[[2]] <- NULL
        list.gr[[1]] <- int.gr
    }
    list.gr[[1]]
}

# 1 -- Load the input data
CNV_grList <- qread(input$gr_list)
CNV_bins_sumexp <- qread(input$bins_sumexp)
CNV_genes_sumexp <- qread(input$genes_sumexp)

# 2 -- Save metadata df and drop all mcols from Granges except the the desired output
feature_col <- params$feature_col
sample_annotations <- names(CNV_grList)
CNV_grList_no_meta <- lapply(CNV_grList, 
    function(gr, feature_col) gr[, feature_col], feature_col=feature_col)

# 3 -- Aggregate segmentation results across all samples
CNV_gRanges_aggregated <- aggregateGr(CNV_grList_no_meta)
colnames(mcols(CNV_gRanges_aggregated)) <- sample_annotations

# 4 -- Get segmentation means per sample per region
CNV_segments <- as.data.frame(mcols(CNV_gRanges_aggregated))
CNV_ranges <- as.data.frame(CNV_gRanges_aggregated)
rownames(CNV_segments) <- paste0(CNV_ranges$seqnames, ':', CNV_ranges$start, 
    '-', CNV_ranges$end)
## drop rows with any
CNV_segments <- CNV_segments[!apply(CNV_segments, MARGIN=1, FUN=anyNA), ]

# 5 -- Convert to data.table
CNV_data <- data.table(CNV_segments, keep.rownames = TRUE)

# 6 -- Conditionally drop sex chromosomes
if (params$drop_sex) CNV_data <- CNV_data[grep("chrX|chrY", rn, invert=TRUE), ]

# 7 -- Calculate row variance
CNV_data[, rowMads := rowMads(as.matrix(.SD[, -'rn']))]
setorderv(CNV_data, col='rowMads', order=(-1))

# 8 -- Save data to disk
fwrite(CNV_data[, -'rowMads'], file=output$ranked_feature_file)

## Save subsets of top most variant features for metaclustering
for (i in seq_along(params$feature_numbers)) {
    fwrite(CNV_data[seq_len(params$feature_numbers[i]), -'rowMads'], 
        file=output$feature_number_files[i])
}