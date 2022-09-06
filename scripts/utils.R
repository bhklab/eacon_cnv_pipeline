#' Create a GenomicRanges object from
#'
#' @param ascn_data `list`
#' @param l2r_data `list`
#'
#' @export
buildGRangesFromASCNAndL2R <- function(ascn_data, l2r_data) {

}



#' Make a GenomicRanges object from the output of the ASCAT R package
#'
#' @param ascn_data `list` Output from `EaCoN::ASCN.ff` function, computed via
#' the `ASCAT` R package. See details for method information.
#'
#' @details
#' Adapted from https://github.com/quevedor2/EaCoN/blob/master/R/output_builder.R.
#' Credit to Rene Quevedo.
#'
#' Information about ASCAT method used to call allele specific copy number (ASCN)
#'   is available in provided references.
#'
#' @references
#' Van Loo et al., 2010. Allele-specific copy number analysis of tumors. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2947907/.
#' Ross et al., 2021. Allele-specific multi-sample copy number segmentation in ASCAT. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8317109/
#'
#' @return `GenomicRanges` object, containing segment genomic co-ordinates and
#'   associated allele specific copy number calls.
#'
#' @authors
#' Rene Quevedo <rene.quevedo at uhnresearch.ca>
#' Christopher Eeles <christopher.eeles at uhnresearch.ca>
#'
#' @md
#' @export
buildGRangesFromASCN <- function(ascn_data) {

}



#' Make a GenomicRanges object from the segment Log2-ratios
#'
#' @param l2r_data `list` Output from `EaCoN::Segment.ff` function, computed via
#'   circular binary segmentation.
#'
#' @details
#' Adapted from https://github.com/quevedor2/EaCoN/blob/master/R/output_builder.R.
#' Credit to Rene Quevedo.
#'
#' @return `GenomicRanges` object, containing segment genomic co-ordinates and
#'   associated log2-ratios as well as gain/loss/normal copy state annotations.
#'
#' @authors
#' Rene Quevedo <rene.quevedo at uhnresearch.ca>
#' Christopher Eeles <christopher.eeles at uhnresearch.ca>
#'
#' @md
#' @export
buildGRangesFromL2R <- function(l2r_data) {
    # extract per segment Log2Ratio
    l2r_segments <- l2r_data$cbs$nocut
    # get gain/loss cutoffs
    gain_cutoff <- l2r_data$meta$eacon[["L2R-segments-gain-cutoff"]]
    loss_cutoff <- l2r_data$meta$eacon[["L2R-segments-loss-cutoff"]]
    # identify segment copy state
    gains <- l2r_segments$Log2Ratio > gain_cutoff
    losses <- l2r_segments$Log2Ratio < loss_cutoff
    l2r_segments$copy_state <- "normal"
    l2r_segments$copy_state[gains] <- "gain"
    l2r_segments$copy_state[losses] <- "loss"
    stopifnot(
        all(unique(l2r_segments$copy_state) %in% c("normal", "gain", "loss"))
    )
    colnames(l2r_segments) <- c("sample", "chr", "start", "end", "probes",
        "log2r", "copy_state")
    rename_chromomsomes <- function(n) switch(n, "23"="X", "24"="Y", n)
    l2r_segments$chr <- vapply(as.character(l2r_segments$chr),
        FUN=rename_chromomsomes, FUN.VALUE=character(1))
    granges <- GenomicRanges::makeGRangesFromDataFrame(l2r_segments,
        keep.extra.columns=TRUE)
    metadata(granges) <- l2r_data$meta[c("basic", "eacon")]
    return(granges)
}



#'
#'
#'
#' @export
annotateGRangesWithTxDB <- function(granges, txdb) {

}