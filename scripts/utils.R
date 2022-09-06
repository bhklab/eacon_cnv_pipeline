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
#' the `ASCAT` R pacakge. See details for method information.
#'
#' @details
#' Adapted from https://github.com/quevedor2/EaCoN/blob/master/R/output_builder.R.
#' Credit to Rene Quevedo.
#'
#' Information about ASCAT method used to call allele specific copy number (ASCN).
#'
#' @references
#' Van Loo et al. Allele-specific copy number analysis of tumors. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2947907/.
#' Ross et al. Allele-specific multi-sample copy number segmentation in ASCAT. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8317109/
#'
#' @return `GenomicRanges` object, containing segment genomic co-ordinates and
#'   associated allele specific copy number calls as well as gain/loss/normal
#'   copy state annotations.
#'
#' @export
buildGRangesFromASCN <- function(ascn_data) {

}


#'
#'
#'
#' @details
#' Adapted from https://github.com/quevedor2/EaCoN/blob/master/R/output_builder.R.
#' Credit to Rene Quevedo
#'
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
    normal <- !gain & !losses

}



#'
#'
#'
#' @export
annotateGRangesWithTxDB <- function(granges, txdb) {

}