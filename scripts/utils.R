#' Create a GenomicRanges object from
#'
#' @param ascn_data `list`
#' @param l2r_data `list`
#'
#' @export
buildGRangesFromASCNAndL2R <- function(ascn_data, l2r_data) {

}

#'
#'
#'
#' @details
#' Adapted from https://github.com/quevedor2/EaCoN/blob/master/R/output_builder.R.
#' Credit to Rene Quevedo
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