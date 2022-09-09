#' Create a GenomicRanges object from the output of segmentation and
#' allele-specific copy number estimation
#'
#' @param ascn_data `list` Output from `EaCoN::ASCN.ff` function, computed via
#'   the `ASCAT` R package.
#' @param l2r_data `list` Output from `EaCoN::Segment.ff` function, computed via
#'   circular binary segmentation.
#'
#' @details
#' See `?buildGRangesFromL2R` and `buildGRangesFromASCN` for details of
#' GenomicRange creation. This helper method simply calls `merge` on the
#' `GRange` objecst returned each function call.
#'
#' @export
buildGRangesFromASCNAndL2R <- function(ascn_data, l2r_data) {
    l2r_granges <- buildGRangesFromL2R(l2r_data)
    ascn_granges <- buildGRangesFromASCN(ascn_data)
    cnv_granges <- merge(l2r_granges, ascn_granges)
    return(cnv_granges)
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
    # build the GenomicRanges object
    ascn_segments <- ascn_data$segments_raw
    colnames(ascn_segments) <- gsub("pos", "", colnames(ascn_segments))
    ascn_segments$TCN <- ascn_segments$nMinor + ascn_segments$nMajor
    granges <- GenomicRanges::makeGRangesFromDataFrame(ascn_segments,
        keep.extra.columns=TRUE)
    # unify chromosome nomenclature
    seqlevels(granges) <- c(1:22, "X", "Y")
    seqlevelsStyle(granges) <- "UCSC"
    return(granges)
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
    # make the GenomicRanges object
    colnames(l2r_segments) <- c("sample", "chr", "start", "end", "probes",
        "log2r", "copy_state")
    rename_chromomsomes <- function(n) switch(n, "23"="X", "24"="Y", n)
    l2r_segments$chr <- vapply(as.character(l2r_segments$chr),
        FUN=rename_chromomsomes, FUN.VALUE=character(1))
    granges <- GenomicRanges::makeGRangesFromDataFrame(l2r_segments,
        keep.extra.columns=TRUE)
    metadata(granges) <- l2r_data$meta[c("basic", "eacon")]
    # unify chromosome nomenclature
    seqlevels(granges) <- c(1:22, "X", "Y")
    seqlevelsStyle(granges) <- "UCSC"
    return(granges)
}



#' Add gene annotations to a GRanges object using a TxDB object
#'
#' @description
#' Annotates the ranges with Ensembl, Entrez and Hugo Symbol based on the
#' genomic co-ordinates in the GRanges object. Ranges with more the one feature
#' are collapsed into a single string delimited with "|".
#'
#' @param granges `GRanges`
#' @param txdb `TxDB`
#' @param keytype `character(1)` The identifier column name from `Org.hs.eg.db`
#' which corresponds to the `gene_id` column of `genes(txdb)`. Default is
#' "ENTREZID", which is the identifier
#'
#' @return `GRanges` Original object with mapped gene annotations added as
#' metadata columns.
#'
#' @import org.Hs.eg.db
#' @export
annotateGRangesWithTxDB <- function(granges, txdb, keytype="ENTREZID", ...) {
    # validate input
    if (!require(org.Hs.eg.db)) stop("Please install org.Hs.eg.db!")
    stopifnot(is(granges, "GRanges") && is(txdb, "TxDb"))
    # extract the gene annotations from the TxDB object
    gene_annot <- genes(txdb)
    # intersect with the genomic ranges object
    olaps <- as.data.frame(findOverlaps(granges, gene_annot))
    # retrieve additional gene annotations from org.Hs.eg.db
    gene_ids <- gene_annot$gene_id[olaps$subjectHits]
    cols <- c("ENSEMBL", "ENTREZID", "SYMBOL")
    gene_labels <- select(org.Hs.eg.db, keys=gene_ids, keytype=keytype,
        columns=cols, multi="first")
    # merge additional annotations with TxDB annotations
    .paste_or <- function(x) if (!all(is.na(x))) paste0(unique(na.omit(x)), collapse="|") else unique(x)
    segment_genes <- data.frame(gene_id=gene_ids, segment=olaps$queryHits)
    segment_genes <- merge(segment_genes, gene_labels, by.x="gene_id",
        by.y=keytype, all.x=TRUE)
    colnames(segment_genes)[colnames(segment_genes)  == "gene_id"] <- keytype
    # collapse multiple features per GRanges segment
    segment_genes <- aggregate(segment_genes[, -2], FUN=.paste_or,
        by=list(segment=segment_genes$segment))
    # add metadata for segments with gene annotations
    mcols(granges)[segment_genes$segment, colnames(segment_genes[, -1])] <- segment_genes[, -1]
    return(granges)
}


#' Uses the chromosome lengths from a BSGenome to divide al chromosomes in the
#' genome into equal sized bins.
#'
#' @param bin_size `integer(1)` Size of the bins to split chromosomes into.
#' Default is `50000L`.
#' @param seq_style `character(1)` Chromosome naming style. Default is "UCSC".
#'
#' @return `GRanges` object with each chromosome divided into bins of size
#' `bin_size`.
#'
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @export
binReferenceGenome <- function(bin_size=50000L, seq_style="UCSC", reference="BSgenome.Hsapiens.UCSC.hg19"){
    if (!require(reference, character.only=TRUE))
        stop("Please install ", reference)
    chrs <- seqlengths(get(reference))[paste0("chr", c(1:22,"X", "Y"))]

    ## Construct intervals across the genome of a certain bin size
    start_points <- seq(1, max(chrs), by=bin_size)
    grList <- lapply(names(chrs), function(chr_id){
        chr <- chrs[chr_id]
        iranges <- IRanges(start=start_points[start_points < chr], width=bin_size)
        end(iranges[length(iranges),]) <- chr
        granges <- GRanges(seqnames = chr_id, iranges)
        granges
    })

    ## Assemble all GRanges and set seq level style
    grList <- as(grList, "GRangesList")
    suppressWarnings(seqlevelsStyle(grList) <- seq_style)
    granges <- unlist(grList)
    return(granges)
}