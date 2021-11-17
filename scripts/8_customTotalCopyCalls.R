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


# 1 -- Read in Bioconductor objects
is_na_input <- is.na(unlist(input))
if (any(!is_na_input)) {
    data_list <- setNames(lapply(input[!is_na_input], qread), 
        names(input)[!is_na_input])
}


# 2 -- Assign custom TCN ranges
tcn_cutoffs <- params$tcn_cutoffs
tcn_cutoffs <- lapply(tcn_cutoffs, function(x) x[order(names(x))])
for (i in seq_along(data_list)) {
    object <- data_list[[i]]
    # parse the object log2r to a data.table
    if (is(object, "GRangesList")) {
        flat_grl <- unlist(object)
        mcol_df <- as.data.table(as.data.frame(mcols(flat_grl)))
    } else {
        mcol_df <- as.data.table(assays(object)$exprs, keep.rownames=TRUE)
        mcol_df <- melt.data.table(mcol_df, id.vars="rn", 
            variable.name="sample", value.name="seg.mean")
    }
    # assign tcn calls to specified rangess
    for (i in seq_along(tcn_cutoffs)) {
        new_col <- paste0("tcn_", names(tcn_cutoffs)[i])
        for (j in seq_along(tcn_cutoffs[[i]])) {
            cut_offs <- tcn_cutoffs[[i]][[j]]
            mcol_df[
                seg.mean > cut_offs[[1]] & seg.mean <= cut_offs[[2]],
                (new_col) := as.numeric(names(tcn_cutoffs[[i]])[j])
            ]
        }
        if (!all(unique(mcol_df[[new_col]]) %in% 
                as.numeric(names(tcn_cutoffs[[i]])))) {
            stop("Custom copy number ranges failed to match some values!")
        }
    }
    # update the object with the new tcn calls
    if (is(object, "GRangesList")) {
        mcols(flat_grl) <- mcol_df
        object <- relist(flat_grl, object)
    } else {
        new_cols <- setdiff(colnames(mcol_df), c("rn", "sample", "seg.mean"))
        for (col in new_cols) {
            new_assay <- dcast(mcol_df, rn ~ sample, value.var=col)
            assays(object)[[col]] <- as.matrix(new_assay[, -1], 
                rownames=new_assay[[1]])
        }
    }
    data_list[[i]] <- object
}


# 3 -- Write new objects to disk
is_na_output <- is.na(unlist(output))
if (any(!is_na_output)) {
    for (i in seq_along(data_list)) {
        qsave(data_list[[i]], file=output[!is_na_output][i])
    }
}

