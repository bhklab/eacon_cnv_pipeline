# 0.1 Load dependencies
library(EaCoN)
library(data.table)

# -- 0.2 Parse snakemake arguments
# input <- snakemake@input
# params <- snakemake@params
# output <- snakemake@output


input <- list(out_dir='procdata')
params <- list(threads=16)
output <- list(results='results')

# -- 1. Read in gamma files

gammaFiles <- list.files(input$out_dir, '.*gammaEval.*txt', recursive=TRUE, 
    full.names=TRUE)

# -- 2. Load pancanPloidy data as reference
data(pancanPloidy.noSegs)
pancan.ploidy <- pancan.obj.segless$breaks$ploidy

all.fits <- lapply(gammaFiles, function(file) {
  list(fit = fread(file, sep = '\t'),
       sample = strsplit(file, '/')[[1]][2]
  )
})
names(all.fits) <- sapply(all.fits, function(x) x$sample)

.fitsToVector <- function(fit.list) {
    newList <- list()
    nameList <- list()
    for (i in seq_along(fit.list)) {
        fit <- as.list(fit.list[[i]]$fit[[2]])
        names(fit) <- tolower(unlist(fit.list[[i]]$fit[[1]]))
        print(fit)
        newList[[i]] <- fit
        nameList <- c(nameList, fit.list[[i]][[2]])
    }
    names(newList) <- unlist(nameList)
    return(newList)
}

vec.fits <- .fitsToVector(all.fits)

# -- 3. Annotated the RDS data associated with the gamma files

# Change to fix error in annotaRDS.Batch
setwd('procdata')

gr.cnv <- EaCoN:::annotateRDS.Batch(
  all.fits,
  'ASCAT',
  nthread = params$threads,
  gamma.method = 'score',
  pancan.ploidy = pancan.ploidy
)

setwd('..')

# Save raw results object to disk
qsave(gr.cnv, file = file.path(input$out_dir, 'optimal_gamma_list.qs'), nthread=params$threads)

## ---- Construct and output eSet objects containing the results and write to disk
buildPSetOut(gr.cnv, 'VDR', file.path(input$out_dir), meta = meta_data)

## ---- Create GRangesList of segmentation results and save them to disk

# Extract segmentation data.frames from the resutls object
segmentation_df_list <- lapply(gr.cnv, function(x) x$seg)

# Convert all data.frames to GRanges and return in a list
list_of_gRanges <- lapply(segmentation_df_list, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = TRUE))

# Convert list of GRanges objects into GRangesList object
cnv_grList <- GRangesList(list_of_gRanges)

# Save GRangesList to disk for downstream analysis
qsave(cnv_grList, file = file.path(output$results, 'fletcher_LMS_grList.qs'), nthread=params$threads)