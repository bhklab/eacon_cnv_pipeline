# 0.1 Load dependencies
library(EaCoN)
library(data.table)

# -- 0.2 Parse snakemake arguments
input <- snakemake@input
params <- snakemake@params
output <- snakemake@output

if (missing(input) < 1) input <- list(out_dir='procdata')
if (missing(params) < 1) params <- list(threads=16)

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
## FIXME:: Why is the function not exporting despite correct roxygen2 docs?

# Change to fix error in annotaRDS.Batch
setwd('procdata')

gr.cnv <- EaCoN:::annotateRDS.Batch(
  all.fits,
  'ASCAT',
  nthread = params$threads,
  gamma.method = 'score',
  pancan.ploidy = pancan.ploidy,
  feature.set = c('bins', 'tads')
)


# gr.cnv <- EaCoN:::annotateRDS(
#   all.fits[[1]]$fit,
#   all.fits[[1]]$sample,
#   'ASCAT',
#   gamma.method = 'score',
#   pancan.ploidy = pancan.ploidy,
#   feature.set = c('bins', 'tads')
# )