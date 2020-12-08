  ###
# Illumina Epic LMS CNV Analysis with EaCoN
# by Christopher Eeles
#
# Adapted from: 
# ""
#
###

###########
# Current unknowns/challenges:
#
# - How would we go about recentering the l2r values based on nd-waviness-sd?
#  - Is this already done within this CNV pipeline?
# - Why do some samples fail to have a model fit? Can this be interepreted as 
#   a QC metric? Is there a way to fix/improve this?
#
###########

### NOTE ####
# - EaCoN package has been modified and must be downladed from ChristopherEeles/EaCoN on GitHub 
#  - run devtools::install_github('ChristopherEeles/EaCoN')

####
# ---- Load Libaries ----
####
library(optparse)
library(devtools)

####
# --- Command Line Options ----
####

### COMMAND LINE IMPLEMENTATION IS NOT YET COMPLETE!
option_list <- list(
  make_option(c('-c', '--cores'), type = 'integer', default = 1, 
              help = 'Number of cores to parallelize over', metavar = 'integer'),
  make_option(c('-i', '--install'), type = 'logical', default = FALSE,
              help = 'Should dependencies be installed? Requires and internet 
                connection', metavar = 'character'),
  make_option(c("-d", "--dir"), type = "character", default = '.',
              help = "Parent directory path"))
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

####
# ---- Dependency Installation ----
####

# Install dependencies for Easy Copy Number Oncoscan CNV Analysis
if (opt$install) {
  #  Install list of CRAN packages
  pkgCheck <- function(requiredPackages) {
    for (p in requiredPackages) {
      if (!require(p,character.only = TRUE)) install.packages(p)
      library(p, character.only = TRUE)
    }
  }
  
  CRAN <- c("devtools","BiocManager", "Rcpp")
  pkgCheck(CRAN)
  
  # Install list of Bioconductor packages
  pkgCheckBioC <- function(requiredPackagesBioc) {
    for (p in requiredPackagesBioc) {
      if (!require(p,character.only = TRUE)) BiocManager::install(p)
      library(p,character.only = TRUE)
    }
  }
  
  BioC <- c("affxparser", "Biostrings", "aroma.light", "BSgenome", 
            "copynumber", "GenomicRanges", "limma", "rhdf5", "sequenza", 
            "TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db")
  pkgCheckBioC(BioC)
  
  # Install list of Github packages
  ##FIXME:: This doesn't work sometimes; further testing needed
  pkgCheckGithub <- function(requiredPackagesGithub) {
    for (p in requiredPackagesGithub){
      if (!require(p,character.only = TRUE)) devtools::install_github(p)
      library(gsub('.*/', '', p), character.only = TRUE)
    }
  }
  
  Git <- c("Crick-CancerGenomics/ascat/ASCAT", "mskcc/facets", 
           "quevedor2/EaCoN")
  pkgCheckGithub(Git)
  
  ## These need to be installed in order
  
  # All Affymetrix Microarrays (Now ThermoFisher)
  install.packages("https://partage.gustaveroussy.fr/pydio_public/083305?dl=true&file=/affy.CN.norm.data_0.1.2.tar.gz", repos = NULL, type = "source")
  
  # ONCOSCAN APT tool
  devtools::install_github("gustaveroussy/apt.oncoscan.2.4.0")
  
  # ONCOSCAN CNV hg19
  install.packages("https://partage.gustaveroussy.fr/pydio_public/cd59c8?dl=true&file=/OncoScanCNV.na33.r2_0.1.0.tar.gz", repos = NULL, type = "source")
  
  # Check for installed genomes for annotation
  if (is.nulBSgenome::installed.genomes()) {
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
  } else if (!grepl(BSgenome::installed.genomes(), pattern = "BSgenome.Hsapiens.UCSC.hg19")) {
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
  }
} else {
  library(EaCoN)
  library(ASCAT)
  library(BSgenome)
  library(sequenza)
  library(limma)
  library(affxparser)
  library(Biostrings)
  library(rhdf5)
  library(parallel)
  library(Biobase)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(PharmacoGx)
  library(qs)
}

####
# ---- Data Import and Preprocessing ----
####

## ---- Configure the data directory
baseDir <- opt$dir
metaDir <- paste0(opt$dir, '/data/metadata')
cores <- opt$cores
dataDir <- paste0(opt$dir, '/data/raw_data')
outDir <- paste0(opt$dir, '/analysis/1_pipeline_results')


## ---- Assemble the pairs.file to instruct batch processing

# Read in meta data file
CNV_metadata <- read.csv(list.files(metaDir, pattern='*csv', full.names=TRUE), 
                         stringsAsFactors = FALSE, header = TRUE)

# Create dataframe with columns ATChannelCel, GCChannelCEL and SampleName
oncoscan_pairs <- data.frame(
                        'ATChannelCel' = file.path(dataDir, CNV_metadata$AT.Array),
                        'GCChannelCel' = file.path(dataDir, CNV_metadata$GC.Array), 
                        'SampleName' = CNV_metadata$Sample.Name)

# Add .CEL to the end of dataframe
oncoscan_pairs[, seq_len(2)] <- sapply(oncoscan_pairs[, seq_len(2)], 
                                       function(x) { paste0(x, ".CEL")})

# Save .tsv file
write.table(oncoscan_pairs, file = paste0(metaDir, "/CEL_file_pairs.tsv"), sep = "\t")

## ---- Normalization ##

## Batch process the .CEL files by sample
## Writes results to results directory
OS.Process.Batch(
  pairs.file = file.path(metaDir, "CEL_file_pairs.tsv"),
  nthread = cores,
  out.dir = outDir,
  force = T
)

## ---- L2R & BAF Segmentation ####
Segment.ff.Batch(
  nthread = cores,
  force = T
)

## ---- Copy-number Estimation ####
ASCN.ff.Batch(
  nthread = cores,
  force = T
)

## ---- Select Optimal Gamma Value per Sample ####

# Get per sample gamma reports
gammaFiles <- list.files(path = outDir, pattern = '.*.gammaEval.txt', recursive = T)

# Load PancanPloidy as reference data
data(pancanPloidy.noSegs)
pancan.ploidy <- pancan.obj.segless$breaks$ploidy

# Assemble the meta data
meta_data <- pancanPloidy.noSegs
  CNV_metadata %>% 
  filter(Sample.Name %in% gsub('/.*', '', gammaFiles)) %>% 
  select(Sample.Name, Well)
colnames(meta_data) <- c('sample', 'TCGA')


# Get a list of all gamme files for each sample with a successful fit
all.fits <- lapply(gammaFiles, function(file) {
  list(fit = read.table(file.path(outDir, file), sep = '\t', header = T, stringsAsFactors = F, check.names = F, fill = F),
       sample = strsplit(file, '/')[[1]][1]
  )
})
names(all.fits) <- sapply(all.fits, function(x) x$sample)

# Change into results directory
old_dir <- getwd()
setwd(outDir) ##FIXME:: Do we still need to change directories if I give the absoluate path?
on.exit(setwd(old_dir))

# Select the optimal gamma value and return the results for that model
gr.cnv <- annotateRDS.Batch(
  all.fits,
  'ASCAT',
  nthread = cores,
  gamma.method = 'score',
  gamma.meta = meta_data,
  pancan.ploidy = pancan.ploidy,
  feature.set = c('bins', 'tads')
)

# Save raw results object to disk
qsave(gr.cnv, file = file.path(outDir, 'optimal_gamma_list.qs'), nthread=opt$cores)

## ---- Construct and output eSet objects containing the results and write to disk
buildPSetOut(gr.cnv, 'VDR', file.path(outDir), meta = meta_data)

## ---- Create GRangesList of segmentation results and save them to disk

# Extract segmentation data.frames from the resutls object
segmentation_df_list <- lapply(gr.cnv, function(x) x$seg)

# Convert all data.frames to GRanges and return in a list
list_of_gRanges <- lapply(segmentation_df_list, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = TRUE))

# Convert list of GRanges objects into GRangesList object
cnv_grList <- GRangesList(list_of_gRanges)

# Save GRangesList to disk for downstream analysis
qsave(cnv_grList, file = file.path(outDir, 'LMS80_CNV_segments.qs'), nthread=opt$cores)
