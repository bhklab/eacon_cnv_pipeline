#!/usr/bin/R

#### A SCRIPT TO INSTALL THE DEPENDENCIES FOR ONCOSCAN CNV ANALYSIS PIPELINE ####

library(devtools)
library(BiocManager)

# Install list of CRAN packages
pkgCheck <- function(requiredPackages) {
  for (p in requiredPackages) {
    if (!require(p, character.only = TRUE))
      install.packages(p)
    library(p, character.only = TRUE)
  }
}


CRAN <- c("devtools", "BiocManager", "Rcpp", "parallel")
pkgCheck(CRAN)

# Install list of Bioconductor packages
pkgCheckBioC <- function(requiredPackagesBioc) {
  for (p in requiredPackagesBioc) {
    if (!require(p, character.only = TRUE))
      BiocManager::install(p)
    library(p, character.only = TRUE)
  }
}

BioC <- c(
  "affxparser",
  "Biostrings",
  "aroma.light",
  "BSgenome",
  "copynumber",
  "GenomicRanges",
  "limma",
  "rhdf5",
  "sequenza"
)
pkgCheckBioC(BioC)

# Install list of Github packages
pkgCheckGithub <- function(requiredPackagesGithub) {
  for (p in requiredPackagesGithub) {
    if (!require(p, character.only = TRUE))
      devtools::install_github(p)
    library(p, character.only = TRUE)
  }
}

Git <- c("Crick-CancerGenomics/ascat/ASCAT",
         "mskcc/facets",
         "ChristopherEeles/EaCoN")
pkgCheckGithub(Git)

## These need to be installed in order

# All Affymetrix Microarrays (Now ThermoFisher)
install.packages(
  "https://partage.gustaveroussy.fr/pydio_public/083305?dl=true&file=/affy.CN.norm.data_0.1.2.tar.gz",
  repos = NULL,
  type = "source"
)

# ONCOSCAN APT tool
devtools::install_github("gustaveroussy/apt.oncoscan.2.4.0")

# ONCOSCAN CNV hg19
install.packages(
  "https://partage.gustaveroussy.fr/pydio_public/cd59c8?dl=true&file=/OncoScanCNV.na33.r2_0.1.0.tar.gz",
  repos = NULL,
  type = "source"
)

# Check installed genomes for annotation
if (is.null(BSgenome::installed.genomes())) {
  BiocManager::install("BSgenome.Hsapeiens.UCSC.hg19")
} else if (!grepl(BSgenome::installed.genomes(), pattern = "BSgenome.Hsapiens.UCSC.hg19")) {
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
}
