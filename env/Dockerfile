FROM bioconductor/bioconductor_docker:devel

ENV R_BIOC_VERSION=3.16

RUN Rscript -e "install.packages('pak')"

# Universal EaCoN dependencies
RUN Rscript -e "pak::pkg_install(c( \
    'sequenza', \
    'affxparser', \
    'Biostrings', \
    'aroma.light', \
    'BSgenome', \
    'copynumber', \
    'GenomicRanges', \
    'limma', \
    'rhdf5'), \
    dependencies=TRUE \
)"
RUN Rscript -e " \
    devtools::install_github('Crick-CancerGenomics/ascat/ASCAT'); \
    devtools::install_github('mskcc/facets') \
"
RUN Rscript -e "devtools::install_github( \
    'gustaveroussy/EaCoN' \
)"


# Affymetrix Array files
RUN Rscript -e "install.packages( \
    'https://zenodo.org/record/5494853/files/affy.CN.norm.data_0.1.2.tar.gz', \
    type='src', \
    repos=NULL \
)"
## Affymetrix OncoScan
RUN Rscript -e "devtools::install_github( \
    'gustaveroussy/apt.oncoscan.2.4.0', \
)"
RUN Rscript -e "install.packages(c( \
    'https://zenodo.org/record/5494853/files/affy.CN.norm.data_0.1.2.tar.gz', \
    'https://zenodo.org/record/5494853/files/OncoScan.na33.r4_0.1.0.tar.gz', \
    'https://zenodo.org/record/5494853/files/OncoScanCNV.na33.r2_0.1.0.tar.gz', \
    'https://zenodo.org/record/5494853/files/OncoScan.na36.r1_0.1.0.tar.gz', \
    'https://zenodo.org/record/5494853/files/OncoScanCNV.na36.r1_0.1.0.tar.gz'), \
    type='src', \
    repos=NULL \
)"
## Affymetrix CytoScan
RUN Rscript -e "devtools::install_github( \
    'gustaveroussy/apt.cytoscan.2.4.0', \
)"
RUN Rscript -e "install.packages(c( \
    'https://zenodo.org/record/5494853/files/CytoScan750K.Array.na33.r4_0.1.0.tar.gz', \
    'https://zenodo.org/record/5494853/files/CytoScanHD.Array.na33.r4_0.1.0.tar.gz', \
    'https://zenodo.org/record/5494853/files/CytoScan750K.Array.na36.r1_0.1.0.tar.gz', \
    'https://zenodo.org/record/5494853/files/CytoScanHD.Array.na36.r1_0.1.0.tar.gz', \
    'https://zenodo.org/record/5494853/files/rcnorm_0.1.5.tar.gz'), \
    type='src', \
    repos=NULL \
)"

# Genomewide SNP6 Array files
RUN Rscript -e "devtools::install_github( \
    'gustaveroussy/apt.snp6.1.20.0' \
)"
RUN Rscript -e "install.packages(c( \
    'https://zenodo.org/record/5494853/files/GenomeWideSNP.6.na35.r1_0.1.0.tar.gz', \
    'https://zenodo.org/record/5494853/files/rcnorm_0.1.5.tar.gz'), \
    type='src', \
    repos=NULL \
)"

# Genomes
RUN Rscript -e "pak::pkg_install(c( \
    'BSgenome.Hsapiens.UCSC.hg19', \
    'BSgenome.Hsapiens.UCSC.hg38', \
    'BSgenome.Hsapiens.1000genomes.hs37d5', \
))"