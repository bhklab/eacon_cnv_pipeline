# -- Build arguments, with default values
ARG R_BIOC_VERSION=3.16


## -- Mandatory EaCoN dependencies
FROM bioconductor/bioconductor_docker:devel AS eacon-base
ENV R_BIOC_VERSION=${R_BIOC_VERSION}
RUN Rscript -e "install.packages('pak')"
# Universal EaCoN dependencies
RUN Rscript -e "pak::pkg_install(c( \
    'XML', \
    'sequenza', \
    'affxparser', \
    'Biostrings', \
    'aroma.light', \
    'BSgenome', \
    'copynumber', \
    'GenomicRanges', \
    'limma', \
    'rhdf5' \
))"
RUN Rscript -e " \
    devtools::install_github('Crick-CancerGenomics/ascat/ASCAT'); \
    devtools::install_github('mskcc/facets') \
"
RUN Rscript -e "devtools::install_github( \
    'gustaveroussy/EaCoN' \
)"
# Pipeline infrastructure dependencies
RUN Rscript -e " \
pak::pkg_install(c( \
    'qs', \
    'TxDb.Hsapiens.UCSC.hg19.knownGene' \
))"
RUN Rscript -e "install.packages('data.table')"


## -- Installed reference genomes
FROM eacon-base AS eacon-genomes
RUN Rscript -e " \
pak::pkg_install(c( \
    'BSgenome.Hsapiens.UCSC.hg19', \
    'BSgenome.Hsapiens.UCSC.hg38', \
    'BSgenome.Hsapiens.1000genomes.hs37d5' \
))"


## -- Shared Affymetrix array dependencies
FROM eacon-genomes AS eacon-affymetrix
RUN Rscript -e " \
options(timeout=Inf); \
install.packages(c('DBI', 'RSQLite')); \
install.packages( \
    'https://zenodo.org/record/5494853/files/affy.CN.norm.data_0.1.2.tar.gz', \
    type='src', \
    repos=NULL \
)"


## -- Dependencies for Affymetrix Oncoscan arrays
FROM eacon-affymetrix AS eacon-oncoscan
RUN Rscript -e " \
options(timeout=Inf); \
devtools::install_github( \
    'gustaveroussy/apt.oncoscan.2.4.0', \
)"
RUN Rscript -e " \
options(timeout=Inf); \
install.packages(c( \
    'https://zenodo.org/record/5494853/files/OncoScan.na33.r4_0.1.0.tar.gz', \
    'https://zenodo.org/record/5494853/files/OncoScanCNV.na33.r2_0.1.0.tar.gz', \
    'ht tps://zenodo.org/record/5494853/files/OncoScan.na36.r1_0.1.0.tar.gz', \
    'https://zenodo.org/record/5494853/files/OncoScanCNV.na36.r1_0.1.0.tar.gz'), \
    type='src', \
    repos=NULL \
)"

FROM eacon-oncoscan AS eacon-retry
RUN Rscript -e "\
    options(timeout=Inf); \
    install.packages( \
        'https://zenodo.org/record/5494853/files/affy.CN.norm.data_0.1.2.tar.gz', \
        type='src', \
        repos=NULL \
    ) \
"

# -- Dependencies for Affymetrix Cytoscan 750k and HD arrays
FROM eacon-affymetrix AS eacon-cytoscan
RUN Rscript -e " \
options(timeout=Inf); \
devtools::install_github( \
    'gustaveroussy/apt.cytoscan.2.4.0', \
)"
RUN Rscript -e " \
options(timeout=Inf); \
install.packages(c( \
    'https://zenodo.org/record/5494853/files/CytoScan750K.Array.na33.r4_0.1.0.tar.gz', \
    'https://zenodo.org/record/5494853/files/CytoScanHD.Array.na33.r4_0.1.0.tar.gz', \
    'https://zenodo.org/record/5494853/files/CytoScan750K.Array.na36.r1_0.1.0.tar.gz', \
    'https://zenodo.org/record/5494853/files/CytoScanHD.Array.na36.r1_0.1.0.tar.gz', \
    'https://zenodo.org/record/5494853/files/rcnorm_0.1.5.tar.gz'), \
    type='src', \
    repos=NULL \
)"


# -- Dependencies for Affymetrix SNP6 arrays
FROM eacon-genomes AS eacon-snp6
RUN Rscript -e "devtools::install_github( \
    'gustaveroussy/apt.snp6.1.20.0' \
)"
RUN Rscript -e " \
options(timeout=Inf); \
install.packages(c( \
    'https://zenodo.org/record/5494853/files/GenomeWideSNP.6.na35.r1_0.1.0.tar.gz', \
    'https://zenodo.org/record/5494853/files/rcnorm_0.1.5.tar.gz'), \
    type='src', \
    repos=NULL \
)"


# -- Dummy rule for Whole Exome Sequencing data, which has no external dependencies
FROM eacon-genomes AS eacon-wes
RUN echo "All dependencies already installed!"