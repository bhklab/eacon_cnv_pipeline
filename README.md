# Easy Copy Number Analysis (EaCoN) Pipeline

This pipeline has been adapted from https://github.com/gustaveroussy/EaCoN. It leverages the `EaCoN` R package to conduct preprocessing and normalization, segmentation and copy number estimation from raw Cytoscan HD .CEL files. The `EaCoN` package also supports copy number estimation from Oncoscan, SNP6 arrays as well as WES data, however, the current pipeline has only implemented support for Cytoscan and Oncoscan based arrays.

## Snakemake

This pipeline leverages the `snakemake` Python package for workflow management. As a result the pipeline and its dependencies are easily
installed from this repository, allowing quick setup, configuration and deployment.

For more information on Snakemake, please see: https://snakemake.readthedocs.io/en/stable/.

## Requirements

Dependency management for this pipeline is handled via `conda` for Python 
and `renv` for R. To get started with setup you can install
miniconda3 using the instructions available here: https://docs.conda.io/en/latest/miniconda.html. If you do not currently have R installed, you can install it via conda using the command: `conda -c conda-forge r-base==3.6`. Please note that this pipleine has not been updated to work with R >= 4.0.

Alternatively you can install it directly from CRAN
as described here: https://cran.r-project.org/.

## Setting Up Your Software Environment

The first step to deploying an analysis pipeline is to install the various
software packages it depends on. We have included the `envs/affymetrix.yml` and `renv.lock` files here to easily accomplish this.

All commands should be executed from the top level directory of this
repository.

### Python and System Dependencies

Conda can be used to install all Python and most OS system dependencies
using:

`conda env create --file envs/environment.yml`

This will take some time to run as it gathers and installs the correct
package versions. The environent it creates should be called `affymetrix`.

If it is not automatically activated after installation please run 
`conda activate affymetrix` before proceeding to the next step.

### R Dependencies

The `renv` package can be used to install all R dependencies (both CRAN and
Bioconductor). R version 3.6.3 and `renv` are included as dependencies in the `environment.yml` file and should be installed automatically when setting up your conda environment.

To initialize this project with renv run:

`Rscript -e 'library(renv); renv::init()'`

If you wish to isolate the R dependencies from your Conda environment R libraries, you can use this command instead:

`Rscript -e 'library(renv); renv::isolate(); renv::init(bare=TRUE)'`

If intialization doesn't trigger dependency installation, you can do so manually using:

`Rscript -e 'renv::restore()'`

For more information on renv and how it can be used to manage dependencies in
your project, please see: https://rstudio.github.io/renv/articles/renv.html.

## Configuring the Pipeline

This pipeline assumes the following directory structure:

```
.
├── envs
├── metadata
├── procdata
├── rawdata
├── renv
├── results
└── scripts
```

Please at minimum create the `rawdata` and `metadata` directories, as they are assumed to hold the raw Cytoscan HD plate data (.CEL) and the pairs file, respectively. For more information on the correct formatting for your pairs file, please see https://github.com/gustaveroussy/EaCoN.
The remaining missing directories will be created automatically as the pipeline runs.

### config.yaml

This file hold the relevant pipeline documentation. Here you can specify the paths
to all the parameters for your current pipeline use case. Documentation is provided
in the `config.yaml` file on what each field should contain.

## Using the Pipeline

### Batch Processing and Normalization

`snakemake --cores 2 batch_process_rawdata`

### Segmentation 

`snakemake --cores 2 segment_processed_data`

### Copy Number Calling

`snakemake --cores 2 estimate_copy_number`

### Determine Optimal Value for Gamma

`snakemake --cores 2 select_optimal_gamma`

### Build Bioconductor SummarizedExperiment Objects

`snakemake --cores 2 build_summarized_experiments`
