# Easy Copy Number Analysis (EaCoN) Pipeline

This pipeline has been adapted from https://github.com/gustaveroussy/EaCoN.
It leverages the `EaCoN` R package to conduct preprocessing and normalization,
segmentation and copy number estimation from specific microarray .CEL files.
The `EaCoN` package supports copy number estimation from Cytoscan, Oncoscan,
SNP6 arrays as well as WES data, however, the current pipeline has only
implemented support for microarray data.

## Snakemake

This pipeline leverages the `snakemake` Python package for workflow management.
As a result the pipeline and its dependencies are easily
installed from this repository, allowing quick setup, configuration and
deployment.

For more information on Snakemake, please see:
https://snakemake.readthedocs.io/en/stable/.

## Software Environment

Dependency management for this pipeline is handled via `conda` for Python
and `renv` for R. To get started with setup you can install
miniconda3 using the instructions available here: https://docs.conda.io/en/latest/miniconda.html.

## Setting Up Your Software Environment

The first step to deploying an analysis pipeline is to install the various
software packages it depends on. We have included the `env/eacon.yml` and `renv.lock` files here to easily accomplish this.

All commands should be executed from the top level directory of this
repository unless otherwise indicated.

### Python and Snakemake

The first step to deploying an analysis pipeline is to install Python,
Snakemake and R via `conda`. We have included the
`env/eacon.yml` which specifies all the requisite dependencies to use
`snakemake` for this pipeline.

You can use `conda` to install all Python and OS system dependencies
using:

`conda env create --file env/eacon.yml`

This will take some time to run as it gathers and installs the correct
package versions. The environment it creates should be called `eacon`.

If it is not automatically activated after installation please run
`conda activate eacon` before proceeding to the next step.

### R Dependencies

R dependencies are handled via `renv`  and defined in `renv.lock`.

To initialize this project with `renv` run:

`Rscript -e 'library(renv); renv::init()'`

Alternatively, if you wish to isolate the R dependencies from your Conda environment R libraries, you can use this command instead:

`Rscript -e 'library(renv); renv::isolate(); renv::init(bare=TRUE)'`

If intialization doesn't trigger dependency installation, you can do so manually using:

`Rscript -e 'renv::restore()'`

For more information on `renv` and how it can be used to manage dependencies in
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

Please at minimum create the `rawdata` and `metadata` directories, as they are assumed to hold the raw microarray plate data (.CEL) and the pairs file, respectively. For more information on the correct formatting for your pairs file, please see https://github.com/gustaveroussy/EaCoN.
The remaining missing directories will be created automatically as the pipeline runs.

### config.yaml

This file hold the relevant pipeline documentation. Here you can specify the paths
to all the parameters for your current pipeline use case. Documentation is provided
in the `config.yaml` file on what each field should contain.

## Using the Pipeline

To run the pipeline end-to-end:
```
snakemake --cores <n_cores> --use-conda
```
Where <n_cores> is the number of cores to parallelize over.

Instructions for running each step individually are presented below.

### Batch Processing and Normalization

`snakemake --cores 2 --use-conda batch_process_rawdata`

### Segmentation

`snakemake --cores 2 --use-conda segment_processed_data`

### Copy Number Calling

`snakemake --cores 2 --use-conda estimate_copy_number`

### Determine Optimal Value for Gamma

`snakemake --cores 2 --use-conda select_optimal_gamma`

### Build Bioconductor SummarizedExperiment Objects

`snakemake --cores 2 --use-conda build_summarized_experiments`

### Filter Samples Based on QC Criteria

`snakemake --cores 2 --use-conda sample_quality_control`

### Use Custom Total Copy Number Cutoffs

`snakemake --cores 2 --use-conda custom_total_copy_calls`

### Select Top Features

`snakemake --cores 2--use-conda select_top_features`
