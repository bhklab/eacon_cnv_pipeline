# Easy Copy Number Analysis (EaCoN) Pipeline

This pipeline has been adapted from https://github.com/gustaveroussy/EaCoN. It leverages the `EaCoN` R package to conduct preprocessing and normalization, segmentation and copy number estimation from specific microarray .CEL files. The `EaCoN` package supports copy number estimation from Cytoscan, Oncoscan, SNP6 arrays as well as WES data, however, the current pipeline has only implemented support for microarray data.

## Snakemake

This pipeline leverages the `snakemake` Python package for workflow management.
As a result the pipeline and its dependencies are easily
installed from this repository, allowing quick setup, configuration and
deployment.

For more information on Snakemake, please see:
https://snakemake.readthedocs.io/en/stable/.

## Software Environment

Dependency management for this pipeline is handled via `conda` for Python
and `singularity` for R. To get started with setup you can install
miniconda3 using the instructions available here:
https://docs.conda.io/en/latest/miniconda.html.

All commands should be executed from the top level directory of this
repository unless otherwise indicated.

### Python and Snakemake

The first step to deploying an analysis pipeline is to install Python,
Snakemake and Singularity via `conda`. We have included the
`envs/eacon.yml` which specifies all the requisite dependencies to use
`snakemake` for this pipeline.

You can use `conda` to install all Python and OS system dependencies
using:

`conda env create --file env/eacon.yml`

This will take some time to run as it gathers and installs the correct
package versions. The environent it creates should be called `eacon`.

If it is not automatically activated after installation please run
`conda activate eacon` before proceeding to the next step.

### R Dependencies

R dependencies are handled via `singularity` and all rules in this
pipeline will run within the container defined in `env/Dockerfile`.
`Snakemake` uses `singularity` to run containerized pipelines, and this
dependency was installed automatically via `conda` in the previous step.

We have already built the container image and pushed it to Dockerhub, so you
shouldn't need to build it yourself locally. To run the pipeline, you do not
need to install `docker`, however you will need it if you want to rebuild the
image.

The images we have versioned for this pipeline are assay specific to minimize
the size our our DockerHub images. The options

If for some reason you cannot get the prebuilt container image to work you can
rebuild the image using the command:

`docker build -f env/Dockerfile -t nci_alamanac`

This will build the image required for this pipeline on your local machine.
You must then update the `container` field in `config.yaml` to point to
the image `nci_almanac`.

## Deployment

To run the pipeline end-to-end:
```
snakemake --cores <n_cores> --use-singularity
```
Where <n_cores> is the number of cores to parallelize over.