# 0.1 Load dependencies
renv::activate()
library(data.table)


# -- 0.2 Parse snakemake arguments
input <- snakemake@input
params <- snakemake@params
output <- snakemake@output

# -- 0.3 Load local utility functions
source(file.path("scripts", "utils.R"))

# -- 1. Read in gamma files

gammaFiles <- list.files(params$out_dir, '.*gammaEval.*txt', recursive=TRUE,
    full.names=TRUE)

gamma_df <- rbindlist(
    setNames(
        lapply(gammaFiles, FUN=fread),
        gsub(".gammaEval.txt", "", basename(gammaFiles))
    ),
    idcol="sample_name"
)

# -- 2. Select the best fits
best_fits <- gamma_df[, .SD[which.max(GoF)], by="sample_name"]

# -- 3. Save the best fit data.frame to csv
fwrite(best_fits, file=output[[1]])