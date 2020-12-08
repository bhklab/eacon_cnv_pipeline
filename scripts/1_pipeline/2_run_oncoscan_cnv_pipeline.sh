#!/bin/bash
#SBATCH --job-name=oncoscan_cnv_pipeline
#SBATCH --mem=30G
#SBATCH -c 32
#SBATCH --output=oncoscan_cnv_pipeline_v1.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ceeles

echo 'Starting clustering'

source /cluster/home/ceeles/.bashrc

baseDir="$bhk/data/vanDeRijn/van_de_Rijn_samples/1c_oncoscan_cnv_analysis"

conda activate rEnv
module load R

Rscript $baseDir/scripts/1_pipeline/1_oncoscan_CNV_pipeline.R -d $baseDir -c 31

echo 'JOB COMPLETE!'