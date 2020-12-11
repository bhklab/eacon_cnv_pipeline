# -- 0.1 Dependeny imports
import os
import re
import glob
import pandas as pd

## TODO:: Check that workflow dependencies are installed

# -- 0.2 Load configuration files
configfile: 'config.yaml'

rawdata = config['rawdata_dir']
metadata = config['metadata_dir']
procdata = config['procdata_dir']
pairs_file = config['pairs_file']

pairs_df = pd.read_csv(os.path.join(f'{metadata}', f'{pairs_file}'), sep='\t')
reference = config['reference']
ref_symbol = reference.split('.')[-1]
nthreads=config['nthreads']


# -- 1. Batch processing of raw CEL or BAM files
rule batch_process_rawdata:
    input:
        pairs_file=os.path.join(metadata, pairs_file)
    params:
        analysis_name=config['analysis_name'],
        array_type=config['array_type'],
        procdata=procdata,
        reference=reference
    threads: nthreads
    output:
        expand('{procdata}/{sample_name}/{sample_name}_{array_type}_Array_{ref_symbol}_processed.RDS', 
            procdata=procdata, sample_name=pairs_df.SampleName, array_type=config['array_type'],
            analysis_name=config['analysis_name'], ref_symbol=ref_symbol)
    script:
        'scripts/1_batchProcessRawdataFiles.R'

# -- 2. L2R and BAF join segmentation
## TODO:: parameterize the segmentation algorithm

rule segment_processed_data:
    input:
        expand('{procdata}/{sample_name}/{sample_name}_{array_type}_Array_{ref_symbol}_processed.RDS', 
            procdata=procdata, sample_name=pairs_df.SampleName, array_type=config['array_type'],
            ref_symbol=ref_symbol)
    params:
        procdata=procdata,
        segmenter=config['segmenter'],
        smoothk=config['smoothk']
    threads: nthreads
    output:
        expand('{procdata}/{sample_name}/{segmenter}/L2R/{sample_name}.SEG.{segmenter}.RDS',
            procdata=procdata, sample_name=pairs_df.SampleName, segmenter=config['segmenter'])
    script:
        'scripts/2_segmentProcessedData.R'

# -- 3. Estimate copy number using the method appropriate to the selected segmenter from the previous step
## TODO:: Try parallelizing in Snakemake instead of in R? Will make deploying on cluster easier
rule estimate_copy_number:
    input:
        expand('{procdata}/{sample_name}/{segmenter}/L2R/{sample_name}.SEG.{segmenter}.RDS',
            procdata=procdata, sample_name=pairs_df.SampleName, segmenter=config['segmenter'])
    params:
        procdata=procdata,
        gamma_range=config['gamma_range']
    threads: nthreads
    output: 
        expand('{procdata}/{sample_name}/{segmenter}/ASCN/gamma_{gamma_value}/{sample_name}.SEG.{segmenter}.RDS',
            procdata=procdata, sample_name=pairs_df.SampleName, segmenter=config['segmenter'],
            gamma_value=range(config['gamma_range'][0], config['gamma_range'][1] + 0.5, 0.5))
    script:
        'scripts/3_estimateCopyNumber.R'