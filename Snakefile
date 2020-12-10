# -- 0.1 Dependeny imports
import os
import re
import glob
import pandas as pd

## TODO:: Check that workflow dependencies are installed

print(os.getcwd())

# -- 0.2 Load configuration files
configfile: 'config.yaml'

rawdata = config['rawdata_dir']
metadata = config['metadata_dir']
procdata = config['procdata_dir']
pairs_file = config['pairs_file']

pairs_df = pd.read_csv(os.path.join(f'{metadata}', f'{pairs_file}'), sep='\t')


# -- 1. Batch processing of raw CEL or BAM files
rule batch_process_rawdata:
    input:
        pairs_file=os.path.join(metadata, pairs_file)
    params:
        analysis_name=config['analysis_name'],
        array_family=config['array_family'],
        procdata=procdata,
        sample_name=pairs_df.SampleName
    threads: config['nthreads']
    output:
        expand('{procdata}/{sample_name}/1.{analysis_name}_{sample_name}_{array_family}_processed.RDS', 
            procdata=procdata, sample_name=pairs_df.SampleName, array_family=config['array_family'],
            analysis_name=config['analysis_name'])
    script:
        'scripts/1_batchProcessRawdataFiles.R'

# -- 2. L2R and BAF join segmentation
rule segment_processed_data:
    input:
        expand('{procdata}/{sample_name}/1.{analysis_name}_{sample_name}_{array_family}_processed.RDS', 
            procdata=procdata, sample_name=pairs_df.SampleName, array_family=config['array_family'],
            analysis_name=config['analysis_name'])
    params:
        analysis_name=config['analysis_name'],
        array_family=config['array_family'],
        procdata=procdata
    threads:
    output:

    script:
        'scripts/2_segmentProcessedData'