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

pairs_df = pd.read_csv(os.path.join(f'{metadata}', f'{pairs_file}'))


# -- 1. Batch processing of raw CEL or BAM files
rule batch_process_rawdata:
    input:
        pairs_file=os.path.join(metadata, pairs_file)
    params:
        array_family=config['array_family'],
        procdata=procdata
    threads: config['nthreads']
    # output:
    #     expand('{procdata}/{sample_name}', 
    #         procdata=procdata, sample_name=pairs_df.SampleName)
    script:
        'scripts/1_batchProcessRawdataFiles.R'

# -- 2. Segment Output from rule 1
