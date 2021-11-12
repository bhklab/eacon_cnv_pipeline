# -- 0.1 Dependeny imports
import os
import re
import glob
import pandas as pd
import numpy as np

## TODO:: Check that workflow dependencies are installed

# -- 0.2 Load configuration files
configfile: "config.yaml"

rawdata = config["rawdata_dir"]
metadata = config["metadata_dir"]
procdata = config["procdata_dir"]
pairs_file = config["pairs_file"]

pairs_df = pd.read_csv(os.path.join(f"{metadata}", f"{pairs_file}"), sep="\t")
reference = config["reference"]
ref_symbol = reference.split(".")[-1]
nthreads=config["nthreads"]


# -- 1. Batch processing of raw CEL or BAM files
## FIXME:: Whether CNV or Array go in the file name is assay dependent!
rule batch_process_rawdata:
    input:
        pairs_file=os.path.join(metadata, pairs_file)
    params:
        analysis_name=config["analysis_name"],
        array_type=config["array_type"],
        procdata=procdata,
        reference=reference,
        rawdata=rawdata
    threads: nthreads
    output:
        expand("{procdata}/{sample_name}/{sample_name}_{array_type}_{ref_symbol}_processed.RDS", 
            procdata=procdata, sample_name=pairs_df.SampleName, array_type=config["array_type"],
            analysis_name=config["analysis_name"], ref_symbol=ref_symbol)
    script:
        "scripts/1_batchProcessRawdataFiles.R"


# -- 2. L2R and BAF join segmentation

rule segment_processed_data:
    input:
        expand("{procdata}/{sample_name}/{sample_name}_{array_type}_{ref_symbol}_processed.RDS", 
            procdata=procdata, sample_name=pairs_df.SampleName, array_type=config["array_type"],
            ref_symbol=ref_symbol)
    params:
        procdata=procdata,
        segmenter=config["segmenter"],
        smoothk=config["smoothk"],
        nrf=config["nrf"],
        SER_pen=config["SER_pen"],
        BAF_filter=config["BAF_filter"]
    threads: nthreads
    output:
        expand("{procdata}/{sample_name}/{segmenter}/L2R/{sample_name}.SEG.{segmenter}.RDS",
            procdata=procdata, sample_name=pairs_df.SampleName, segmenter=config["segmenter"])
    script:
        "scripts/2_segmentProcessedData.R"


# -- 3. Estimate copy number using the method appropriate to the selected segmenter from the previous step
## TODO:: Try parallelizing in Snakemake instead of in R? Will make deploying on cluster easier
rule estimate_copy_number:
    input:
        expand("{procdata}/{sample_name}/{segmenter}/L2R/{sample_name}.SEG.{segmenter}.RDS",
            procdata=procdata, sample_name=pairs_df.SampleName, segmenter=config["segmenter"])
    params:
        procdata=procdata,
        gamma_range=config["gamma_range"]
    threads: nthreads
    output: 
        # expand("{procdata}/{sample_name}/{segmenter}/ASCN/gamma{gamma_value}/{sample_name}.ASCN.{segmenter}.RDS",
        #     procdata=procdata, sample_name=pairs_df.SampleName, segmenter=config["segmenter"],
        #     gamma_value=np.round(np.arange(config["gamma_range"][0], config["gamma_range"][1] + 0.05, 0.05), 2))
    script:
        "scripts/3_estimateCopyNumber.R"


# -- 4. Select the optimal gamma value for each sample
analysis_name = config["analysis_name"]
results_dir = config["results_dir"]

rule select_optimal_gamma:
    input:
        out_dir=procdata
    params:
        nthreads=nthreads,
        analysis_name=analysis_name,
        results=results_dir
    output:
        f"{results_dir}/{analysis_name}_grList.qs",
        f"{procdata}/{analysis_name}_optimal_gamma_list.qs"
    script:
        "scripts/4_selectOptimalGamma.R"


# -- 5. Build SummarizedExperiment
rule build_summarized_experiments:
    input: 
        gr_cnv=f"{procdata}/{analysis_name}_optimal_gamma_list.qs",
        pairs_file=os.path.join(metadata, pairs_file)
    params:
        nthreads=nthreads,
        analysis_name=analysis_name,
        results=results_dir
    output:
        [f"{results_dir}/{analysis_name}_{feature}_SumExp.qs" for feature 
            in ["bins", "gene"]]
    script:
        "scripts/5_buildSummarizedExperiments.R"


# -- 6. QC filter samples
rule sample_quality_control:
    input:
        procdata=procdata,
        summarized_experiments=[
            f"{results_dir}/{analysis_name}_{feature}_SumExp.qs" 
                for feature in ["bins", "gene"]
        ]
    params:
        mapd=config["mapd"],
        ndwavinesssd=config["ndwavinesssd"],
        snpqc=config["snpqc"]
    output:
        qc_csv=os.path.join(procdata, "sample_qc.csv")
    script:
        "scripts/6_sampleQualityControl.R"

# -- 7. Select top variant features (CNV regions)
feature_numbers = config["feature_numbers"]
drop_sex = config["drop_sex"]

rule select_top_variant_features:
    input:
        gr_list=f"{results_dir}/{analysis_name}_grList.qs",
        bins_sumexp=f"{results_dir}/{analysis_name}_bins_SumExp.qs",
        genes_sumexp=f"{results_dir}/{analysis_name}_gene_SumExp.qs",
    params:
        feature_numbers=feature_numbers,
        drop_sex=drop_sex
    output:
        ranked_feature_file=f"{results_dir}/{analysis_name}_features_sorted_by_mad.csv",
        feature_number_files=expand(
            "{results_dir}/{analysis_name}_{feature_number}_most_variant_{feature_type}.csv",
            results_dir=results_dir, analysis_name=analysis_name, 
            feature_number=feature_numbers, feature_type="regions"
        )
    script:
        "scripts/6_selectTopFeatures.R"