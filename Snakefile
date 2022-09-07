# -- 0.1 Dependeny imports
import os
import re
import glob
import pandas as pd
import numpy as np


# -- 0.2 Load configuration files
configfile: "config.yaml"

rawdata = config["rawdata_dir"]
metadata = config["metadata_dir"]
procdata = config["procdata_dir"]
pairs_file = config["pairs_file"]
analysis_name = config["analysis_name"]
results_dir = config["results_dir"]

# Load the pairs file or create one if it doesn't exist
## Note: this must be done in pipeline preamble to ensure sample names can
## be used to define output files! Should be in .tsv format.
if (os.path.exists(os.path.join(f""))):
    pairs_df = pd.read_csv(os.path.join(metadata, pairs_file), sep="\t")
else:
    if re.match("snp6", config["array_type"]):
        cel_paths = glob.glob(f"{rawdata}/**/*CEL", recursive=True)
        pairs_df = pd.DataFrame({
            "cel_files": cel_paths,
            "SampleName": [re.split(r"\/|\\", path)[1] for path in cel_paths]
        })
        pairs_df.to_csv(os.path.join(metadata, pairs_file), sep="\t")
    else:
        except OSError:
            print(
                f"No pairs file found. Please create one at:\n"
                f"\t{os.path.join(metadata, pairs_file)}"
            )


reference = config["reference"]
ref_symbol = reference.split(".")[-1]
nthreads=config["nthreads"]

tcn_cutoffs = config["tcn_cutoffs"]
feature_numbers = config["feature_numbers"]
drop_sex = config["drop_sex"]
feature_col = config["feature_col"]


# -- 0.3 All rule, defines the final output of this pipeline and runs all necessary steps
rule all:
    input:
        ranked_feature_file=f"{results_dir}/{analysis_name}_features_sorted_by_mad.csv",
        feature_number_files=expand(
            "{results_dir}/{analysis_name}_{feature_number}_most_variant_{feature_type}.csv",
            results_dir=results_dir, analysis_name=analysis_name,
            feature_number=feature_numbers, feature_type="regions"
        )

# -- 1. Batch processing of raw CEL or BAM files
## FIXME:: Whether CNV or Array go in the file name is assay dependent!
rule batch_process_rawdata:
    input:
        pairs_file=os.path.join(metadata, pairs_file)
    params:
        analysis_name=config["analysis_name"],
        array_type=config["array_type"],
        cluster_type=config["cluster_type"]
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
rule estimate_copy_number:
    input:
        expand("{procdata}/{sample_name}/{segmenter}/L2R/{sample_name}.SEG.{segmenter}.RDS",
            procdata=procdata, sample_name=pairs_df.SampleName, segmenter=config["segmenter"])
    params:
        procdata=procdata,
        gamma_range=config["gamma_range"]
    threads: nthreads
    output:
        touch(f"{procdata}/estimate_copy_number.done")
    script:
        "scripts/3_estimateCopyNumber.R"


# -- 4. Select the optimal gamma value for each sample
rule select_optimal_gamma:
    input:
        f"{procdata}/estimate_copy_number.done"
    params:
        out_dir=procdata,
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
        cnv_objects=[
            *[f"{results_dir}/{analysis_name}_{feature}_SumExp.qs"
                for feature in ["bins", "gene"]],
            f"{results_dir}/{analysis_name}_grList.qs"
        ]
    params:
        procdata=procdata,
        mapd=config["mapd"],
        ndwavinesssd=config["ndwavinesssd"],
        snpqc=config["snpqc"],
        cellularity=config["cellularity"]
    output:
        qc_csv=os.path.join(procdata, "sample_qc.csv"),
        cnv_objects=[
            *[f"{results_dir}/{analysis_name}_{feature}_SumExp_passed_qc.qs"
                for feature in ["bins", "gene"]],
            f"{results_dir}/{analysis_name}_grList_pass_qc.qs"
        ]
    script:
        "scripts/6_sampleQualityControl.R"


# -- 7. Use custom log2r cut-offs for calling TCN, useful if cellularity is known
rule custom_total_copy_calls:
    input:
        gr_list=f"{results_dir}/{analysis_name}_grList_pass_qc.qs"
            if len(tcn_cutoffs) > 0 else None,
        bins_sumexp=f"{results_dir}/{analysis_name}_bins_SumExp_passed_qc.qs"
            if len(tcn_cutoffs) > 0 else None,
        genes_sumexp=f"{results_dir}/{analysis_name}_gene_SumExp_passed_qc.qs"
            if len(tcn_cutoffs) > 0 else None
    params:
        tcn_cutoffs=tcn_cutoffs
    output:
        gr_list=f"{results_dir}/{analysis_name}_grList_passed_qc_custom_tcn.qs"
            if len(tcn_cutoffs) > 0 else None,
        bins_sumexp=f"{results_dir}/{analysis_name}_bins_SumExp_passed_qc_custom_tcn.qs"
            if len(tcn_cutoffs) > 0 else None,
        genes_sumexp=f"{results_dir}/{analysis_name}_gene_SumExp_passed_qc_custom_tcn.qs"
            if len(tcn_cutoffs) > 0 else None
    script:
        "scripts/7_customTotalCopyCalls.R"


# -- 8. Select top variant features (CNV regions)
rule select_top_variant_features:
    input:
        gr_list=f"{results_dir}/{analysis_name}_grList_passed_qc.qs"
            if len(tcn_cutoffs) == 0
            else f"{results_dir}/{analysis_name}_grList_passed_qc_custom_tcn.qs",
        bins_sumexp=f"{results_dir}/{analysis_name}_bins_SumExp_passed_qc.qs"
            if len(tcn_cutoffs) == 0
            else f"{results_dir}/{analysis_name}_bins_SumExp_passed_qc_custom_tcn.qs",
        genes_sumexp=f"{results_dir}/{analysis_name}_gene_SumExp_passed_qc.qs"
            if len(tcn_cutoffs) == 0
            else f"{results_dir}/{analysis_name}_gene_SumExp_passed_qc_custom_tcn.qs"
    params:
        feature_numbers=feature_numbers,
        drop_sex=drop_sex,
        feature_col=feature_col
    output:
        ranked_feature_file=f"{results_dir}/{analysis_name}_features_sorted_by_mad.csv",
        feature_number_files=expand(
            "{results_dir}/{analysis_name}_{feature_number}_most_variant_{feature_type}.csv",
            results_dir=results_dir, analysis_name=analysis_name,
            feature_number=feature_numbers, feature_type="regions"
        )
    script:
        "scripts/8_selectTopFeatures.R"