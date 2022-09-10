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
array_type=config["array_type"]

tcn_cutoffs = config["tcn_cutoffs"]
feature_col = config["feature_col"]


# -- 0.3 All rule, defines the final output of this pipeline and runs all necessary steps
rule all:
    input:
        f"{results_dir}/{analysis_name}_RagExp_pass_qc_custom_tcn.qs" if
            len(tcn_cutoffs) else f"{results_dir}/{analysis_name}_RagExp_pass_qc.qs"


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
        f"{procdata}/{analysis_name}_optimal_gamma.csv"
    script:
        "scripts/4_selectOptimalGamma.R"


# -- 5. Build RaggedExperiment
rule build_ragged_experiment:
    input:
        f"{procdata}/{analysis_name}_optimal_gamma.csv"
    params:
        out_dir=procdata,
        nthreads=nthreads,
        analysis_name=analysis_name,
        results=results_dir
    output:
        f"{results_dir}/{analysis_name}_RagExp.qs"
    script:
        "scripts/5_buildRaggedExperiment.R"


# -- 6. QC filter samples
rule sample_quality_control:
    input:
        ragged_exp=f"{results_dir}/{analysis_name}_RagExp.qs"
    params:
        procdata=procdata,
        nthreads=nthreads,
        mapd=config["mapd"],
        ndwavinesssd=config["ndwavinesssd"],
        snpqc=config["snpqc"],
        cellularity=config["cellularity"]
    output:
        qc_csv=os.path.join(procdata, "sample_qc.csv"),
        ragged_exp=f"{results_dir}/{analysis_name}_RagExp_pass_qc.qs"
    script:
        "scripts/6_sampleQualityControl.R"


# -- 7. Use custom log2r cut-offs for calling TCN, useful if cellularity is known
rule custom_total_copy_calls:
    input:
        f"{results_dir}/{analysis_name}_RagExp_pass_qc.qs"
    params:
        tcn_cutoffs=tcn_cutoffs
    output:
        f"{results_dir}/{analysis_name}_RagExp_pass_qc_custom_tcn.qs" if
            len(tcn_cutoffs) else f"{results_dir}/{analysis_name}_RagExp_pass_qc.qs"
    script:
        "scripts/7_customTotalCopyCalls.R"
