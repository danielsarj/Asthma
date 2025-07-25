import os
import pandas as pd
import re

configfile: "/project/lbarreiro/USERS/daniel/Asthma/05_QTL_mapping/05c_saigeqtl_config.yaml"

# Retrieve PLINK files from the config file
PLINK_PREFIX=config["plink_prefix"]
PLINK_IN=config["plink_in"]

# Create dictionary for input files per celltype
with open(config["input_list"]) as f:
    INPUT_LIST={
        line.strip().split()[0]: {
            "gene_list": line.strip().split()[1],
            "pheno_file": line.strip().split()[2],
        }
        for line in f if line.strip()}

# Retrieve covariates per celltype
covs_celltype={}
for celltype, paths in INPUT_LIST.items():
    pheno_path=paths["pheno_file"]
    # Read header only 
    df=pd.read_csv(pheno_path, sep="\t", nrows=0)
    columns=list(df.columns)
    # Find index of the last "PC#" column
    last_pc_idx=max(i for i, col in enumerate(columns) if re.match(r"PC\d+$", col))
    # Get column names from 3rd column up to and including last PC
    selected_columns=columns[2:last_pc_idx + 1]
    # Join as comma-separated string
    pc_string=",".join(selected_columns)
    covs_celltype[celltype]=pc_string

# Get gene-chr list for each celltype
gene_chr_by_celltype={}
for celltype, paths in INPUT_LIST.items():
    gene_file=paths["gene_list"]
    # Read file without headers
    df=pd.read_csv(gene_file, sep="\t", header=None)
    # Assign column: first=gene, second=chromosome
    df.columns=["gene", "chromosome"]
    # Create dictionary: gene_name -> chromosome
    gene_chr_dict=pd.Series(df["chromosome"].values, index=df["gene"]).to_dict()
    gene_chr_by_celltype[celltype]=gene_chr_dict

cell_gene_pairs=[
    (celltype, gene)
    for celltype, gene_dict in gene_chr_by_celltype.items()
    for gene in gene_dict
]

rule all:
    input:
        expand(
            "/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/SAIGE_results/STEP3/{celltype}_{gene}.genePval.txt",
            celltype=[ct for ct, _ in cell_gene_pairs],
            gene=[g for _, g in cell_gene_pairs]
        )

rule step1:
    output:
        rda="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/SAIGE_results/STEP1/{celltype}_{gene}.rda",
        variance_ratio="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/SAIGE_results/STEP1/{celltype}_{gene}.varianceRatio.txt"
    conda:
        "saigeqtl_env"
    params:
        covars=lambda wildcards: covs_celltype[wildcards.celltype],
        sample_covars="V1,V2,V3,V4,age,gender", 
        sample_id_col="IDs",
        cell_id_col="cell_ID",
        pheno_file=lambda wildcards: INPUT_LIST[wildcards.celltype]["pheno_file"],
        plink_prefix=PLINK_PREFIX 
    shell:
        """
        step1_fitNULLGLMM_qtl.R \
            --useSparseGRMtoFitNULL=FALSE \
            --useGRMtoFitNULL=FALSE \
            --phenoFile={params.pheno_file} \
            --phenoCol={wildcards.gene} \
            --covarColList={params.covars} \
            --sampleCovarColList={params.sample_covars} \
            --sampleIDColinphenoFile={params.sample_id_col} \
            --cellIDColinphenoFile={params.cell_id_col} \
            --traitType=count \
            --outputPrefix=/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/SAIGE_results/STEP1/{wildcards.celltype}_{wildcards.gene} \
            --skipVarianceRatioEstimation=FALSE \
            --isRemoveZerosinPheno=FALSE \
            --isCovariateOffset=FALSE \
            --isCovariateTransform=TRUE \
            --skipModelFitting=FALSE \
            --tol=0.00001 \
            --plinkFile={params.plink_prefix} \
            --IsOverwriteVarianceRatioFile=TRUE
        """

rule step2:
    input:
        bed=PLINK_IN + ".bed",
        bim=PLINK_IN + ".bim",
        fam=PLINK_IN + ".fam",
        model="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/SAIGE_results/STEP1/{celltype}_{gene}.rda",
        variance_ratio="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/SAIGE_results/STEP1/{celltype}_{gene}.varianceRatio.txt"
    output:
        assoc="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/SAIGE_results/STEP2/{celltype}_{gene}.SAIGE.txt",
        index="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/SAIGE_results/STEP2/{celltype}_{gene}.SAIGE.txt.index"
    conda:
        "saigeqtl_env"
    params:
        minMAF="0.05",
        cutoff="2",
        chunk_size="10000",
        loco="FALSE",
        chr=lambda wildcards: gene_chr_by_celltype[wildcards.celltype][wildcards.gene]
    shell:
        """
        step2_tests_qtl.R \
            --bedFile={input.bed} \
            --bimFile={input.bim} \
            --famFile={input.fam} \
            --SAIGEOutputFile=/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/SAIGE_results/STEP2/{wildcards.celltype}_{wildcards.gene}.SAIGE.txt \
            --chrom={params.chr} \
            --minMAF={params.minMAF} \
            --GMMATmodelFile={input.model} \
            --SPAcutoff={params.cutoff} \
            --LOCO={params.loco} \
            --varianceRatioFile={input.variance_ratio} \
            --markers_per_chunk={params.chunk_size}
        """

rule step3:
    input:
        assoc="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/SAIGE_results/STEP2/{celltype}_{gene}.SAIGE.txt"
    output:
        geneassoc="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/SAIGE_results/STEP3/{celltype}_{gene}.genePval.txt"
    conda:
        "saigeqtl_env"
    shell:
        """
        step3_gene_pvalue_qtl.R \
            --assocFile={input.assoc} \
            --geneName={wildcards.gene} \
            --genePval_outputFile=/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/SAIGE_results/STEP3/{wildcards.celltype}_{wildcards.gene}.genePval.txt
        """