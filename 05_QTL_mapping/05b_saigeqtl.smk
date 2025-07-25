import os
import pandas as pd
import re

configfile: "/project/lbarreiro/USERS/daniel/Asthma/05_QTL_mapping/05c_saigeqtl_config.yaml"

# Load gene list
with open(config["gene_list"]) as f:
    GENE_CHR = {
        line.strip().split()[0]: line.strip().split()[1]
        for line in f if line.strip()
    }
GENES = list(GENE_CHR.keys())
CELLTYPE = config["celltype"]
PHENO_FILE = config["pheno_file"]
PLINK_PREFIX = config["plink_prefix"]
PLINK_IN = config["plink_in"]

# Retrieve covariates
df=pd.read_csv(PHENO_FILE, sep="\t", nrows=0)
columns=list(df.columns)
# Find index of the last "PC#" column
last_pc_idx=max(i for i, col in enumerate(columns) if re.match(r"PC\d+$", col))
# Get column names from 3rd column up to and including last PC
selected_columns=columns[2:last_pc_idx + 1]
# Join as comma-separated string
CELL_COVS=",".join(selected_columns)

rule all:
    input:
        expand("/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/SAIGE_results/STEP3/{celltype}_{gene}.genePval.txt", celltype=CELLTYPE, gene=GENES)

rule step1:
    output:
        rda="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/SAIGE_results/STEP1/{celltype}_{gene}.rda",
        variance_ratio="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/SAIGE_results/STEP1/{celltype}_{gene}.varianceRatio.txt"
    conda:
        "saigeqtl_env"
    params:
        sample_covars="V1,V2,V3,V4,age,gender", 
        sample_id_col="IDs",
        cell_id_col="cell_ID"
    shell:
        """
        step1_fitNULLGLMM_qtl.R \
            --useSparseGRMtoFitNULL=FALSE \
            --useGRMtoFitNULL=FALSE \
            --phenoFile={PHENO_FILE} \
            --phenoCol={wildcards.gene} \
            --covarColList={CELL_COVS} \
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
            --plinkFile={PLINK_PREFIX} \
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
        chr=lambda wildcards: GENE_CHR[wildcards.gene]
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