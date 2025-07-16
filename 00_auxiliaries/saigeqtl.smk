import os

configfile: "saigeqtl_config.yaml"

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

rule all:
    input:
        expand("/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/outputs/{celltype}/{gene}.SAIGE.txt", celltype=CELLTYPE, gene=GENES)

rule step1:
    output:
        rda="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/outputs/{celltype}/{celltype}_{gene}.rda",
        variance_ratio="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/outputs/{celltype}/{celltype}_{gene}.varianceRatio.txt"
    conda:
        "saigeqtl_env"
    params:
        inv_norm="FALSE",
        covars="age_Scale,YRI_Scale,PC1,PC2,PC3",
        sample_covars="age_Scale,YRI_Scale",
        sample_id_col="SOC_indiv_ID",
        cell_id_col="cell_ID"
    shell:
        """
        step1_fitNULLGLMM_qtl.R \
            --useSparseGRMtoFitNULL=FALSE \
            --useGRMtoFitNULL=FALSE \
            --invNormalize={params.inv_norm} \
            --phenoFile={PHENO_FILE} \
            --phenoCol={wildcards.gene} \
            --covarColList={params.covars} \
            --sampleCovarColList={params.sample_covars} \
            --sampleIDColinphenoFile={params.sample_id_col} \
            --cellIDColinphenoFile={params.cell_id_col} \
            --traitType=count \
            --outputPrefix=/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/outputs/{wildcards.celltype}/{wildcards.celltype}_{wildcards.gene} \
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
        bed= PLINK_IN + ".bed",
        bim= PLINK_IN + ".bim",
        fam= PLINK_IN + ".fam",
        model="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/outputs/{celltype}/{celltype}_{gene}.rda",
        variance_ratio="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/outputs/{celltype}/{celltype}_{gene}.varianceRatio.txt"
    output:
        assoc="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/outputs/{celltype}/{gene}.SAIGE.txt",
        index="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/outputs/{celltype}/{gene}.SAIGE.txt.index"
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
            --SAIGEOutputFile=/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/outputs/{wildcards.celltype}/{wildcards.gene}.SAIGE.txt \
            --chrom={params.chr} \
            --minMAF={params.minMAF} \
            --GMMATmodelFile={input.model} \
            --SPAcutoff={params.cutoff} \
            --LOCO={params.loco} \
            --varianceRatioFile={input.variance_ratio} \
            --markers_per_chunk={params.chunk_size}
        """
