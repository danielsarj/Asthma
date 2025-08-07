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
CIS_REGIONS = config["cis_regions"]


rule all:
    input:
        expand("/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step3/outputs/{celltype}/exp_PCs1/{gene}.genePval.txt", celltype=CELLTYPE, gene=GENES)

rule step1:
    output:
        rda="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/outputs/{celltype}/exp_PCs1/{celltype}_{gene}.rda",
        variance_ratio="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/outputs/{celltype}/exp_PCs1/{celltype}_{gene}.varianceRatio.txt"
    conda:
        "saigeqtl_env"
    params:
        inv_norm="FALSE",
        covars="age_Scale,YRI_Scale,PC1",
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
            --outputPrefix=/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/outputs/{wildcards.celltype}/exp_PCs1/{wildcards.celltype}_{wildcards.gene} \
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
        bed=PLINK_IN+".bed",
        bim=PLINK_IN+".bim",
        fam=PLINK_IN+".fam",
        model="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/outputs/{celltype}/exp_PCs1/{celltype}_{gene}.rda",
        variance_ratio="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/outputs/{celltype}/exp_PCs1/{celltype}_{gene}.varianceRatio.txt"
    output:
        assoc="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/outputs/{celltype}/exp_PCs1/{gene}.SAIGE.txt",
        index="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/outputs/{celltype}/exp_PCs1/{gene}.SAIGE.txt.index",
        ranges=temp("/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/tmp/{celltype}_{gene}_ranges.txt")
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
        awk -v gene="{wildcards.gene}" '$1 == gene {{print $2, $3, $4}}' OFS='\t' {CIS_REGIONS} > {output.ranges}

        step2_tests_qtl.R \
            --bedFile={input.bed} \
            --bimFile={input.bim} \
            --famFile={input.fam} \
            --SAIGEOutputFile=/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/outputs/{wildcards.celltype}/exp_PCs1/{wildcards.gene}.SAIGE.txt \
            --chrom={params.chr} \
            --rangestoIncludeFile={output.ranges} \
            --minMAF={params.minMAF} \
            --GMMATmodelFile={input.model} \
            --SPAcutoff={params.cutoff} \
            --LOCO={params.loco} \
            --varianceRatioFile={input.variance_ratio} \
            --markers_per_chunk={params.chunk_size}
        """

rule step3:
    input:
        assoc="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/outputs/{celltype}/exp_PCs1/{gene}.SAIGE.txt"
    output:
        geneassoc="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step3/outputs/{celltype}/exp_PCs1/{gene}.genePval.txt"
    conda:
        "saigeqtl_env"
    shell:
        """
        step3_gene_pvalue_qtl.R \
            --assocFile={input.assoc} \
            --geneName={wildcards.gene} \
            --genePval_outputFile=/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step3/outputs/{wildcards.celltype}/exp_PCs1/{wildcards.gene}.genePval.txt
        """