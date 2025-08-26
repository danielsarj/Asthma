import os
import pandas as pd
import random

configfile: "saigeqtl_config.yaml"

# Load gene list
with open(config["gene_list"]) as f:
    GENE_CHR = {
        line.strip().split()[0]: line.strip().split()[1]
       for line in f if line.strip()
    }

# load files
GENES = list(GENE_CHR.keys())
CELLTYPE = config["celltype"]
PHENO_FILE = config["pheno_file"]
PLINK_PREFIX = config["plink_prefix"]
PLINK_IN = config["plink_in"]
CIS_REGIONS = config["cis_regions"]
N_PERMS = config["n_perms"]
PERMS = ["no_perm"] + [f"perm{i}" for i in range(1, N_PERMS + 1)]

# retrieve batch columns
df=pd.read_csv(PHENO_FILE, sep="\t", nrows=0)
columns=list(df.columns)
selected_columns = [col for col in columns if "batchID" in col]
BATCHES=",".join(selected_columns)

rule all:
    input:
        expand(
            "/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step3/outputs/{celltype}/perms/{perm}/{gene}.genePval.txt",
            celltype=CELLTYPE, gene=GENES, perm=PERMS
        )

rule permute_plink:
    input:
        bed = PLINK_IN + ".bed",
        bim = PLINK_IN + ".bim",
        fam = PLINK_IN + ".fam"
    output:
        bed = "/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/inputs/{perm}.bed",
        bim = "/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/inputs/{perm}.bim",
        fam = "/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/inputs/{perm}.fam"
    conda:
        "saigeqtl_env"
    shell:
        """
        if [ "{wildcards.perm}" = "no_perm" ]; then
            cp {input.bed} {output.bed}
            cp {input.bim} {output.bim}
            cp {input.fam} {output.fam}
        else
            awk '{{print $1, $2}}' {input.fam} | shuf | \
            awk 'NR==FNR{{a[NR]=$2; next}} {{print $1, $2, $1, a[FNR]}}' - {input.fam} \
            > /project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/inputs/perm_ids.{wildcards.perm}.txt

            plink --bfile {PLINK_IN} \
                  --update-ids /project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/inputs/perm_ids.{wildcards.perm}.txt \
                  --make-bed \
                  --out /project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/inputs/{wildcards.perm}
        fi
        """

rule step1:
    output:
        rda="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/outputs/{celltype}/perms/no_perm/{celltype}_{gene}.rda",
        variance_ratio="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/outputs/{celltype}/perms/no_perm/{celltype}_{gene}.varianceRatio.txt"
    conda:
        "saigeqtl_env"
    params:
        inv_norm="FALSE",
        covars=lambda wildcards: f"age_Scale,{BATCHES},YRI_Scale,percent.mt,PC1,PC2,PC3,PC4",
        sample_covars=lambda wildcards: f"age_Scale,YRI_Scale,{BATCHES}",
        offset_col="log_total_counts",
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
            --offsetCol={params.offset_col} \
            --traitType=count \
            --outputPrefix=/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/outputs/{wildcards.celltype}/perms/no_perm/{wildcards.celltype}_{wildcards.gene} \
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
        bed="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/inputs/{perm}.bed",
        bim="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/inputs/{perm}.bim",
        fam="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/inputs/{perm}.fam",
        model="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/outputs/{celltype}/perms/no_perm/{celltype}_{gene}.rda",
        variance_ratio="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/outputs/{celltype}/perms/no_perm/{celltype}_{gene}.varianceRatio.txt"
    output:
        assoc="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/outputs/{celltype}/perms/{perm}/{gene}.SAIGE.txt",
        index="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/outputs/{celltype}/perms/{perm}/{gene}.SAIGE.txt.index",
        ranges=temp("/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/tmp/{celltype}_{gene}_{perm}_ranges.txt")
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
            --SAIGEOutputFile={output.assoc} \
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
        assoc="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/outputs/{celltype}/perms/{perm}/{gene}.SAIGE.txt"
    output:
        geneassoc="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step3/outputs/{celltype}/perms/{perm}/{gene}.genePval.txt"
    conda:
        "saigeqtl_env"
    shell:
        """
        step3_gene_pvalue_qtl.R \
            --assocFile={input.assoc} \
            --geneName={wildcards.gene} \
            --genePval_outputFile={output.geneassoc}
        """