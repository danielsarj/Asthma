#!/bin/sh
#SBATCH --time=36:00:00
#SBATCH --mem=100G
#SBATCH --partition=caslake
#SBATCH --account=pi-lbarreiro

module load plink
cd  /project/lbarreiro/USERS/daniel/asthma_project/genotypes/imputed_vcfs

# bgzip and tabix imputated vcfs 
for chr in {1..22}; do
       gunzip chr${chr}.dose.vcf.gz
      /project/lbarreiro/USERS/daniel/SOFTWARE/htslib-1.19/bgzip chr${chr}.dose.vcf
       /project/lbarreiro/USERS/daniel/SOFTWARE/htslib-1.19/tabix -p vcf chr${chr}.dose.vcf.gz
done

# concatenate vcfs
/project/lbarreiro/USERS/daniel/SOFTWARE/bcftools-1.19/bcftools concat chr*.dose.vcf.gz -Oz -o merged.vcf.gz 
/project/lbarreiro/USERS/daniel/SOFTWARE/htslib-1.19/tabix -p vcf merged.vcf.gz

# remove indels and multi allelic sites
/project/lbarreiro/USERS/daniel/SOFTWARE/bcftools-1.19/bcftools view -m2 -M2 -v snps merged.vcf.gz -Oz -o reduced_merged.vcf.gz

# convert vcf to plink format
plink --chr 1-22 --keep-allele-order --vcf-idspace-to - --make-bed --out merged --vcf reduced_merged.vcf.gz

# LD prune
plink --bfile merged --indep-pairwise 500 50 0.5 --out merged
plink --bfile merged --extract merged.prune.in --make-bed --out merged_pruned

# perform genotype QC 
plink --bfile merged_pruned --maf 0.05 --hwe 1e-6 --make-bed --out filtered_merged

# change snp IDs
awk '{print $1, $1":"$4, $3, $4, $5, $6}' filtered_merged.bim > updated.bim
mv updated.bim filtered_merged.bim

# create dosage files per chr
for chr in {1..22}; do
	plink --bfile filtered_merged --chr ${chr} --out imputed_chr${chr}_dosage --recode A
done
