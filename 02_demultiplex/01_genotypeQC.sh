module load plink

# convert vcf to plink format
plink --allow-extra-chr --chr 1-22 --keep-allele-order --vcf-idspace-to - --make-bed --out /project/lbarreiro/USERS/daniel/asthma_project/genotypes/all_genotypes --vcf /project/lbarreiro/USERS/daniel/asthma_project/genotypes/all_genotypes.GRCh38.vcf.gz

# perform genotype QC 
plink --bfile /project/lbarreiro/USERS/daniel/asthma_project/genotypes/all_genotypes --geno 0.1 --maf 0.05 --hwe 1e-6 --make-bed --out /project/lbarreiro/USERS/daniel/asthma_project/genotypes/filtered_genotypes

# back to vcf
plink --bfile /project/lbarreiro/USERS/daniel/asthma_project/genotypes/filtered_genotypes --keep-allele-order --recode vcf --out /project/lbarreiro/USERS/daniel/asthma_project/genotypes/filtered_genotypes

gzip project/lbarreiro/USERS/daniel/asthma_project/genotypes/*.vcf
