module load plink
cd  /project/lbarreiro/USERS/daniel/asthma_project/genotypes

# convert vcf to plink format
plink --allow-extra-chr --chr 1-22 --keep-allele-order --vcf-idspace-to - --make-bed --out all_genotypes --vcf all_genotypes.GRCh38.vcf.gz

# perform genotype QC 
plink --bfile all_genotypes --geno 0.1 --maf 0.05 --hwe 1e-6 --make-bed --out filtered_genotypes

# back to vcf
plink --bfile filtered_genotypes --keep-allele-order --recode vcf --out filtered_genotypes

# rename chr 
/project/lbarreiro/USERS/daniel/SOFTWARE/bcftools-1.19/bcftools annotate --rename-chrs chr_names.txt filtered_genotypes.vcf -Oz -o filtr_renamed_genotypes.vcf
