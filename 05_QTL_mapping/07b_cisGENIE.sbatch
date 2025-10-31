#!/bin/sh
#SBATCH --time=5:00:00
#SBATCH --mem=50G
#SBATCH --partition=caslake
#SBATCH --account=pi-lbarreiro

# load plink and change directory
module load plink
cd /project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/GENIE/

# load variables
CONDITION=$1
CELLTYPE=$2
GENE_LOCAL="inputs/${CONDITION}_${CELLTYPE}_gene_annotation.txt"
GENIE_BIN="/project/lbarreiro/USERS/daniel/SOFTWARE/GENIE/build/GENIE"

# read file line by line, skipping header
tail -n +2 "$GENE_LOCAL"  | head -n 10 |  while IFS=$'\t' read -r GENE_NAME CHROMOSOME START END _rest; do

    # calculate +/- 500kb window
    FROM=$(( START - 500000 ))
    TO=$(( END + 500000 ))

    # subset PLINK file for each gene (+/- 500kb window)
    plink --bfile filtered_merged --chr "$CHROMOSOME" --from-bp "$FROM" --to-bp "$TO" --make-bed --out "GENIE_${CONDITION}_${CELLTYPE}_cis_${GENE_NAME}"

    # count number of SNPs in the .bim file and create annotation file
    NUM_SNPS=$(wc -l < "GENIE_${CONDITION}_${CELLTYPE}_cis_${GENE_NAME}.bim")
    yes 1 | head -n "$NUM_SNPS" > "inputs/${CONDITION}_${CELLTYPE}_${GENE_NAME}_anno.txt"

    # run cisGENIE for asthma
    $GENIE_BIN -g GENIE_${CONDITION}_${CELLTYPE}_cis_${GENE_NAME} -p inputs/${CONDITION}_${CELLTYPE}_${GENE_NAME}.txt -a inputs/${CONDITION}_${CELLTYPE}_${GENE_NAME}_anno.txt -c inputs/COV_albuterol.txt -o outputs/${CONDITION}_${CELLTYPE}_${GENE_NAME}_asthma.txt -e inputs/ENV_asthma.txt -m G+GxE -k 5 -jn 10

    # run cisGENIE for income
    $GENIE_BIN -g GENIE_${CONDITION}_${CELLTYPE}_cis_${GENE_NAME} -p inputs/${CONDITION}_${CELLTYPE}_${GENE_NAME}.txt -a inputs/${CONDITION}_${CELLTYPE}_${GENE_NAME}_anno.txt -o outputs/${CONDITION}_${CELLTYPE}_${GENE_NAME}_income.txt -e inputs/ENV_income.txt -m G+GxE -k 5 -jn 10

    # remove temporary files
    rm "GENIE_${CONDITION}_${CELLTYPE}_cis_${GENE_NAME}".*

done