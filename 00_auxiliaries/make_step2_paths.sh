#!/bin/bash

# define output directory
OUTPUT_DIR="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step2/inputs"

# list of cell types to process
for CELL_TYPE in NK_NI B_NI CD4_T_NI CD8_T_NI monocytes_NI; do

    # define gene list and output directory
    GENE_LIST="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/inputs/${CELL_TYPE}_gene_list.txt"
    OUTDIR="/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step1/outputs/${CELL_TYPE}"

    # define the output file for this cell type
    OUTPUT_FILE="${OUTPUT_DIR}/${CELL_TYPE}_gene_paths.tsv"
    > "$OUTPUT_FILE"  # clear it if it exists

    # process genes
    while IFS= read -r gene; do
        rda_path="${OUTDIR}/${CELL_TYPE}_${gene}.rda"
        varratio_path="${OUTDIR}/${CELL_TYPE}_${gene}.varianceRatio.txt"
        echo -e "${gene}\t${rda_path}\t${varratio_path}" >> "$OUTPUT_FILE"
    done < "$GENE_LIST"
done
