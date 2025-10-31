#!/bin/bash
# define conditions and celltypes
CONDITIONS=("NI" "IVA" "RV")
CELLTYPES=("B" "T-CD4" "T-CD8" "Mono" "NK")

# path to sbatch script
SBATCH_SCRIPT="/project/lbarreiro/USERS/daniel/Asthma/05_QTL_mapping/07b_cisGENIE.sbatch"

# loop over conditions and celltypes
for CONDITION in "${CONDITIONS[@]}"; do
    for CELLTYPE in "${CELLTYPES[@]}"; do
        echo "Submitting sbatch for $CONDITION $CELLTYPE"
        sbatch "$SBATCH_SCRIPT" "$CONDITION" "$CELLTYPE"
    done
done