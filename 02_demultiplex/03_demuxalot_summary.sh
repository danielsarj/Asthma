#!/bin/sh
module load singularity

InputIDs=("B1-NI" "B1-RV" "B1-IVA" "B2-NI" "B2-RV" "B2-IVA" "B3-NI" "B3-RV" "B3-IVA" "B4-NI" "B4-RV" "B4-IVA")

for ID in ${InputIDs[@]}; do

    echo ${ID}

    DEMUXALOT_OUTDIR=/project/lbarreiro/USERS/daniel/asthma_project/alignment/${ID}_GRCh38/demuxalot

    singularity exec --bind /project/lbarreiro/USERS/daniel /project/lbarreiro/USERS/daniel/SOFTWARE/Demuxafy.sif bash demuxalot_summary.sh $DEMUXALOT_OUTDIR/assignments_refined.tsv.gz > $DEMUXALOT_OUTDIR/demuxalot_summary.tsv

done
