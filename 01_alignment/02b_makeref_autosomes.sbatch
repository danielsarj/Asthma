#!/bin/sh
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --partition=caslake
#SBATCH --account=pi-lbarreiro

# make cellranger ref
cd /project/lbarreiro/USERS/daniel/asthma_project/alignment
module load cellranger
cellranger mkref --genome=GRCh38 --fasta=GRCh38.autosomes.fna --genes=GRCh38.autosomes.gtf
