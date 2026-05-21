#!/bin/bash
#SBATCH --partition=caslake
#SBATCH --account=pi-lbarreiro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --mem=10G
#SBATCH --job-name=saige_qtl
#SBATCH --output=saige_qtl.out
#SBATCH --error=saige_qtl.err

cd /project/lbarreiro/USERS/daniel/Asthma/00_auxiliaries
conda activate saigeqtl_env

snakemake -s saigeqtl.smk \
  --jobs 1000 \
  --latency-wait 30 \
  --keep-going \
  --use-conda \
  --rerun-incomplete \
  --use-envmodules \
  --cluster "sbatch \
    --partition=caslake \
    --account=pi-lbarreiro \
    --time=3:00:00 \
    --mem=32G \
    --cpus-per-task=1 \
    -e /project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/slurm_logs/{rule}.{wildcards}.err \
    -o /project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/slurm_logs/{rule}.{wildcards}.out"