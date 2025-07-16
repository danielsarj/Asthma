conda activate saigeqtl_env
snakemake -s saigeqtl.smk \
  --jobs 1000 \
  --latency-wait 30 \
  --keep-going \
  --rerun-incomplete \
  --use-envmodules \
  --cluster "sbatch \
    --partition=caslake \
    --account=pi-lbarreiro \
    --time=1:00:00 \
    --mem=32G \
    --cpus-per-task=1 \
    -e /project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/slurm_logs/{rule}.{wildcards}.err \
    -o /project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/slurm_logs/{rule}.{wildcards}.out"