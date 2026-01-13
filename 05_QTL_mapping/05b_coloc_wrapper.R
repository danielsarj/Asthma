"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/colocalization/jobs')

# arguments
input_prefix <- c('IVA_B_17','NI_B_13','RV_B_17','IVA_CD4-T_8','NI_CD4-T_4','RV_CD4-T_5',
                  'IVA_CD8-T_1','NI_CD8-T_6','RV_CD8-T_2','IVA_Mono_1','NI_Mono_0','RV_Mono_20','IVA_NK_5','NI_NK_2','RV_NK_2')
gwas_files <- c('FerreiraMAR_COA.h.tsv.gz', 'SakaueS_COA.h.tsv.gz')

# sbtach file topper
sbatch_topper <- '#!/bin/sh\n' %&% '#SBATCH --time=20:00:00\n' %&%
  '#SBATCH --mem=50G\n' %&% '#SBATCH --partition=caslake\n' %&%
  '#SBATCH --account=pi-lbarreiro\n\n' %&% 'conda init\n' %&%
  'conda activate seurat_env\n\n' %&%
  'export R_LIBS_USER=$CONDA_PREFIX/lib/R/library\n' %&%
  'export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH\n\n' %&%
  'module load R/4.3.1'

for (f in input_prefix){
  for (g in gwas_files){
  
    # command line
      command_line <- 'Rscript /project/lbarreiro/USERS/daniel/Asthma/05_QTL_mapping/05a_coloc.R ' %&%
    '--gwas ' %&% g %&% ' --eqtl ' %&% f
    
    # create sbatch file
    cat(sbatch_topper%&%'\n\n'%&%command_line,
        file=f%&%'_'%&%g%&%'_coloc.sbatch', append=F)
    
    #submit job
    system('sbatch '%&%f%&%'_'%&%g%&%'_coloc.sbatch')
    
  }
}