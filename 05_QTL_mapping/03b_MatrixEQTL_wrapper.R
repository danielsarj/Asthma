"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/jobs')

# arguments
conditions <- c('NI', 'RV', 'IVA')
celltypes <- c('B', 'T-CD4', 'T-CD8', 'Mono', 'NK')
pcs <- c(1:20)

# sbtach file topper
sbatch_topper <- '#!/bin/sh\n' %&% '#SBATCH --time=2:00:00\n' %&%
  '#SBATCH --mem=50G\n' %&% '#SBATCH --partition=caslake\n' %&%
  '#SBATCH --account=pi-lbarreiro\n\n' %&% 'conda init\n' %&%
  'conda activate seurat_env\n\n' %&%
  'export R_LIBS_USER=$CONDA_PREFIX/lib/R/library\n' %&%
  'export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH\n\n' %&%
  'module load R/4.3.1'

for (cond in conditions){
  for (ctype in celltypes){
    for (pc in pcs){
      # command line
        command_line <- 'Rscript /project/lbarreiro/USERS/daniel/Asthma/05_QTL_mapping/03a_runMatrixEQTL.R ' %&%
      '--cond ' %&% cond %&% ' --ctype ' %&% ctype %&% ' --pcs ' %&% pc
    
      # create sbatch file
      cat(sbatch_topper%&%'\n\n'%&%command_line,
          file=cond%&%'_'%&%ctype%&%'_'%&%pc%&%'_MatrixEQTL.sbatch', append=F)
    
      #submit job
      system('sbatch '%&%cond%&%'_'%&%ctype%&%'_'%&%pc%&%'_MatrixEQTL.sbatch')
    }
  }
}