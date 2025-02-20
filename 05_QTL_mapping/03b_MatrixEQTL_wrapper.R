"%&%" <- function(a,b) paste(a,b, sep = "")

# arguments
conditions <- c('NI', 'RV', 'IVA')
celltypes <- c('B', 'CD4-T', 'CD8-T', 'Mono', 'NK')
chromosomes <- seq(1:22)

# sbtach file topper
sbatch_topper <- '#!/bin/sh\n' %&% '#SBATCH --time=36:00:00\n' %&%
  '#SBATCH --mem=100G\n' %&% '#SBATCH --partition=caslake\n' %&%
  '#SBATCH --account=pi-lbarreiro\n\n' %&% 'conda init\n' %&%
  'conda activate seurat_env\n\n' %&%
  'export R_LIBS_USER=$CONDA_PREFIX/lib/R/library\n' %&%
  'export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH\n\n' %&%
  'module load R/4.3.1'

for (cond in conditions){
  for (ctype in celltypes){
    for (chrom in chromosomes){
      # command line
      command_line <- 'Rscript /project/lbarreiro/USERS/daniel/Asthma/05_QTL_mapping/03a_runMatrixEQTL.R ' %&%
        '--cond ' %&% cond %&% ' --ctype ' %&% ctype %&% ' --chrom ' %&% chrom
    
      # create sbatch file
      cat(sbatch_topper%&%'\n\n'%&%command_line,
          file=cond%&%'_'%&%ctype%&%'_'%&%chrom%&%'_MatrixEQTL.sbatch', append=F)
    
      #submit job
      system('sbatch '%&%cond%&%'_'%&%ctype%&%'_'%&%chrom%&%'_MatrixEQTL.sbatch')
    }
  }
}