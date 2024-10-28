"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/alignment')

# arguments
batches <- c('B1','B2','B3','B4')
conditions <- c('NI', 'RV', 'IVA')

# sbtach file topper
sbatch_topper <- '#!/bin/sh\n' %&% '#SBATCH --time=36:00:00\n' %&% 
  '#SBATCH --mem=180G\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=20\n' %&% '#SBATCH --partition=caslake\n' %&% 
  '#SBATCH --account=pi-lbarreiro\n'

for (b in 1:length(batches)){
  for (c in 1:length(conditions)){
    # command line
    command_line <- '/project/lbarreiro/SHARED/PROGRAMS/cellranger-7.0.1/cellranger count ' %&%
      '--id='%&% batches[b] %&%'-' %&% conditions[c] %&%
      ' --fastqs=/project/lbarreiro/USERS/jbatista/asthma_project/20240206_LH00315_0085_A22H7FKLT3-LB-MP-2s-Astma-pl1-4lns/FastX/,'%&%
      '/project/lbarreiro/USERS/jbatista/asthma_project/20240213_LH00315_0092_A22HGWMLT3-LB-MP-2s-Astma-pl1-8lns/FastX '%&%
      '--sample=LB-MP-pl1-Asthma-'%&% batches[b] %&%'-'%&% conditions[c] %&%
      ' --transcriptome=/project/lbarreiro/USERS/daniel/asthma_project/alignment/GRCh38_with_IAV_and_RV'

    # create sbatch file
    cat(sbatch_topper%&%'\n\n'%&%command_line, 
      file=batches[b]%&%'_'%&% conditions[c]%&%'_cellrangercount.sbatch', append=F)

    #submit job
    system('sbatch '%&%batches[b]%&%'_'%&%conditions[c]%&%'_cellrangercount.sbatch')
  }
}
