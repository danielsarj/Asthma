"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/alignment')

# paths and arguments
project_folder <- '/project/lbarreiro/USERS/daniel/'
batches <- c('B1','B2','B3','B4')
conditions <- c('NI', 'RV', 'IVA')

# sbtach file topper
sbatch_topper <- '#!/bin/sh\n' %&% '#SBATCH --time=36:00:00\n' %&%
'#SBATCH --mem=120G\n' %&% '#SBATCH --partition=caslake\n' %&%
'#SBATCH --account=pi-lbarreiro\n\nmodule load singularity'

for (b in 1:length(batches)){
for (c in 1:length(conditions)){

# create output directory
outdir <- '/project/lbarreiro/USERS/daniel/asthma_project/alignment/'%&% batches[b] %&%'-' %&% conditions[c] %&%'_GRCh38/demuxalot'
system('mkdir ' %&% outdir)

# retrieve paths to input files
bam_file <- '/project/lbarreiro/USERS/daniel/asthma_project/alignment/'%&% batches[b] %&%'-' %&% conditions[c] %&%'_GRCh38/outs/possorted_genome_bam.bam'
barcode_file <- '/project/lbarreiro/USERS/daniel/asthma_project/alignment/'%&% batches[b] %&%'-' %&% conditions[c] %&%'_GRCh38/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
vcf_file <- '/project/lbarreiro/USERS/daniel/asthma_project/genotypes/'%&% batches[b] %&%'_genotypes.vcf'
inds <- '/project/lbarreiro/USERS/daniel/asthma_project/genotypes/'%&% batches[b] %&%'_indvs.txt'

# command line
command_line <- 'singularity exec --bind ' %&% project_folder %&%
' /project/lbarreiro/USERS/daniel/SOFTWARE/Demuxafy.sif Demuxalot.py -a ' %&%
bam_file %&% ' -n ' %&% inds %&% ' -v ' %&% vcf_file %&% ' -b ' %&% barcode_file %&%
' -o ' %&% outdir %&% ' -r True'

# create sbatch file
cat(sbatch_topper%&%'\n'%&%command_line, file=batches[b]%&%'_'%&% conditions[c]%&%'_GRCh38_demuxalot.sbatch', append=F)

#submit job
system('sbatch '%&%batches[b]%&%'_'%&%conditions[c]%&%'_GRCh38_demuxalot.sbatch')
}
}
