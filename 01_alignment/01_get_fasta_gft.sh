cd /project/lbarreiro/USERS/daniel/asthma_project/alignment

# Human
cp /project/lbarreiro/USERS/cecily/REFERENCES/GRCh38.subset.primary_assembly.genome.fa .
mv GRCh38.subset.primary_assembly.genome.fa GRCh38.subset.primary_assembly.genome.fna
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz

# H1N1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/343/785/GCF_001343785.1_ViralMultiSegProj274766/GCF_001343785.1_ViralMultiSegProj274766_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/343/785/GCF_001343785.1_ViralMultiSegProj274766/GCF_001343785.1_ViralMultiSegProj274766_genomic.gtf.gz

# RV
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/816/835/GCF_002816835.1_ASM281683v1/GCF_002816835.1_ASM281683v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/816/835/GCF_002816835.1_ASM281683v1/GCF_002816835.1_ASM281683v1_genomic.gtf.gz

# unzip them
gunzip *.gz

# filter human GTF
sed -i ’s/chrM/chrMT/g’ gencode.v47.basic.annotation.gtf

# concatenate all FASTAs and GTFs
cat *.fna > merged_reference.fna
cat *.gtf > merged_reference.gtf
