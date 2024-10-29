cd /project/lbarreiro/USERS/daniel/asthma_project/alignment

# Human
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
cp /project/lbarreiro/USERS/cecily/REFERENCES/GRCh38.subset.primary_assembly.genome.fa .
mv GRCh38.subset.primary_assembly.genome.fa GRCh38.subset.primary_assembly.genome.fna
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz

# H1N1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/343/785/GCF_001343785.1_ViralMultiSegProj274766/GCF_001343785.1_ViralMultiSegProj274766_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/343/785/GCF_001343785.1_ViralMultiSegProj274766/GCF_001343785.1_ViralMultiSegProj274766_genomic.gtf.gz

# RV
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/816/835/GCF_002816835.1_ASM281683v1/GCF_002816835.1_ASM281683v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/816/835/GCF_002816835.1_ASM281683v1/GCF_002816835.1_ASM281683v1_genomic.gtf.gz

# unzip them
gunzip *.gz

# rename human contigs in GTF
sed -i ’s/NC_000001.11/chr1/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000002.12/chr2/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000003.12/chr3/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000004.12/chr4/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000005.10/chr5/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000006.12/chr6/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000007.14/chr7/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000008.11/chr8/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000009.12/chr9/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000010.11/chr10/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000011.10/chr11/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000012.12/chr12/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000013.11/chr13/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000014.9/chr14/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000015.10/chr15/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000016.10/chr16/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000017.11/chr17/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000018.10/chr18/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000019.10/chr19/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000020.11/chr20/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000021.9/chr21/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000022.11/chr22/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000023.11/chrX/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_000024.10/chrY/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
sed -i ’s/NC_012920.1/chrMT/g’ GCF_000001405.40_GRCh38.p14_genomic.gtf
grep -v -E "^chr[0-9]+_" GCF_000001405.40_GRCh38.p14_genomic.gtf > GRCh38.cannonical_chrs.gtf
mv GCF_000001405.40_GRCh38.p14_genomic.gtf GCF_000001405.40_GRCh38.p14_genomic.gtfn 

# concatenate all FASTAs and GTFs
cat *.fna > merged_reference.fna
cat *.gtf > merged_reference.gtf

