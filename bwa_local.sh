#!/bin/bash

cd "$(dirname "$0")"

source activate bwa #activate bwa env.

path='bwa_results'# lead all output file to this folder

# index the reference
bwa index Lai_56601.fasta

# align the contigs to reference 
echo $path
bwa mem -t2 Lai_56601.fasta JQPC01.fasta > $path/JQPC01.sam
# -o output, convert to binary, bam, format
samtools sort -o $path/JQPC01.bam $path/JQPC01.sam
#index the alignment file
samtools index $path/JQPC01.bam
# producing genotype likelihoods in VCF or BCF format
bcftools mpileup -Ou -f Lai_56601.fasta $path/JQPC01.bam > $path/JQPC01.bcf
bcftools call -mv -Ob $path/JQPC01.bcf > $path/JQPC01_final.bcf
bcftools view -i '%QUAL>=20' $path/JQPC01_final.bcf > $path/JQPC01_final.vcf




