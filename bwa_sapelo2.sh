#!/bin/bash
#PBS -q batch                                                            
#PBS -N lepto_snps                                            
#PBS -l nodes=1:ppn=2 -l mem=20gb                                        
#PBS -l walltime=10:00:00                                                
#PBS -M rx32940@uga.edu                                                  
#PBS -m abe                                                              
#PBS -o /scratch/rx32940/bwa_results                        
#PBS -e /scratch/rx32940/bwa_results                        
#PBS -j oe     


cd $PBS_O_WORKDIR

module load BWA/0.7.15-foss-2016b
module load SAMtools/1.9-foss-2016b
module load BCFtools/1.9-foss-2016b 

seq_path='/scratch/rx32940/lepto_wgs_seq'# lead all output file to this folder
echo $seq_path
echo "seq"
o_path='/scratch/rx32940/bwa_results'
echo $o_path
echo "o"

for file in $seq_path; do
    echo "in loop"
    echo $file
    if [$file != "Lai_56601.fasta"]
    then
        echo "in if"
        # index the reference
        echo $seq
        bwa index $seq_path/Lai_56601.fasta
        # align the contigs to reference 
        bwa mem -t2 $seq_path/Lai_56601.fasta $seq_path/$file.fasta > $o_path/$file.sam
        # -o output, convert to binary, bam, format
        samtools sort -o $o_path/$file.bam $o_path/$file.sam
        #index the alignment file
        samtools index $o_path/$file.bam
        # producing genotype likelihoods in VCF or BCF format
        bcftools mpileup -Ou -f $seq_path/Lai_56601.fasta $o_path/$file.bam > $o_path/$file.bcf
        bcftools call -mv -Ob $o_path/$file.bcf > $o_path/${file}_final.bcf
        bcftools view -i '%QUAL>=20' $o_path/${file}_final.bcf > $o_path/${file}_final.vcf
    fi
done



