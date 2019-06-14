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


seq_path='/scratch/rx32940/lepto_wgs_seq'
o_path='/scratch/rx32940/bwa_results'

# index the reference
bwa index $seq_path/Lai_56601.fasta

for file in $seq_path/*.fasta; do
    echo "in loop"
    echo $file
    isolate=$(echo $file | awk -F'[/.]' '{print $4}')
    echo $isolate
    if [ "$isolate" != "Lai_56601" ]; then
        echo "in if"
        # align the contigs to reference 
        bwa mem -t2 $seq_path/Lai_56601.fasta $seq_path/$isolate.fasta > $o_path/$isolate.sam
        # -o output, convert to binary, bam, format
        samtools sort -o $o_path/$isolate.bam $o_path/$isolate.sam
        #index the alignment file
        samtools index $o_path/$isolate.bam
        # producing genotype likelihoods in VCF or BCF format
        bcftools mpileup -Ou -f $seq_path/Lai_56601.fasta $o_path/$isolate.bam > $o_path/$isolate.bcf
        bcftools call -mv -Ob $o_path/$isolate.bcf > $o_path/${isolate}_final.bcf
        bcftools view -i '%QUAL>=20' $o_path/${isolate}_final.bcf > $o_path/${isolate}_final.vcf
    fi
done



