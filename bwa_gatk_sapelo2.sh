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
module load picard/2.16.0-Java-1.8.0_144 
module load GATK/4.0.11.0-foss-2016b-Python-2.7.14


seq_path='/scratch/rx32940/lepto_wgs_seq' # path to genomes
o_path='/scratch/rx32940/bwa_results' # path to output files

# index the reference
bwa index $seq_path/Lai_56601.fasta

# align to reference indiviually
for file in $seq_path/*.fasta; do
    echo "in loop"
    echo $file
    isolate=$(echo $file | awk -F'[/.]' '{print $5}')
    echo $isolate
    if [ "$isolate" != "Lai_56601" ]; then
        echo "in if"
        # align the contigs to reference 
        bwa mem -t2 $seq_path/Lai_56601.fasta $seq_path/$isolate.fasta > $o_path/$isolate.sam
        # -o output, convert to binary, bam, format
        samtools sort -o $o_path/$isolate.bam $o_path/$isolate.sam
        # mark duplicates, this step isn't working. necessary for contigs?
        time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates I=$o_path/$isolate.bam O=$o_path/${isolate}_marked_dup.bam M=$o_path/${isolate}_md_metrics.txt
        #index the alignment file (bam)
        samtools index $o_path/${isolate}_marked_dup.bam
    fi
done

time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar CreateSequenceDictionary R=$seq_path/Lai_56601.fasta O=$seq_path/Lai_56601.dict

# GATK cohort variant calling workflow: https://software.broadinstitute.org/gatk/documentation/article.php?id=3893
for file in $o_path/*_marked_dup.bam; do

    isolate=$(echo $file | awk -F'[/._]' '{print $6}')
    echo $isolate

    time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=$file O=$o_path/${isolate}_sorted.bam RGID=$isolate RGLB=WGS RGPL=illumina RGPU=na RGSM=$isolate
    samtools index $o_path/${isolate}_sorted.bam
    gatk HaplotypeCaller -R $seq_path/Lai_56601.fasta -I $o_path/${isolate}_sorted.bam -O $o_path/$isolate.vcf -ERC GVCF

done