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

module load GATK/4.0.11.0-foss-2016b-Python-2.7.14

seq_path='/scratch/rx32940/lepto_wgs_seq' # path to genomes
o_path='/scratch/rx32940/bwa_results' # path to output files

gatk HaplotypeCaller -R $seq_path/Lai_56601.fasta -I $o_path/JQOL00000000_sorted.bam -O $o_path/JQOL00000000_test.vcf -ERC GVCF