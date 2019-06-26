#!/bin/bash
#PBS -q batch                                                            
#PBS -N lepto_fastq                                            
#PBS -l nodes=1:ppn=2 -l mem=20gb                                        
#PBS -l walltime=10:00:00                                                
#PBS -M rx32940@uga.edu                                                  
#PBS -m abe                                                              
#PBS -o /scratch/rx32940/fastq_seq                        
#PBS -e /scratch/rx32940/fastq_seq                        
#PBS -j oe     

cd $PBS_O_WORKDIR

module load SRA-Toolkit/2.9.1-centos_linux64

cat SRR_Acc_List.txt | xargs -I{} fastq-dump -O /scratch/rx32940/fastq_seq {}


