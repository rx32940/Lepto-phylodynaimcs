#!/bin/bash
#PBS -q batch                                                            
#PBS -N fastqc_minion                                            
#PBS -l nodes=1:ppn=2 -l mem=20gb                                        
#PBS -l walltime=20:00:00                                                
#PBS -M rx32940@uga.edu                                                  
#PBS -m abe                                                              
#PBS -o /scratch/rx32940                        
#PBS -e /scratch/rx32940                        
#PBS -j oe 

cd $PBS_O_WORKDIR

module load FastQC/0.11.8-Java-1.8.0_144

path="/scratch/rx32940/minion_blood_simulation/test_runs/data/demultiplex2"

cat $path/barcode01/*fastq > $path/barcode01.fastq
fastqc $path/barcode01.fastq -o $path