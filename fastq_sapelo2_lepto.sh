    

cd $PBS_O_WORKDIR

module load SRA-Toolkit/2.9.1-centos_linux64

cat /scratch/rx32940/SRR_Acc_List.txt | xargs -I{} fastq-dump --gzip -O /scratch/rx32940/LeptoFastqSRA {}


