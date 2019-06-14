#!/bin/bash
cd "$(dirname "$0")" # change to the working directory where bash file located
# run mummer(nucmer) http://mummer.sourceforge.net/examples/#showsnpsnucmer
# contig based approach

# nucmer -maxmatch -c 100 -p nucmer /Users/rx32940/Desktop/Xu_WGS_lepto/Lai_56601.fasta JQPH01.fasta
# -C option in show-snps assures that only SNPs found in uniquely aligned sequence will be reported, thus excluding SNPs contained in repeats

# alternative approach, filter with delta filter

for file in ./lepto_wgs_seq/*.fasta; do
    
    echo $file 
    isolate=$(echo $file | awk -F'[/.]' '{print $4}')
    echo $isolate
    mkdir ./nucmer/results/$isolate
    nucmer -maxmatch -p ./nucmer/results/$isolate/$isolate ./Lai_56601.fasta ./lepto_wgs_seq/$isolate.fasta
    delta-filter -1 ./nucmer/results/$isolate/$isolate.delta > ./nucmer/results/$isolate/${isolate}_filtered.delta # -r and -q option only require the match set to be consistent with respect to either the reference or query respectively
# -1 is useful for applications that require a 1-to-1 mapping, such as SNP finding (intersection between -r and -q option)
    echo "$isolate alignment done"
    show-coords -r -c -l ./nucmer/results/$isolate/${isolate}_filtered.delta > ./nucmer/results/$isolate/$isolate.coords # summary of the delta file
    mummerplot --png --fat -p ./nucmer/results/$isolate/$isolate ./nucmer/results/$isolate/${isolate}_filtered.delta -R ./Lai_56601.fasta -Q ./lepto_wgs_seq/$isolate.fasta
    show-snps -CTlq -x 1 ./nucmer/results/$isolate/${isolate}_filtered.delta > ./nucmer/results/$isolate/${isolate}_labeled.snps
   
done
# -C, no snps from ambiguous mapping
# -l sequence length information from output
# -q sort output by query ID and SNP position
# - H, no header, -T tab-delimited, both for formatiing purpose

# convert Mummer snps to vcf
#python MUMmerSNPs2VCF.py nucmer_labeled.snps nucmer_labeled.vcf


