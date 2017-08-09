#!/bin/sh

#load modules
module load hisat2/2.0.5
module load samtools/1.3.1
module load gnu_parallel/201612222


### Change directory to the scratch directory
cd /scratch/CHANGEME/DaphniaPunk

#Only way I could get this to run is when the index files are copied to DaphniaPunk scratch directory by hand. 
##NOTE: it actually works just fine as long as you include the entire path to where the index is location after -x flag: ~/class_shared/DaphniaPunk/Daphnia_pulex_INDEX3 - cheers,Breanna
#parallelize mapping
ls | grep "paired_threads.fastq" |cut -d "_" -f 1,2 | sort | uniq | time parallel -j+0 --eta hisat2 -p 10 --dta -x Daphnia_pulex_INDEX3 -1 {1}_All_R1_paired_threads.fastq {1}_All_R2_paired_threads.fastq -S {1}.sam ::: ${i}


while read i
do

samtools sort -@ 10 -o "$i".bam "$i".sam

done < list
