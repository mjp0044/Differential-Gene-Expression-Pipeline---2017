#! /bin/sh
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load stringtie/1.3.3
module load gnu_parallel/201612222

#run with 
#10cores
#2 hours
#20gb 

#change into scratch directory with all the files
cd /scratch/CHANGME/DaphniaPunk

#copy over annotation file daphnia_genes2010_beta3.gtf before running script!!

#Assemble transcripts for each sample
ls | grep ".bam" | cut -d "_" -f 1 | sort | uniq | time parallel -j+0 --eta stringtie -p 10 -G daphnia_genes2010_beta3.gtf -o {1}assembled.gtf -l {1} {1}*.bam ::: ${i}

#Get mergelist of assembled transcripts, letter 'a' added to transcript filename to allow grab of .gtf files ignoring index
ls *assembled.gtf | sort | uniq > mergelist.txt

#Merge transcripts from all samples
stringtie --merge -p 10 -G daphnia_genes2010_beta3.gtf -o stringtie_merged.gtf mergelist.txt

#Make appropriate directory heirarchy for ballgown
mkdir ballgown

#Estimate transcript abundances and create table counts for Ballgown
ls | grep ".bam" | cut -d "_" -f 1 | sort | uniq | time parallel -j+0 --eta stringtie -p 10 -G stringtie_merged.gtf -o ballgown/daph{1}/daph{1}.gtf -B -e {1}*.bam ::: ${i}

#Make tarball to bring over to R
tar -cvzf ballgown.tgz ballgown

