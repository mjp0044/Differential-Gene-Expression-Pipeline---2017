#!/bin/sh
#
# load required modules
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load hisat2/2.0.5
module load cufflinks/2.2.1
module load python/3.3.2

#genome location: /apps/bio/unzipped/jgi/Dappu1
	#files include: (frozen gene catalogue are all manually curated genes as of 2007)
	#Daphnia_pulex.allmasked, Daphnia_pulex.fasta, FrozenGeneCatalog20110204.CDS.fasta
	#FrozenGeneCatalog20110204.gff, FrozenGeneCatalog20110204.gff3,
	#FrozenGeneCatalog20110204.proteins.fasta, FrozenGeneCatalog20110204.transcripts.fasta

#make scratch directory
#mkdir /scratch/DaphniaPunk

#move to scratch directory
cd /scratch/DaphniaPunk

#convert gff3 to gtf (have to convert the gff3 for it to work for some reason)
gffread /apps/bio/unzipped/jgi/Dappu1/FrozenGeneCatalog20110204.gff3 -T -o FrozenGeneCatalog20110204.gtf

#extract splice site information from gene annotationexon files
hisat2_extract_splice_sites.py FrozenGeneCatalog20110204.gtf >FrozenGeneCatalog20110204.ss

#extract exon information from gene annotation files
hisat2_extract_exons.py FrozenGeneCatalog20110204.gtf >FrozenGeneCatalog20110204.exon

#build an indexed reference genome
#for human, 160GB of memory, 8 cores, and 2hrs (-p allows it to use multiple cores)
hisat2-build -p 16 --ss FrozenGeneCatalog20110204.ss --exon FrozenGeneCatalog20110204.exon /apps/bio/unzipped/jgi/Dappu1/Daphnia_pulex.fasta Daphnia_pulex_INDEX2

#ran with 16 cores and 20gb memory and only took about 3min