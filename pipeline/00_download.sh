#!/usr/bin/bash
#SBATCH -p short

# FIX THIS TO SPECIES YOU WANT FROM FUNGIDB 
# OR REWRITE THIS SCRIPT TO DOWNLOAD WHAT YOU WANT
SPECIES=CimmitisRS
source config.txt

VERSION=49
DB=db
mkdir -p $DB
URL=https://fungidb.org/common/downloads/release-${VERSION}
MRNA=${SPECIES}_mRNA.fasta
GENOME=${SPECIES}_Genome.fasta
GFF=${SPECIES}_Genes.gff3

if [ ! -f $DB/$MRNA ]; then
	curl -o $DB/$MRNA $URL/$SPECIES/fasta/data/FungiDB-${VERSION}_${SPECIES}_AnnotatedTranscripts.fasta 
fi

if [ ! -f $DB/$GENOME ]; then
	 curl -o $DB/$GENOME $URL/$SPECIES/fasta/data/FungiDB-${VERSION}_${SPECIES}_Genome.fasta
fi

if [ ! -f $DB/$GFF ]; then
	curl -o $DB/$GFF $URL/$SPECIES/gff/data/FungiDB-${VERSION}_${SPECIES}.gff
fi
