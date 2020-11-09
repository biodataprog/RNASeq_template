#!/usr/bin/bash
#SBATCH -p intel --mem 32gb -N 1 -n 24 --out logs/STAR_all.log -J STARall

module load STAR
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

INDIR=fastq
OUTDIR=results/STAR_all
IDX=db/STAR
SAMPLEFILE=samples.csv
SPECIES=CimmitisRS
GENOME=db/${SPECIES}_Genome.fasta
GFF=db/${SPECIES}_Genes.gff3
GTF=db/${SPECIES}_Genes.gtf

source config.txt

if [ ! -f $GTF ]; then
	grep -P "\texon\t" $GFF | perl -p -e 's/ID=[^;]+;Parent=([^;]+);/gene_id $1/' > $GTF
fi
if [ ! -d $IDX ]; then
	STAR --runThreadN $CPU --runMode genomeGenerate --genomeDir $IDX --genomeFastaFiles $GENOME \
		--sjdbGTFfile $GTF --genomeSAindexNbases 11
fi
mkdir -p $OUTDIR
OUTNAME=$SPECIES
FWD=$(cut -f4 -d, $SAMPLEFILE | tail -n +2 | perl -p -e "s/(\S+)\n/$INDIR\/\$1\_$FWDEXT.$FASTQEXT,/g" | perl -p -e 's/,$//')
REV=$(cut -f4 -d, $SAMPLEFILE | tail -n +2 | perl -p -e "s/(\S+)\n/$INDIR\/\$1\_$REVEXT.$FASTQEXT,/g" | perl -p -e 's/,$//')
echo "$FWD $REV"
SAMHDR=$(cut -f1 -d, $SAMPLEFILE | tail -n +2 | perl -p -e 's/(\S+)\n/ID:$1 , /' | perl -p -e 's/ , $//')

STAR    --outSAMstrandField intronMotif --runThreadN $CPU \
	--outMultimapperOrder Random --twopassMode Basic \
	--genomeDir $IDX --outFileNamePrefix "$OUTDIR/$OUTNAME." \
	--readFilesCommand zcat  --outSAMattrRGline "$SAMHDR" \
	--readFilesIn "$FWD" "$REV"
