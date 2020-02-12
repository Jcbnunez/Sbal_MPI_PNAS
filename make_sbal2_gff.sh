#!/bin/bash

#SBATCH -J Sbal2_RNA_informed_predictions
#SBATCH -n 2
#SBATCH --mem=60GB
#SBATCH -t 100:00:00

# Load modules
module load samtools/1.9
module load bamtools/2.4.1
module load augustus/3.3 
module load bwa/0.7.15 
module load bedtools/2.26.0
module load hisat2/2.1.0 
module load exonerate/2.2.0 
module load picard-tools/2.9.2  #Sort and remove duplicates

#User defined variables

PROJECT=Sbal2_RNAprior
reference=~/data/S_balanoides_genomics_resources/Analyses/Nunez_et_al_2018_Barnacle_Chapter/1_Genome_Masking/Sbal2_masked.fasta
BAM=~/data/S_balanoides_genomics_resources/Analyses/Nunez_et_al_2018_Barnacle_Chapter/6_RNA_confrimation/flt.rmdp.srt.ME_RNA_cat_SBAL2.bam
CPU=2

##### Extract RNA evidence #####

samtools sort -n $BAM > $PROJECT.accept_hits.bam
# filter alignments with filterBam
filterBam --uniq --in $PROJECT.accept_hits.bam --out $PROJECT.acchits.flt.bam
samtools view -H $PROJECT.acchits.flt.bam > $PROJECT.header.txt
samtools sort $PROJECT.acchits.flt.bam > $PROJECT.both.ssf.bam
bam2hints --intronsonly --in=$PROJECT.both.ssf.bam --out=$PROJECT.hints.gff


##### Augustus Raw 1

augustus --species=fly --extrinsicCfgFile=extrinsic.M.RM.E.W.P.cfg --alternatives-from-evidence=true --hintsfile=$PROJECT.hints.gff --allow_hinted_splicesites=atac --introns=on --genemodel=complete $reference > $PROJECT.round1.gff

##### Create an exon-exon junction database

cat $PROJECT.round1.gff | tee $PROJECT.aug.prelim.gff | grep -P "\tintron\t" > $PROJECT.aug1.introns.gff

# in case hints.gff also contains other hints than "intron", you need to filter for "intron", first!
cat $PROJECT.hints.gff $PROJECT.aug1.introns.gff > $PROJECT.introns.lst
intron2exex.pl --introns=$PROJECT.introns.lst --seq=$reference --exex=$PROJECT.exex.fa --map=$PROJECT.map.psl

# Extract RNA reads from bam file
bedtools bamtofastq -i $BAM -fq $PROJECT.FR.reads.fq 

#Map with BWA mem
bwa index $PROJECT.exex.fa
bwa mem -t $CPU $PROJECT.exex.fa $PROJECT.FR.reads.fq > $PROJECT.exonjunlib.sam
samtools view -S -F 4 $PROJECT.exonjunlib.sam > flt.$PROJECT.exonjunlib.sam
samMap.pl flt.$PROJECT.exonjunlib.sam $PROJECT.map.psl > $PROJECT.BWAMEM.global.sam

# discard intron containing alignments from the original bam file
operation_N_filter=~/data/AUGUSTUS_LOCAL/3.0.2/augustus-3.0.2/auxprogs/auxBamFilters/operation_N_filter.txt
bamtools filter -in $PROJECT.accept_hits.bam -out $PROJECT.accept_hits.noN.bam -script $operation_N_filter

# create a bam file with header from the bowtie.global.sam file
cat $PROJECT.header.txt $PROJECT.BWAMEM.global.sam > $PROJECT.BWAMEM.global.h.sam
samtools view -bS -o $PROJECT.BWAMEM.global.h.bam $PROJECT.BWAMEM.global.h.sam

# join bam files
bamtools merge -in $PROJECT.BWAMEM.global.h.bam -in $PROJECT.accept_hits.noN.bam -out $PROJECT.both.Evidence.bam
samtools sort -n $PROJECT.both.Evidence.bam > srt.$PROJECT.both.Evidence.bam
filterBam --uniq --in srt.$PROJECT.both.Evidence.bam --out flt.srt.$PROJECT.both.Evidence.bam
 
# Creating intron hints (step 2)
samtools sort flt.srt.$PROJECT.both.Evidence.bam > flt.srt2.$PROJECT.both.Evidence.bam
bam2hints --intronsonly --in=flt.srt2.$PROJECT.both.Evidence.bam --out=$PROJECT.hints.Iterated.gff

####################################################

#Predict with Augusus
augustus --species=fly --extrinsicCfgFile=extrinsic.M.RM.E.W.P.cfg --alternatives-from-evidence=true --hintsfile=$PROJECT.hints.Iterated.gff --allow_hinted_splicesites=atac --gff3=on $reference > $PROJECT.iteratedHints.gff3
