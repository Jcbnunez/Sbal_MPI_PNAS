#!/bin/bash

#SBATCH -J D_exon
#SBATCH -n 1
#SBATCH --mem=90GB
#SBATCH -t 12:00:00

module load samtools
#program
Var_at_Pos=~/data/S_balanoides_genomics_resources/Misc_resources/popoolation1/Variance-at-position.pl

#vars
NAME=ME_2011_sbal2
BAM=~/data/S_balanoides_genomics_resources/Analyses/Nunez_et_al_2018_Barnacle_Chapter/3_BWA_mem_BAM_files/ME/flt.rmdp.srt.ME_2011_Sbal2_masked.bam
reference=~/data/S_balanoides_genomics_resources/Analyses/Nunez_et_al_2018_Barnacle_Chapter/1_Genome_Masking/Sbal2_masked.fasta


GTF=~/data/S_balanoides_genomics_resources/Analyses/Nunez_et_al_20XX_MPI_paper/13.Revision2020/3.PI_D_goodSbal2scaffs/GTF_sbal2_exons.gtf


#=#samtools mpileup -q 35 -Q 20 -f $reference $BAM > $NAME.pile

MEASURE=D
OUT=$NAME.exons.D.txt

perl $Var_at_Pos --pileup $NAME.pile --gtf $GTF --output $OUT --measure $MEASURE --pool-size 40 --fastq-type sanger --min-count 2  --min-covered-fraction 0.51 --min-coverage 10


