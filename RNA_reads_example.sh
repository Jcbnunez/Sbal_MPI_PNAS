#!/bin/bash

#SBATCH -J HISAT2_LowB_RNAexp
#SBATCH -n 12
#SBATCH --mem=120GB
#SBATCH -t 30:00:00

#########################################################
#########################################################

##############
#SYSTEM CHECK#
##############
echo "-----------------------------"
echo ">>>>The current folder is<<<<"
pwd
echo ">>>>Operating System Information<<<<"
uname -a
echo ">>>>I/O Information<<<<"
ulimit -a
echo ">>>>Computing Node is<<<<"
hostname
echo ">>>>Domain Name is<<<<"
domainname
echo ">>>>System DNS is<<<<"
dnsdomainname
echo "-----------------------------"

#########################################################

####################
#LOAD OSCAR MODULES#
####################

module load hisat2/2.1.0 #Align reads
module load samtools/1.9 # Manipulate SAM/BAM files
module load bcftools/1.9
module load picard-tools/2.9.2  #Sort and remove duplicates
module load bedtools/2.26.0 
module load bbmap/38.23
module load fgbio/0.6.1 
module load tabix/0.2.6 
module load vcftools/0.1.16
module load seqtk/1.3
module load gatk/4.0.9.0 

#Additional features & Programs
qualimap=~/data/Jcbn/Software/qualimap_v2.1/qualimap
HAPCUT=/users/jnunez/software/HapCUT2/build

#########################################################
#Script variables
CPU=12
JAVAMEM=110g
QUAL=20
PRIOR=2.8e-09
#########################################################

#########################
#DECLARE INPUT VARIABLES#
#########################

#Input Files
F=~/data/S_balanoides_genomics_resources/Reads/LowB.F.trim.rna.fastq.gz
R=~/data/S_balanoides_genomics_resources/Reads/LowB.R.trim.rna.fastq.gz
#head -n 1 $F $R
ProjectName=LowB_RNA_exp

#References
reference=~/data/S_balanoides_genomics_resources/Analyses/MPI_paper/4.DNA_RNA_AA_References/allozymes_as_genes_mtDNA.fasta
#########################################################

########
#SCRIPT#
########
echo "Starting Process ..."
date

#-##Index reference
#-#hisat2-build $reference MPI_ref
#-#
#-## call SNPs with sam tools
#-## the --dta flag is mandatory for differential expression analysis
#-#hisat2 -x MPI_ref -1 $F -2 $R -p $CPU --dta -S $ProjectName.sam
#-#
#-## Raw mapping statistics
#-#samtools flagstat --threads $CPU $ProjectName.sam > flagstats_raw_$ProjectName.txt
#-#
#-## reformat to Bam
#-#samtools view -q $QUAL -F 0x0004 --threads $CPU -b $ProjectName.sam > flt.$ProjectName.bam
#-#
#-## Remove non-mapping reads, fix paired SAM flags added by BWA
#-#samtools flagstat --threads $CPU flt.$ProjectName.bam > flagstats_filter_$ProjectName.txt
#-#
#-## Sort with picard(DNA/RNA)
#-#java -Xmx$JAVAMEM -jar $PICARD SortSam I=flt.$ProjectName.bam O=flt.srt.$ProjectName.bam SO=coordinate VALIDATION_STRINGENCY=SILENT
#-#
#-##flagstats
#-#samtools flagstat --threads $CPU flt.srt.$ProjectName.bam > flagstats_final_$ProjectName.txt
#-#
#-##index
#-#samtools index flt.srt.$ProjectName.bam
#-#
#-##Quality Check
#-#$qualimap bamqc -bam flt.srt.$ProjectName.bam  -outdir ./Qualimap_$ProjectName --java-mem-size=$JAVAMEM
#-#
#-#rm $ProjectName.sam 
#-#rm $ProjectName.bam
#-#rm flt.$ProjectName.bam

###########################
#Begin SNP calling to VCF #
###########################

bcftools mpileup --min-MQ $QUAL -Ou -f $reference -A flt.srt.$ProjectName.bam | bcftools call -m --prior $PRIOR --skip-variants indels  > $ProjectName.allsites.vcf
bgzip $ProjectName.allsites.vcf
tabix $ProjectName.allsites.vcf.gz

########################################################
# Making a reference | due to distance to ME reference #
########################################################

#Otherwise all fixed differences may be lost after phasing and huge biases will ensue!
# It is important to exclude indels in the previous steps, otherwise other individual files will not properly line up

bcftools consensus -H 1 -I -M N -f $reference  $ProjectName.allsites.vcf.gz > $ProjectName.consensus_ref.fasta
samtools faidx $ProjectName.consensus_ref.fasta
 
###############################
# Phasing with PE information #
###############################

gatk SelectVariants -V $ProjectName.allsites.vcf.gz  --exclude-non-variants true -O $ProjectName.variants.vcf.gz

gunzip $ProjectName.variants.vcf.gz
rm $ProjectName.variants.vcf.gz.tbi

$HAPCUT/extractHAIRS --bam flt.srt.$ProjectName.bam --VCF $ProjectName.variants.vcf --out $ProjectName.fragment_file

$HAPCUT/HAPCUT2 --verbose 1 --fragments $ProjectName.fragment_file --vcf $ProjectName.variants.vcf --output $ProjectName.hapcut2.txt

java -jar $FGBIO HapCutToVcf -v $ProjectName.variants.vcf -i $ProjectName.hapcut2.txt -o $ProjectName.phased.vcf

bgzip $ProjectName.variants.vcf
tabix $ProjectName.variants.vcf.gz

rm $ProjectName.fragment_file
#rm $ProjectName.hapcut2.txt
#rm $ProjectName.vcf

####################################
# Compress and tabix the final VCF #
####################################

bgzip $ProjectName.phased.vcf
tabix $ProjectName.phased.vcf.gz

######################
# Extract haplotypes #
######################

# Some extra filetering -- Not Done: only get the phased variants. 
vcftools --gzvcf $ProjectName.phased.vcf.gz --min-alleles 2 --max-alleles 2 --minQ $QUAL --minDP 6  --phased  --recode --recode-INFO-all --out $ProjectName.phasedOnly

bgzip $ProjectName.phasedOnly.recode.vcf
tabix $ProjectName.phasedOnly.recode.vcf.gz

#$vcf2fasta -f $ProjectName.consensus_ref.fasta -P 2 $ProjectName.phasedOnly.recode.vcf

bcftools consensus -H 1 -f $ProjectName.consensus_ref.fasta $ProjectName.phasedOnly.recode.vcf.gz > $ProjectName.haplotype0.fasta
bcftools consensus -H 2 -f $ProjectName.consensus_ref.fasta $ProjectName.phasedOnly.recode.vcf.gz > $ProjectName.haplotype1.fasta

sed -i "s/^>/>$ProjectName.0./g" $ProjectName.haplotype0.fasta 
sed -i "s/^>/>$ProjectName.1./g" $ProjectName.haplotype1.fasta 

echo "Finishing Process"
date
echo "cheers JCBN"
