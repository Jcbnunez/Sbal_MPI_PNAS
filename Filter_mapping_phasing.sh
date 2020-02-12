#Part 2: DNA mapping, filtering, and phasing (Example code)
#For datasets 1, 2, 3, and S. cariosus, we mapped DNA using BWA-MEM (Li 2013) and RNA using HiSat2. Post-processing was done using samtools v 1.9, bcftools v 1.9, and picard-tools. VCF files were made with bcftools using the multiallelic caller (-m; autosomal mutation rate prior of 2.8x10-9). Short read phasing was done using HAPCUT2 (Edge, Bafna et al. 2017). More details in Text S3. Sanger amplicons were phased using the PHASE “MR0” model in DnaSP v 5.1 (Librado and Rozas 2009). All analyses were conducted at Brown University’s supercomputer OSCAR. Using the following software packages

module load bwa/0.7.15 
module load samtools/1.9
module load bcftools/1.9
module load picard-tools/2.9.2  
module load bbmap/38.23
module load fgbio/0.6.1 
module load tabix/0.2.6 
module load gatk/4.0.9.0

#All reads were filtered to remove low quality sequences using bbmap (Bushnell 2016), we enforced a minimum length of 80, trimmed 5 bp on each side, and a minimum mean quality (phred) of 20:

bbduk.sh in=F.fq in2=R.fq out=F.trim.fq out2=R.trim.fq outs=S.trim.fq qtrim=rl trimq=20 maq=20 minlen=80 ftl=5 ftr=145

#We mapped DNA using BWA-MEM (Li 2013). 

bwa mem -M -t $CPU $reference $F $R > $ProjectName.sam

# Where $CPU, $F, and $R are user defined variables for CPU number, forward and reverse reads respectively. $ProjectName is a global variable carried over throughout the pipeline. $reference is a fasta file with the allozymes indexed with bwa.
#We mapped RNA reads using HiSat2. 

hisat2 -x indexed_ref -1 $F -2 $R -p $CPU -S $ProjectName.sam

#Post-processing was done using samtools v 1.9 (Li, Handsaker et al. 2009), bamtools v 1.9 enforcing a minimum mapping quality of 35 ($QUAL) and the following Sam flags: -f 0x0002 (keep reads mapped in proper pairs), -F 0x0004 (remove unmapped reads), -F 0x0008 (remove reads with unmapped mates). 

samtools view -b -q $QUAL -f 0x0002 -F 0x0004 -F 0x0008 --threads $CPU $ProjectName.sam > flt.$ProjectName.bam

#PCR duplicates were removed and bam files sorted by coordinate using picard-tools (https://broadinstitute.github.io/picard/). 
# Sort with picard
java -Xmx$JAVAMEM -jar $PICARD SortSam I=flt.$ProjectName.bam O=flt.srt.$ProjectName.bam SO=coordinate VALIDATION_STRINGENCY=SILENT

# Remove duplicates with picard (For DNA samples)
java -Xmx$JAVAMEM -jar $PICARD MarkDuplicates I=flt.srt.$ProjectName.bam O=flt.rmdp.srt.$ProjectName.bam M=$ProjectName.dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

# where JAVAMEM is pre-allocated memory for java.
#Indexing was done with samtools:

samtools index flt.rmdp.srt.$ProjectName.bam

#The reads deposited to SRA, which corresponds to allozymes and mtDNA, where collected at this stage using:

samtools view -F 12 -b flt.rmdp.srt.$ProjectName.bam | samtools sort -n -O BAM | java -jar $PICARD SamToFastq INTERLEAVE=false I=/dev/stdin F=.$ProjectName.F.mapped.fq F2=.$ProjectName.R.mapped.fq
 
#Reads were uploaded to NCBI using aspera/3.8.1
#Variant calling was performed with bcftools using the multiallelic caller algorithm (-m) and a mutation rate prior (μ_P) of 2.8x10-9 ($PRIOR). This was done exclusively for individual samples. Pooled samples were analyzed using PoPoolation software (see Population genetics, and allozyme scoring, additional info).
# NOTE: The command below will generate a VCF with all sites, both variant and invariant. For variant only use GATK (see below) or the –v flag in bcftools call. Here, $reference is the fasta file containing all allozyme loci.

bcftools mpileup --min-MQ $QUAL -Ou -f $reference -A flt.rmdp.srt.$ProjectName.bam | bcftools call -m --prior $PRIOR --skip-variants indels  > $ProjectName.allsites.vcf

#Filter vcf by quality parameters such as minimum coverage ($MINCOV) and $QUAL. One could always add more filters for greater stringency.

bcftools filter –e '%QUAL< $QUAL || %MAX(DP)<=$MINCOV' $ProjectName.allsites.vcf > $ProjectName.filter.allsites.vcf

bgzip $ProjectName.filter.allsites.vcf
tabix $ProjectName.filter.allsites.vcf.gz

# NOTE: Here, we extracted variant sites using GATK.

gatk SelectVariants -V $ProjectName.filter.allsites.vcf.gz --exclude-non-variants true -O $ProjectName.variants.vcf.gz

#High-throughput sequencing libraries were phased using HAPCUT2 (Edge, Bafna et al. 2017) and parsed with fgbio (https://github.com/fulcrumgenomics/fgbio). 

$HAPCUT/extractHAIRS --bam flt.rmdp.srt.$ProjectName.bam --VCF $ProjectName.variants.vcf --out $ProjectName.fragment_file

$HAPCUT/HAPCUT2 --verbose 1 --fragments $ProjectName.fragment_file --vcf $ProjectName.variants.vcf --output $ProjectName.hapcut2.txt

java -jar $FGBIO HapCutToVcf -v $ProjectName.variants.vcf -i $ProjectName.hapcut2.txt -o $ProjectName.phased.vcf

#The final phased fasta files (for analyses) were generated using the bcftools consensus function.
