#Part 3: Population genetics, and allozyme scoring, additional info
#For all datasets, we filtered all singleton loci, and enforced a 5% minimum allele frequency cut-off. Pool-seq datasets were analyzed using the PoPoolation1 (Kofler, Orozco-terWengel et al. 2011) and PoPoolation2 (Kofler, Pandey et al. 2011) tool-kits. PoPoolation1 was used to estimate population genetic parameters such as π, θ, and Tajima’s D. For these analyses, we enforced the flags: --min-count 2 --min-coverage 10 --max-coverage 200. We subsample coverage to a set value. For Popoolation1 analysis Maine and Rhode Island pool this value was 100X. For other populations, this value was 20X. Please note that the table provided in File S1, was not normalized to uniform coverage, instead it shows the raw coverage values for each site.  Here is an example of our pipeline for PoPoolation2:

#-## Create a pileup file
samtools mpileup -f $reference -B BAM1 BAM2 … > $ProjectName.mpileup

#-##filter indels part1: create a map file
perl  $popoolation2/indel_filtering/identify-indel-regions.pl --indel-window 5 --min-count 2 --input $ProjectName.mpileup --output $ProjectName.gtf

#-##filter indels part2: filter index via gft map
perl $popoolation2/indel_filtering/filter-sync-by-gtf.pl --input $ProjectName.mpileup --gtf $ProjectName.gtf --output $ProjectName.noindels.mpileup

#-## Make synchronized pileup
java -ea -Xmx$JAVAMEM -jar $popoolation2/mpileup2sync.jar --input $ProjectName.noindels.mpileup --output $ProjectName.noindels.sync --fastq-type sanger --min-qual 35 --threads $CPU

#-## Eliminate Singleton calls in the sync file
#-## Its impossible to differentiate true singeltons/doubletons from sequencing errors. Thus, we will remove them from the Sync File.
cat $ProjectName.noindels.sync | sed -E 's|:1:|:0:|g' | sed -E 's|\t1:|\t0:|g' | sed -E 's|:2:|:0:|g' | sed -E 's|\t2:|\t0:|g' > $ProjectName.noindels.NoSigt.sync

#-##Subsample to uniform coveraage
perl ~/data/Jcbn/Software/popoolation2/subsample-synchronized.pl --input $ProjectName.noindels.sync --output $ProjectName.noindels.$targetCOV.sync --target-coverage $targetCOV  --max-coverage $MAXCOV --method withoutreplace

#the analysis table was generated using …/popoolation2/snp-frequency-diff.pl on the sync file.