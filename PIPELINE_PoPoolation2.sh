identify_indel_regions=./popoolation2/indel_filtering/identify-indel-regions.pl
filter_sync_by_gtf=./popoolation2/indel_filtering/filter-sync-by-gtf.pl
mpileup2sync=./popoolation2/mpileup2sync.jar
snp_frequency_diff=./popoolation2/snp-frequency-diff.pl

Project="<define a name>"

samtools mpileup -f $reference -B ME_HI_UH ME_FI_UH RI_D_UH RI_C_UH ME_HI_LC ME_FI_LC RI_D_LC RI_C_LC ICE_MT NOR_MT UKW_MT WCAN_MT > $Project.mpileup 
 
#“$reference” is the path to the reference genes. And “$ME_HI_XX” and/or “$ICE_pool” … are the corresponding microhabitat/regional pools. HI=Hodgsons Island. FI=Farmer’s Island (see Table S1 for their corresponding SRRs and names in NCBI). MT are mid tide samples for biogeographic reference (Iceland, Norway, United Kingdom Wales, Western Canada).

perl $identify_indel_regions --indel-window 5 --min-count 2 --input $Project.mpileup --output $Project.gtf

perl $filter_sync_by_gtf --input $Project.mpileup --gtf $Project.gtf --output $Project.noindels.mpileup

java -ea -Xmx$JAVAMEM -jar $mpileup2sync --input $Project.noindels.mpileup --output $Project.noindels.sync --fastq-type sanger --min-qual 35

perl $snp_frequency_diff --input $Project.noindels.sync --output $Project.noindels.NoSigt.AlleleCount.txt --min-count $MCount --min-coverage $MinX --max-coverage $MaxX
