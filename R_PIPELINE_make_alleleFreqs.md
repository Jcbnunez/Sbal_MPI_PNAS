#Load libraries
```{r}
library(tidyverse)
library(magrittr)
library(reshape2)
library(stringr)
library(vcfR)
library(GenomicFeatures)
library(VariantAnnotation)
```

# Load pool-seq data
```{r}
allele_counts_raw= read.table("~/data/S_balanoides_genomics_resources/Analyses/Nunez_et_al_20XX_MPI_paper/13.Revision2020/4.regulatory_evo/FileS1.txt", sep= "\t", header = F,  na.strings = "na")

allele_counts_raw[,c(1:3,4:5,8,9:c(9+12))] -> allele_counts
#"ME_HI_UH", "ME_FI_UH", "RI_D_UH", "RI_C_UH", "ME_HI_LC", "ME_FI_LC", "RI_D_LC", "RI_C_LC", "ICE_MT", "NOR_MT", "UKW_MT", "WCAN_MT"

names(allele_counts)= c("chr","pos","rc","allele_count","allele_states","major_alleles.maa","minor_alleles.mia","maa_ME_HI_UH", "maa_ME_FI_UH", "maa_RI_D_UH", "maa_RI_C_UH", "maa_ME_HI_LC", "maa_ME_FI_LC", "maa_RI_D_LC", "maa_RI_C_LC", "maa_ICE_MT", "maa_NOR_MT", "maa_UKW_MT", "maa_WCAN_MT")

allele_counts  %>%  separate(major_alleles.maa, into = c("ME_HI_UH_ma", "ME_FI_UH_ma", "RI_D_UH_ma", "RI_C_UH_ma", "ME_HI_LC_ma", "ME_FI_LC_ma", "RI_D_LC_ma", "RI_C_LC_ma", "ICE_MT_ma", "NOR_MT_ma", "UKW_MT_ma", "WCAN_MT_ma"), sep = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), remove = F) %>%  separate(minor_alleles.mia, into = c("ME_HI_UH_mi", "ME_FI_UH_mi", "RI_D_UH_mi", "RI_C_UH_mi", "ME_HI_LC_mi", "ME_FI_LC_mi", "RI_D_LC_mi", "RI_C_LC_mi", "ICE_MT_mi", "NOR_MT_mi", "UKW_MT_mi", "WCAN_MT_mi"), sep = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), remove = F) %>% separate(maa_ME_HI_UH, into = c("ME_HI_UH_count","ME_HI_UH_cov"), sep="/") %>% separate(maa_ME_FI_UH, into = c("ME_FI_UH_count","ME_FI_UH_cov"), sep="/") %>% separate(maa_RI_D_UH, into = c("RI_D_UH_count","RI_D_UH_cov"), sep="/") %>% separate(maa_RI_C_UH, into = c("RI_C_UH_count","RI_C_UH_cov"), sep="/") %>% separate(maa_ME_HI_LC, into = c("ME_HI_LC_count","ME_HI_LC_cov"), sep="/") %>% separate(maa_ME_FI_LC, into = c("ME_FI_LC_count","ME_FI_LC_cov"), sep="/") %>% separate(maa_RI_D_LC, into = c("RI_D_LC_count","RI_D_LC_cov"), sep="/") %>% separate(maa_RI_C_LC, into = c("RI_C_LC_count","RI_C_LC_cov"), sep="/") %>% separate(maa_ICE_MT, into = c("ICE_MT_count","ICE_MT_cov"), sep="/") %>% separate(maa_NOR_MT, into = c("NOR_MT_count","NOR_MT_cov"), sep="/") %>% separate(maa_UKW_MT, into = c("UKW_MT_count","UKW_MT_cov"), sep="/") %>% separate(maa_WCAN_MT, into = c("WCAN_MT_count","WCAN_MT_cov"), sep="/") -> allele_counts_sep
```



## add ancestral state
```{r}
CARvcf=read.vcfR("~/data/S_balanoides_genomics_resources/Analyses/Nunez_et_al_20XX_MPI_paper/5.BAM_files/CAR-sr/Alloz_CAR_ind.allsites.vcf.gz")

CARvcf@fix %>% .[,c(1,2,4:6)] -> CAR_validation_raw
CARvcf@gt -> CAR_validation_gt

cbind(CAR_validation_raw,CAR_validation_gt) -> CAR_Ancestral_info
CAR_Ancestral_info=as.data.frame(CAR_Ancestral_info)
names(CAR_Ancestral_info)[c(1:4,7) ] = c("chr","pos","ref_mapping","ancestral_alt","car_Genotype")

CAR_Ancestral_info$chr = as.character(CAR_Ancestral_info$chr)
CAR_Ancestral_info$car_Genotype = as.character(CAR_Ancestral_info$car_Genotype)
CAR_Ancestral_info$pos = as.numeric(as.character(CAR_Ancestral_info$pos))

CAR_Ancestral_info %<>% separate(ancestral_alt, into = c("ancestral_alt1","ancestral_alt2"), sep = ",")
CAR_Ancestral_info$ancestral_alt1 = as.character(CAR_Ancestral_info$ancestral_alt1)
CAR_Ancestral_info$ancestral_alt2 = as.character(CAR_Ancestral_info$ancestral_alt2)
CAR_Ancestral_info$ref_mapping = as.character(CAR_Ancestral_info$ref_mapping)


CAR_Ancestral_info %<>% mutate(ancetral_rc = ifelse(is.na(CAR_Ancestral_info$ancestral_alt1), print(CAR_Ancestral_info$ref_mapping), print(CAR_Ancestral_info$ancestral_alt1)))

left_join(allele_counts_sep,CAR_Ancestral_info ) -> allele_counts_sep_polarized

allele_counts_sep_polarized %<>% separate(car_Genotype, into = c("car_Genotype","car_gt_info"), sep = ":")

#allele_counts_sep_polarized %>% .[which(.$car_Genotype %in% c("0/1","1/0")),]
```

## add missing data
```{r}
allele_counts_sep_polarized %>% mutate(Missing_sites = str_count(major_alleles.maa,"N")) %>% mutate(Nonpolymorphic_sites = str_count(minor_alleles.mia,"N")) -> allele_counts_sep_polarized
```

## call p
```{r}
data_in = allele_counts_sep_polarized

data_in$ancetral_rc = as.character(data_in$ancetral_rc)

counts = c("ME_HI_UH_count","ME_FI_UH_count","RI_D_UH_count","RI_C_UH_count","ME_HI_LC_count","ME_FI_LC_count","RI_D_LC_count","RI_C_LC_count","ICE_MT_count","NOR_MT_count","UKW_MT_count","WCAN_MT_count")
covs = c("ME_HI_UH_cov","ME_FI_UH_cov","RI_D_UH_cov","RI_C_UH_cov","ME_HI_LC_cov","ME_FI_LC_cov","RI_D_LC_cov","RI_C_LC_cov","ICE_MT_cov","NOR_MT_cov","UKW_MT_cov","WCAN_MT_cov")

data_in[counts] = lapply(data_in[counts], as.numeric)
data_in[covs] = lapply(data_in[covs], as.numeric)
lapply(data_in[counts], class)
lapply(data_in[covs], class)

data_in %>% mutate(MEHIUH_p = ifelse(data_in$ME_HI_UH_ma == data_in$ancetral_rc, print(data_in$ME_HI_UH_count/data_in$ME_HI_UH_cov), print(1- (data_in$ME_HI_UH_count/data_in$ME_HI_UH_cov))))  %>% mutate(MEFIUH_p = ifelse(data_in$ME_FI_UH_ma == data_in$ancetral_rc, print(data_in$ME_FI_UH_count/data_in$ME_FI_UH_cov), print(1- (data_in$ME_FI_UH_count/data_in$ME_FI_UH_cov)))) %>% mutate(MEHILC_p = ifelse(data_in$ME_HI_LC_ma == data_in$ancetral_rc, print(data_in$ME_HI_LC_count/data_in$ME_HI_LC_cov), print(1- (data_in$ME_HI_LC_count/data_in$ME_HI_LC_cov)))) %>% mutate(MEFILC_p = ifelse(data_in$ME_FI_LC_ma == data_in$ancetral_rc, print(data_in$ME_FI_LC_count/data_in$ME_FI_LC_cov), print(1- (data_in$ME_FI_LC_count/data_in$ME_FI_LC_cov)))) %>% mutate(RICLC_p = ifelse(data_in$RI_C_LC_ma == data_in$ancetral_rc, print(data_in$RI_C_LC_count/data_in$RI_C_LC_cov), print(1- (data_in$RI_C_LC_count/data_in$RI_C_LC_cov)))) %>% mutate(RICUH_p = ifelse(data_in$RI_C_UH_ma == data_in$ancetral_rc, print(data_in$RI_C_UH_count/data_in$RI_C_UH_cov), print(1- (data_in$RI_C_UH_count/data_in$RI_C_UH_cov)))) %>% mutate(RIDUH_p = ifelse(data_in$RI_D_UH_ma == data_in$ancetral_rc, print(data_in$RI_D_UH_count/data_in$RI_D_UH_cov), print(1- (data_in$RI_D_UH_count/data_in$RI_D_UH_cov)))) %>% mutate(RIDLC_p = ifelse(data_in$RI_D_LC_ma == data_in$ancetral_rc, print(data_in$RI_D_LC_count/data_in$RI_D_LC_cov), print(1- (data_in$RI_D_LC_count/data_in$RI_D_LC_cov)))) %>% mutate(ICE_p = ifelse(data_in$ICE_MT_ma == data_in$ancetral_rc, print(data_in$ICE_MT_count/data_in$ICE_MT_cov), print(1- (data_in$ICE_MT_count/data_in$ICE_MT_cov)))) %>% mutate(UKW_p = ifelse(data_in$UKW_MT_ma == data_in$ancetral_rc, print(data_in$UKW_MT_count/data_in$UKW_MT_cov), print(1- (data_in$UKW_MT_count/data_in$UKW_MT_cov)))) %>% mutate(NOR_p = ifelse(data_in$NOR_MT_ma == data_in$ancetral_rc, print(data_in$NOR_MT_count/data_in$NOR_MT_cov), print(1- (data_in$NOR_MT_count/data_in$NOR_MT_cov))))%>% mutate(CAN_p = ifelse(data_in$WCAN_MT_ma == data_in$ancetral_rc, print(data_in$WCAN_MT_count/data_in$WCAN_MT_cov), print(1- (data_in$WCAN_MT_count/data_in$WCAN_MT_cov)))) -> data_out

allele_counts_sep_p = data_out

```

## Add GFF mapping
```{r}
makeTxDbFromGFF("~/data/S_balanoides_genomics_resources/Analyses/Nunez_et_al_20XX_MPI_paper/3.Gene_predictions_Flowerdew/RNA_seq_pred.RNA_Official.gff3") -> Sbal2_gff_mpi
#Create a FaFi file
allo_FaFI <- FaFile("~/data/S_balanoides_genomics_resources/Analyses/Nunez_et_al_20XX_MPI_paper/4.DNA_RNA_AA_References/allozymes_as_genes_mtDNA.fasta")

#### === > Beging Mapping -- General Mapping
makeGRangesFromDataFrame(allele_counts_sep_p, start.field="pos", end.field="pos") -> allele_counts_raw_asVCF3

print("Predicting All vars")
All_var3 <- locateVariants(allele_counts_raw_asVCF3, Sbal2_gff_mpi, AllVariants())
print("Done Predicting All vars")

All_var_df3 = data.frame(All_var3)
names(All_var_df3)[1:2] = c("chr","pos")

left_join(allele_counts_sep_p, All_var_df3) -> allele_counts_sep_p_annot

###### Functional SNPS

allele_counts_sep_p -> tmp_snps_poly

tmp_snps_poly %<>% separate(allele_states, into = c("state1","state2"), sep = "\\/")
tmp_snps_poly$state1 = as.character(tmp_snps_poly$state1)
tmp_snps_poly$state2 = as.character(tmp_snps_poly$state2)
tmp_snps_poly$rc = as.character(tmp_snps_poly$rc)

tmp_snps_poly[,c("chr","pos","state1","state2","rc")] -> tmp_snps_poly_asVCF
tmp_snps_poly_asVCF$rc = as.character(tmp_snps_poly_asVCF$rc)
makeGRangesFromDataFrame(tmp_snps_poly_asVCF, start.field="pos", end.field="pos") -> tmp_snps_poly_GR

GR1_VCF <- GRanges(tmp_snps_poly_GR,
SNP=DNAStringSet(ifelse(tmp_snps_poly$state1 == tmp_snps_poly_asVCF$rc, paste(tmp_snps_poly$state2), paste(tmp_snps_poly$state1)))
)
GR1_VCF
print("GR1_VCF made")

print("Predicting Coding Consequences")
coding_poly <- predictCoding(GR1_VCF, Sbal2_gff_mpi, allo_FaFI, GR1_VCF$SNP )
print("Done Predicting Coding Consequences")

coding_df = data.frame(coding_poly)
names(coding_df)[1:2] = c("chr","pos")

left_join(allele_counts_sep_p_annot, coding_df[,c("chr","pos","CDSLOC.start", "CDSLOC.end","CDSLOC.width","PROTEINLOC","CONSEQUENCE","REFCODON","VARCODON",     "REFAA","VARAA")]) -> allele_counts_sep_p_annot

write.table(allele_counts_sep_p_annot, file = "./allele_counts_sep_p_annot.txt", sep = "\t",quote = F ,row.names = F )
```