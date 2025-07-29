#GDC_patient_mutation.r

#2025.7.28
### aims:
#(1) check effect of mutations
### conclusion:
#(1) 
### reference:
#(1) 





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### 1: INITIALIZE #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


### library
library("GenomicRanges")
library("tidyr")
library("dplyr")
library("ggplot2")
library("stringr")
library(vcfR)


### output parameters

## output directories
out_path <- "output" ###!!!
fig_path <- "figs" ###!!!
system("mkdir -p output") ###!!!
system("mkdir -p figs") ###!!!



## aesthetic setting for plots
# SETTING FIGURE SIZE
one.c <- 90   #single column
one.5c <- 140 #1.5 column
two.c <- 190  #full width
# SETIING TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### 2: annotate gdc data by mutation type #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

### readVcF
library(VariantAnnotation)
vcf <- readVcf("gdc.gatk4_wgs.chou_orf143.norm.vcf",yieldSize=100)
vcf


### basic checking

## readVcf_showheader
header(vcf)

## headeraccessors
samples(header(vcf))
geno(header(vcf))

## readVcf_rowRanges
# head(rowRanges(vcf), 3)

## readVcf_fixed
ref(vcf)[1:5]
qual(vcf)[1:5]

### readVcf_ALT
alt(vcf)[1:5]


### does the data include intron?
# test with the intron of c16norep142
a <- ranges(rowRanges(vcf))|>start()>89716905
b <- ranges(rowRanges(vcf))|>start()<89717343
rowRanges(vcf)[a&b]
## Note: 


### predictCoding

## make TxDB
library(GenomicFeatures)
transcripts <- data.frame(
  tx_id=1:4,
  tx_chrom=c("chr16","chr16","chr16","chr16"),
  tx_strand=c("+","+","+","+"),
  tx_start=c(89711856,89711856,89711856,89711856),
  tx_end=c(89718165,89718165,89718165,89718165))
splicings <-  data.frame(
  tx_id=c(1L, 1L, 2L,3L,4L,4L),
  exon_rank=c(1, 2, 1,1,1,2), #  'splicings' must contain unique (tx_id, exon_rank) pairs
  exon_start=c(89716582, 89717344,89716634,89717451,89716888,89717344),
  exon_end=c(89716637, 89717386,89716891,89718032,89716904,89717386),
  # cds_id=c(1L,1L, 2L),
  cds_start=c(89716582, 89717344,89716634,89717451,89716888,89717344),
  cds_end=c(89716637, 89717386,89716891,89718032,89716904,89717386),
  cds_phase=c(NA, NA, NA,NA,NA,NA))
genes <-  data.frame(
  tx_id=c(1L, 1L, 2L,3L,4L,4L),
  gene_id=c("c16riboseqorf141","c16riboseqorf141", "c16riboseqorf142","c16riboseqorf143","c16norep142","c16norep142")
)
# assemble
txdb <- makeTxDb(transcripts, splicings,genes)

# predict coding
library(BSgenome.Hsapiens.UCSC.hg38)
coding <- predictCoding(vcf, txdb, seqSource=Hsapiens)
coding


### summarize the impact on protein

## coding info
# summarize the table
coding[duplicated(names(coding))]
# GRanges object with 2 ranges and 17 metadata columns:
#   seqnames    ranges strand | paramRangeID            REF                ALT      QUAL      FILTER
# <Rle> <IRanges>  <Rle> |     <factor> <DNAStringSet> <DNAStringSetList> <numeric> <character>
#   chr16:89717346_C/T    chr16  89717346      + |           NA              C                  T        NA           .
# chr16:89717361_G/A    chr16  89717361      + |           NA              G                  A        NA           .
# varAllele    CDSLOC    PROTEINLOC   QUERYID        TXID         CDSID      GENEID   CONSEQUENCE
# <DNAStringSet> <IRanges> <IntegerList> <integer> <character> <IntegerList> <character>      <factor>
#   chr16:89717346_C/T              T        20             7        17           4             4 c16norep142 nonsynonymous
# chr16:89717361_G/A              A        35            12        18           4             4 c16norep142 nonsynonymous
# REFCODON       VARCODON         REFAA         VARAA
# <DNAStringSet> <DNAStringSet> <AAStringSet> <AAStringSet>
#   chr16:89717346_C/T            CCT            CTT             P             L
# chr16:89717361_G/A            CGT            CAT             R             H
# -------
#   seqinfo: 2779 sequences from an unspecified genome
coding[names(coding)=="chr16:89717346_C/T"]
# note: two mutation duplicated because they are in two ORFs

## make a df
PH01<- data.frame(ranges(coding))
coding_info <- data.frame(#mut_name = names(coding),
                          GENEID   = coding$GENEID,
                          CONSEQUENCE = coding$CONSEQUENCE)
coding_info<-cbind(coding_info,PH01)





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### 3: have quick look of the mutations #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


### is there are any mutations in TIS ranges of each ORF (which covers the start codon)

## create TIF ranges: -6 to +2
orf_TIS <- GRanges(
  seqnames=c("chr16", "chr16", "chr16", "chr16"),
  ranges=IRanges(c(89716582,89716634,89717451,89716888)-6, width=11, names=c("c16riboseqorf141_TIS","c16riboseqorf142_TIS","c16riboseqorf143_TIS","c16norep142_TIS")),
  strand="+"
)
orf_TIS

findOverlaps(orf_TIS,rowRanges(vcf), type = "any")

# Hits object with 0 hits and 0 metadata columns:
## no overlap. I added a positive control and the script works. There's no hit

## see if there are mutations in a bigger window: +- 10 range
orf_TIS <- GRanges(
  seqnames=c("chr16", "chr16", "chr16", "chr16"),
  ranges=IRanges(c(89716582,89716634,89717451,89716888)-10, width=23, names=c("c16riboseqorf141_TIS","c16riboseqorf142_TIS","c16riboseqorf143_TIS","c16norep142_TIS")),
  strand="+"
)
orf_TIS
ov <- findOverlaps(orf_TIS,rowRanges(vcf))
orf_TIS[queryHits(ov)]
# GRanges object with 1 range and 0 metadata columns:
#   seqnames            ranges strand
# <Rle>         <IRanges>  <Rle>
#   c16riboseqorf142_TIS    chr16 89716624-89716646      +
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths


### type of consequences in terms of impact of ORF
table(coding$CONSEQUENCE)


### check the what is the consequence of frameshift
coding[coding$CONSEQUENCE=="frameshift"]

# frameshift      nonsense nonsynonymous    synonymous 
# 9             2            39            24 

### check the consequence of nonsense mutation

subset(coding,CONSEQUENCE=="nonsense")
# GRanges object with 2 ranges and 17 metadata columns:
#   seqnames    ranges strand | paramRangeID            REF                ALT      QUAL      FILTER
# <Rle> <IRanges>  <Rle> |     <factor> <DNAStringSet> <DNAStringSetList> <numeric> <character>
#   chr16:89716700_C/T    chr16  89716700      + |           NA              C                  T        NA           .
# chr16:89717784_C/T    chr16  89717784      + |           NA              C                  T        NA           .
# varAllele    CDSLOC    PROTEINLOC   QUERYID        TXID         CDSID           GENEID CONSEQUENCE
# <DNAStringSet> <IRanges> <IntegerList> <integer> <character> <IntegerList>      <character>    <factor>
#   chr16:89716700_C/T              T        67            23         9           2             2 c16riboseqorf142    nonsense
# chr16:89717784_C/T              T       334           112        53           3             5 c16riboseqorf143    nonsense
# REFCODON       VARCODON         REFAA         VARAA
# <DNAStringSet> <DNAStringSet> <AAStringSet> <AAStringSet>
#   chr16:89716700_C/T            CAG            TAG             Q             *
#   chr16:89717784_C/T            CGA            TGA             R             *
#   -------
#   seqinfo: 2779 sequences from an unspecified genome





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### 4: frequency of the mutations: gdc #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#



### read the mutation loc and count info
in_exons_indel <- read.table("in_exons_indel.tsv", sep="\t",header=1)
in_exons_indel_adeno <- read.table("in_exons_indel_adeno.tsv", sep="\t",header=1)
in_exons_SNP <- read.table("in_exons_snp.tsv", sep="\t",header=1)
in_exons_SNP_adeno <- read.table("in_exons_snp_adeno.tsv", sep="\t",header=1)


##
in_exons_indel$mutation_type <- "INDEL"
in_exons_SNP$mutation_type <- "SNP"
#
in_exons_indel_adeno$mutation_type <- "INDEL"
in_exons_SNP_adeno$mutation_type <- "SNP"


##
in_exons_indel$out_orf <- ""
mut_loc_count <- rbind(in_exons_SNP,in_exons_indel)
#
in_exons_indel_adeno$out_orf <- ""
mut_loc_count_adeno <- rbind(in_exons_SNP_adeno,in_exons_indel_adeno)



## the table was made with python. change to 1-indexed
mut_loc_count$Start <- mut_loc_count$Start+1
mut_loc_count_adeno$Start <- mut_loc_count_adeno$Start+1
#
mut_loc_count$index<-NULL
mut_loc_count_adeno$index<-NULL
mut_loc_count$Chromosome<-NULL
mut_loc_count_adeno$Chromosome<-NULL

#
names(mut_loc_count) <- tolower(names(mut_loc_count))
names(mut_loc_count_adeno) <- tolower(names(mut_loc_count_adeno))


### summarize duplicated line
## mut_loc_count: same positions with multiple alt allele are split in multiple line
## coding_info: in addition to that, mutations affecting multiple ORFs are split in multiple line

## mut_loc_count
mut_loc_count<-
mut_loc_count %>%
  group_by(start,
           end,
           in_orf,
           out_orf,
           in_orf_no_tr,
           mutation_type)%>%
  summarise(total_count = sum(count))
mut_loc_count_adeno<-
  mut_loc_count_adeno %>%
  group_by(start,
           end,
           in_orf,
           out_orf,
           in_orf_no_tr,
           mutation_type)%>%
  summarise(total_count = sum(count))

## coding_info
coding_info<-
coding_info %>%
  group_by(start,
           end,
           names)%>%
  mutate(GENEID = paste0(GENEID, collapse = ";"),
         CONSEQUENCE = paste0(CONSEQUENCE, collapse = ";"))


### merge
coding_info <- left_join(coding_info,mut_loc_count,by=c("start","end"))
# 
names(mut_loc_count_adeno)[length(names(mut_loc_count_adeno))] <- "total_count_adeno"
coding_info <- left_join(coding_info,mut_loc_count_adeno)





## check 
# any 0 count?
coding_info[coding_info$total_count==0,]
# Note: nope! good
# total_count >= total_count_adeno
coding_info[coding_info$total_count<coding_info$total_count_adeno,]
# Note: nope. good





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### 6: identify truncation #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# reference: home/note/hmp01/hmp01.patient_mutation_20240715.r


### check the what is the consequence of frameshift
to_check<-coding[coding$CONSEQUENCE=="frameshift"]
to_check_df<-data.frame(to_check)


# prep the DNA sequences (include 100 bp downstream seq in case there's elongation)
c16riboseqorf141_DNA_exon1 <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, "chr16", 89716582, 89716637))
c16riboseqorf141_DNA_exon2 <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, "chr16", 89717344, 89717386+100))
c16riboseqorf141_DNA <- paste(c16riboseqorf141_DNA_exon1,c16riboseqorf141_DNA_exon2,sep="")
c16riboseqorf142_DNA <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, "chr16", 89716634, 89716891+100))
c16riboseqorf143_DNA <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, "chr16", 89717451, 89718032+100))
#
wt_seq <-
  data.frame(
    orf_name = c("c16riboseqorf141", "c16riboseqorf142", "c16riboseqorf143"),
    wt_seq = c(
      c16riboseqorf141_DNA,
      c16riboseqorf142_DNA,
      c16riboseqorf143_DNA
    )
  )

## make mutant seq
for (i in 1:nrow(to_check_df))
{
  # get seq
  mut_seq <- wt_seq[wt_seq$orf_name==to_check_df$GENEID[i],"wt_seq"]
  # replace some char (nucleotides) in the string the
  stringr::str_sub(mut_seq,to_check_df$CDSLOC.start[i],to_check_df$CDSLOC.end[i]) <- toString(to_check_df$ALT[[i]])
  # translate
  mut_seq_aa <- translate(DNAString(mut_seq))
  to_check_df$mut_seq[i] <- as.character(mut_seq_aa)
}

## mut cds length
# Find the position of the first occurrence of "*"
position <- str_locate(to_check_df$mut_seq, pattern='\\*')
to_check_df$mut_cds_length <-  position[,1]*3

## wt cds length
to_check_df[to_check_df$GENEID=="c16riboseqorf142","wt_cds_length"] <-258
to_check_df[to_check_df$GENEID=="c16riboseqorf143","wt_cds_length"] <-582

## merge with other info
coding_info <- left_join(coding_info,to_check_df[,c("GENEID","start","end","mut_cds_length", "wt_cds_length")])


### nonsense mutation
coding[coding$CONSEQUENCE=="nonsense"]

# GRanges object with 2 ranges and 17 metadata columns:
#   seqnames    ranges strand | paramRangeID            REF                ALT      QUAL      FILTER
# <Rle> <IRanges>  <Rle> |     <factor> <DNAStringSet> <DNAStringSetList> <numeric> <character>
#   chr16:89716700_C/T    chr16  89716700      + |           NA              C                  T        NA           .
# chr16:89717784_C/T    chr16  89717784      + |           NA              C                  T        NA           .
# varAllele    CDSLOC    PROTEINLOC   QUERYID        TXID         CDSID           GENEID CONSEQUENCE
# <DNAStringSet> <IRanges> <IntegerList> <integer> <character> <IntegerList>      <character>    <factor>
#   chr16:89716700_C/T              T        67            23         9           2             2 c16riboseqorf142    nonsense
# chr16:89717784_C/T              T       334           112        53           3             5 c16riboseqorf143    nonsense
# REFCODON       VARCODON         REFAA         VARAA
# <DNAStringSet> <DNAStringSet> <AAStringSet> <AAStringSet>
#   chr16:89716700_C/T            CAG            TAG             Q             *
#   chr16:89717784_C/T            CGA            TGA             R             *
#   -------
#   seqinfo: 2779 sequences from an unspecified genome

## both alt allele are at the first base of the stop codon -> pos of mutation +2 is the end of mutatant cds


## chr16:89716700_C/T :
(89716700+2)-89716634+1
# [1] 69

## chr16:89717784_C/T
(89717784+2)-89717451+1
# [1] 336

## store info
coding_info[coding_info$CONSEQUENCE=="nonsense"&
            coding_info$GENEID=="c16riboseqorf143","mut_cds_length"] <- 336
coding_info[coding_info$CONSEQUENCE=="nonsense"&
              coding_info$GENEID=="c16riboseqorf142","mut_cds_length"] <- 69




### export
write.csv(coding_info,
          sprintf("%s/gdc_coding_mutation_info.csv",out_path),
          quote=F,
          row.names=F)








#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### 7: do mutaitons overlap binding sites #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

### make range object for feature: binding sites
bs_df <- data.frame(chr="chr16", 
                    bs=c("miR-184",
                         "miR-184",
                         "miR-4739",
                         "miR-4739",
                         "miR-4739",
                         "miR-4739",
                         "miR-4739",
                         "miR-491-5p",
                         "miR-214-3p",
                         "miR-214-3p",
                         "miR-214-3p",
                         "miR-214-3p",
                         "miR-361-3p",
                         "miR-361-3p",
                         "miR-361-3p",
                         "miR-361-3p",
                         "miR-377-3p",
                         "miR-187-3p",
                         "miR-324-5p",
                         "TR"),
                       start=c(89712141,
                               89712158,
                               89717543,
                               89717550,
                               89717553,
                               89717558,
                               89717561,
                               89717660,
                               89717551,
                               89717557,
                               89717563,
                               89717565,
                               89716563,
                               89716566,
                               89716571,
                               89716575,
                               89717836,
                               89717671,
                               89716690,
                               89717829), 
                       end=c(89712151,
                             89712164,
                             89717547,
                             89717551,
                             89717554,
                             89717559,
                             89717567,
                             89717666,
                             89717555,
                             89717557,
                             89717563,
                             89717571,
                             89716564,
                             89716566,
                             89716572,
                             89716580,
                             89717842,
                             89717679,
                             89716697,
                             89718127),
                       strand=".")


bs_gr <- makeGRangesFromDataFrame(bs_df[,c("chr","start","end","strand")])
names(bs_gr) <- bs_df$bs

### check overlap
PH01 <- findOverlaps(bs_gr,rowRanges(vcf), type = "any")


###
write.csv(data.frame(mutation=names(rowRanges(vcf))[subjectHits(PH01)],
                     feature = names(bs_gr)[queryHits(PH01)]),
          sprintf("%s/gdc_bs_mutation_info.csv",out_path),
          quote=F,
          row.names=F)