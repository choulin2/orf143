#transcript_vis.r
# Lin Chou
#2023.10.14
### background:
#(1) 
### aims:
#(1)  visualize VPS9D1-AS1 transcripts
### conclusion:
#(1) 
### reference:
#(1) 





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### 1: INITIALIZE #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


### library
library("tidyr")
library("dplyr")
library("ggplot2")
library("stringr")



### read data
coord_bs <- read.csv("transcript_features.csv")


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
##### 2: manually organize data #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


### format
coord_bs<-
  coord_bs %>%
  dplyr::select(id,strand,substrate,cancer,reference,hg38_start,hg38_end)
names(coord_bs)<-c("id","strand","substrate","cancer","reference","start","end")
# add feature name
coord_bs$name <- coord_bs$substrate
coord_bs[coord_bs$name=="","name"] <- coord_bs[coord_bs$name=="","id"]
# width of features
coord_bs$length <-(coord_bs$end-coord_bs$start)+1
coord_bs[coord_bs$strand=="","strand"]<- "*"
# remove MYU
coord_annotation <- coord_bs%>%filter(name!="MYU")

## separate features into two panels
# BS
coord_anno_bs <- coord_annotation %>% filter(strand!="+")%>%filter(id!="TR")
# orfs
coord_anno_orf <- coord_annotation %>% filter(strand=="+")
#
coord_anno_TR <- coord_annotation %>% filter(id=="TR")




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### 3: prep for gene model #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

### set up range
MYU_start <- 89711856
MYU_end <- 89718165
# https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000261373;r=16:89711856-89718165;t=ENST00000562866


### make gene models
library(biomaRt)
library(Gviz)
library(GenomicRanges)
ensembl = useMart( #host="aug2020.archive.ensembl.org",
  biomart='ENSEMBL_MART_ENSEMBL'    ,
  dataset='hsapiens_gene_ensembl'   )
## make tracks
itrack  = IdeogramTrack(chromosome=16, "hg19")
gtrack  = GenomeAxisTrack()
getrack = BiomartGeneRegionTrack(start   = MYU_start,
                                 end     = MYU_end,
                                 biomart = ensembl,
                                 chrom   = "chr16",
                                 genome  = "hg38",
                                 name    = "genes",
                                 geneSymbol = TRUE
                                 # ,filters=list(hgnc_symbol="VPS9D1-AS1" )
)

# BS
aTrack_bs <- AnnotationTrack(start = coord_anno_bs$start, width = coord_anno_bs$length, 
                             chromosome = "chr16", 
                             strand = coord_anno_bs$strand,
                             group = coord_anno_bs$name,
                             genome = "hg38", name = "binding sites")
# orfs
aTrack_orf <- AnnotationTrack(start = coord_anno_orf$start, width = coord_anno_orf$length, 
                              chromosome = "chr16", 
                              strand = coord_anno_orf$strand,
                              group = coord_anno_orf$name,
                              genome = "hg38", name = "nORFs")
# repeat
aTrack_repeat <- AnnotationTrack(start = coord_anno_TR$start, width = coord_anno_TR$length,
                                 chromosome = "chr16",
                                 strand = c("*"),
                                 group = c("Repeats"),
                                 genome = "hg38", name = "TR")





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### 4: make gene model #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


### overview of the entire MYU
pdf(file = sprintf("%s/gene_model_MYU.pdf",fig_path),
    width = one.c*0.0393701, # The width of the plot in inches
    height = one.c*0.0393701) # The height of the plot in inches
# plot
plotTracks(list(itrack,
                # sTrack,
                gtrack,
                getrack,
                aTrack_bs,
                aTrack_orf,
                aTrack_repeat
                ),
           type="l",
           groupAnnotation = "group",
           transcriptAnnotation  = "symbol",
           # noLetters = T,
           chromosome = 16,
           from=MYU_start-10,to=MYU_end+10
)
dev.off()