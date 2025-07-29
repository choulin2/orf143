#why_orf143.r

#2025.3.12
### aims:
#(1)  show why we chose orf143:
# genome-wdie id of TR and comparison with nORFs.
# check length of orf143





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### 1: INITIALIZE #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


### library
library("RMariaDB")
library("tidyr")
library("plyr")
library("dplyr")
library("ggplot2")
library("ggrepel")
library("stringr")
library("GenomicRanges")


### read ORF data
nORF_coord <- read.csv("/home/lic130/home/data/hmp01/nORF_coord", sep = "\t")


### read de novo nORF based on Mudge 2022/ Sandmann 2023 (two papers have same set of ORFs)
# read file
# https://www.cell.com/molecular-cell/fulltext/S1097-2765(23)00075-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1097276523000758%3Fshowall%3Dtrue#supplementaryMaterial
evo <- read.csv("Sandmann_2023_table_s1.csv")
names(evo)[1] <- "orf_name"
nrow(evo)


### read data
TR_raw<-NULL
for (i in c(1:22,"X","Y")){
  # read table
  my_table<-read.csv(sprintf("/TRF_results/chr%s.fasta.2.7.7.80.10.50.500.dat",i),skip = 15,sep=" ",header = F)
  # assign chr
  my_table$chr <- i
  # assemble
  TR_raw <- rbind(TR_raw, my_table)
  print(i)
}

# names
names(TR_raw) <- c("start",
                   "end",
                   "period_size",
                   "copy_number",
                   "consensus_size",
                   "percent_matches",
                   "percent_indels",
                   "score",
                   "A",
                   "C",
                   "G",
                   "T",
                   "entropy",
                   "seq_consensus",
                   "seq_full",
                   "chr")


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
##### 2: make ranges #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


### nORF

## exclude introns by separate collapsed columns into multiple rows
nORF_coord<-separate_rows(nORF_coord, starts, ends, convert = TRUE)
sum(grepl(";",nORF_coord$starts))

## add evo info
nORF_coord <-
  left_join(nORF_coord, evo[,c("orf_name", "Conservation.ORF", "Modes.of.evolution")])

## granges
nORF_range <-makeGRangesFromDataFrame(nORF_coord,
                                      seqnames.field="chrm",
                                      start.field="starts",
                                      end.field="ends",
                                      keep.extra.columns=T)




### grnages of TR
# make ranges
TR_range <-makeGRangesFromDataFrame(TR_raw,
                                    seqnames.field="chr",
                                    start.field="start",
                                    end.field="end")
# check
head(TR_range)
unique(TR_range@seqnames)
length(TR_range)
# merge overlapping ranges
TR_range <- reduce(TR_range)





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### 3:nORFs that overlap with TR #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


### number of ORF overlaps with TR

## make a function
get_overlap <- function(my_ranges){
  ## total length
  length_total_ORF <- sum(width(my_ranges))
  
  ## total number of nt overlapping with TR
  # find overlap
  hits <-findOverlaps(query = my_ranges,
                      subject = TR_range,
                      type= "any",
                      select = "all",
                      ignore.strand=T
  )
  overlaps <- pintersect(my_ranges[queryHits(hits)], 
                         TR_range[subjectHits(hits)],
                         ignore.strand=T)
  return(overlaps)}


## nORF
overlap <- get_overlap(nORF_range)
overlap <-
  overlap[, c("orf_name", "Modes.of.evolution")] %>%
  as.data.frame()
# summarize
overlap <-
  overlap %>%
  group_by(orf_name)%>%
  # total length overlapping with TR (some orfs have intron)
  dplyr::summarise(total_TR_overlap=sum(width))
# number of ORFs that overlap with TR
nrow(overlap)
# [1] 168


### plot the nORFs that overlap with TR

## data
overlap <- get_overlap(nORF_range)
overlap <-
  overlap[, c("orf_name", "Modes.of.evolution")] %>%
  as.data.frame()
# summarize
overlap <- overlap %>%
  group_by(orf_name,Modes.of.evolution)%>%
  # total length overlapping with TR (some orfs have intron)
  dplyr::summarise(total_TR_overlap=sum(width))%>%
  arrange(desc(total_TR_overlap))
overlap$index <- 1:nrow(overlap)
overlap[grepl("norep",overlap$orf_name),"group"] <- "single_study"
overlap[grepl("riboseqorf",overlap$orf_name),"group"] <- "multiple_studies"

## length of CDS
nORF_length <-
  as.data.frame(nORF_range) %>% group_by(orf_name) %>% dplyr::summarise(orf_length = sum(width))
overlap <- left_join(overlap,nORF_length,by="orf_name")
sum((overlap$total_TR_overlap/overlap$orf_length)>1)
# [1] 0

## plot
p <-
  overlap %>%
  ggplot(aes(x=index,y=total_TR_overlap,color=group,label=orf_name))+
  geom_point()+
  geom_label_repel(data         = subset(overlap, total_TR_overlap>200),
                   nudge_x       = subset(overlap, total_TR_overlap>200)$index*20, # x position
  )+
  ggtitle("nORFs that overlap with TR") +
  xlab("arbitrary index") + 
  ylab("Size of overlap (bp)")+
  scale_color_manual("Translation detected in",values=c("black","grey"))+
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text.y = element_text(size = txt),
        axis.text.x = element_text(size = txt, angle = 0),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))

p
## export
pdf(file = sprintf("%s/nORF_with_TR.pdf",fig_path),
    width = one.c*0.0393701, # The width of the plot in inches
    height = one.c*0.0393701) # The height of the plot in inches
# plot
p
dev.off()





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### 4: how does the length of orf143 compare to others? #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


### descriptive analysis
summary(nORF_length$orf_length)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.0    72.0   105.0   133.1   165.0  5904.0 

# 4 times longer than mean
582/133.1
# [1] 4.372652


## percentile rank
nORF_length %>%
  arrange(desc(orf_length))%>%
  mutate(percent_rank = rank(orf_length)/length(orf_length),
         Rank = 1:nrow(nORF_length))%>%
  head(n=20)

# # A tibble: 20 × 4
# orf_name         orf_length percent_rank  Rank
# <chr>                 <int>        <dbl> <int>
#   1 c6norep158             5904        1         1
# 2 c11riboseqorf108       1050        1.00      2
# 3 c5norep107              822        1.00      3
# 4 c17riboseqorf88         741        1.00      4
# 5 c1riboseqorf57          726        0.999     5
# 6 c6norep53               681        0.999     6
# 7 c19norep137             660        0.999     7
# 8 c19riboseqorf4          630        0.999     8
# 9 c15norep24              624        0.999     9
# 10 c17norep138             624        0.999    10
# 11 c7riboseqorf164         618        0.999    11
# 12 c13riboseqorf11         609        0.998    12
# 13 c12norep105             597        0.998    13
# 14 c6riboseqorf8           597        0.998    14
# 15 c14riboseqorf26         594        0.998    15
# 16 c5norep173              594        0.998    16
# 17 c16riboseqorf59         585        0.998    17
# 18 c16riboseqorf143        582        0.998    18
# 19 c5norep205              582        0.998    19
# 20 c8norep145              579        0.997    20