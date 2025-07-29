#cancer_stats.r
# Lin Chou
#2025.7.17
### aims:
#(1) testing whether there are more mutation than expected in ORFs.
#(2) use the counts from cancer_count_indel.ipynb and cancer_count_snp.ipynb





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### 1: INITIALIZE #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


### library
library("GenomicRanges")
library("tidyr")
library("dplyr")
library("ggplot2")
library("stringr")
library("vcfR")




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### 2: are exonic mutations biased toward ORFs? #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


### set up lengths
exon_ORF_length = 948
exon_noORF_length = 805
exon_noTR_ORF_length = 744
exon_noTR_noORF_length =710

### GDC+HMF
# SNV
binom.test(67+35,67+3+35+28,p=exon_ORF_length/(exon_ORF_length+exon_noORF_length))

# Exact binomial test
# 
# data:  67 + 35 and 67 + 3 + 35 + 28
# number of successes = 102, number of trials = 133, p-value = 1.077e-07
# alternative hypothesis: true probability of success is not equal to 0.5407872
# 95 percent confidence interval:
#   0.6857673 0.8358066
# sample estimates:
#   probability of success 
# 0.7669173 

# INDEL
binom.test(15+8,15+0+8+0,p=exon_ORF_length/(exon_ORF_length+exon_noORF_length))

# Exact binomial test
# 
# data:  15 + 8 and 15 + 0 + 8 + 0
# number of successes = 23, number of trials = 23, p-value = 1.197e-06
# alternative hypothesis: true probability of success is not equal to 0.5407872
# 95 percent confidence interval:
#   0.8518149 1.0000000
# sample estimates:
#   probability of success 
# 1 

### GDC
# SNV
binom.test(67,67+3,p=exon_ORF_length/(exon_ORF_length+exon_noORF_length))

# Exact binomial test
# 
# data:  67 and 67 + 3
# number of successes = 67, number of trials = 70, p-value = 8.073e-15
# alternative hypothesis: true probability of success is not equal to 0.5407872
# 95 percent confidence interval:
#   0.8798214 0.9910731
# sample estimates:
#   probability of success 
# 0.9571429 

# INDEL
binom.test(15,15+0,p=exon_ORF_length/(exon_ORF_length+exon_noORF_length))

# Exact binomial test
# 
# data:  15 and 15 + 0
# number of successes = 15, number of trials = 15, p-value = 0.0001075
# alternative hypothesis: true probability of success is not equal to 0.5407872
# 95 percent confidence interval:
#   0.7819806 1.0000000
# sample estimates:
#   probability of success 
# 1 

## Adenomas and Adenocarcinomas
# SNV
binom.test(34,34+2,p=exon_ORF_length/(exon_ORF_length+exon_noORF_length))

# Exact binomial test
# 
# data:  34 and 34 + 2
# number of successes = 34, number of trials = 36, p-value = 2.046e-07
# alternative hypothesis: true probability of success is not equal to 0.5407872
# 95 percent confidence interval:
#   0.8133633 0.9931997
# sample estimates:
#   probability of success 
# 0.9444444 

# INDEL
binom.test(12,12+0,p=exon_ORF_length/(exon_ORF_length+exon_noORF_length))

# Exact binomial test
# 
# data:  12 and 12 + 0
# number of successes = 12, number of trials = 12, p-value = 0.0007136
# alternative hypothesis: true probability of success is not equal to 0.5407872
# 95 percent confidence interval:
#   0.7353515 1.0000000
# sample estimates:
#   probability of success 
# 1

### HMF
# SNV
binom.test(35,35+28,p=exon_ORF_length/(exon_ORF_length+exon_noORF_length))

# Exact binomial test
# 
# data:  35 and 35 + 28
# number of successes = 35, number of trials = 63, p-value = 0.8996
# alternative hypothesis: true probability of success is not equal to 0.5407872
# 95 percent confidence interval:
#   0.4248876 0.6808269
# sample estimates:
#   probability of success 
# 0.5555556 

# INDEL
binom.test(8,8+0,p=exon_ORF_length/(exon_ORF_length+exon_noORF_length))

# Exact binomial test
# 
# data:  8 and 8 + 0
# number of successes = 8, number of trials = 8, p-value = 0.009292
# alternative hypothesis: true probability of success is not equal to 0.5407872
# 95 percent confidence interval:
#   0.6305834 1.0000000
# sample estimates:
#   probability of success 
# 1 





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### 3: are exonic mutations biased toward ORFs? - no TR #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


### set up
exon_ORF_length = 948
exon_noORF_length = 805
exon_noTR_ORF_length = 744
exon_noTR_noORF_length =710

### GDC+HMF
# SNV
binom.test(48+23,48+0+23+27,p=exon_noTR_ORF_length/(exon_noTR_ORF_length+exon_noTR_noORF_length))

# Exact binomial test
# 
# data:  48 + 23 and 48 + 0 + 23 + 27
# number of successes = 71, number of trials = 98, p-value = 2.634e-05
# alternative hypothesis: true probability of success is not equal to 0.5116919
# 95 percent confidence interval:
#   0.6250085 0.8099460
# sample estimates:
#   probability of success 
# 0.7244898 

# INDEL
binom.test(10+5,10+0+5+0,p=exon_noTR_ORF_length/(exon_noTR_ORF_length+exon_noTR_noORF_length))

# Exact binomial test
# 
# data:  10 + 5 and 10 + 0 + 5 + 0
# number of successes = 15, number of trials = 15, p-value = 6.456e-05
# alternative hypothesis: true probability of success is not equal to 0.5116919
# 95 percent confidence interval:
#   0.7819806 1.0000000
# sample estimates:
#   probability of success 
# 1 


### GDC
# SNV
binom.test(48,48+0,p=exon_noTR_ORF_length/(exon_noTR_ORF_length+exon_noTR_noORF_length))

# Exact binomial test
# 
# data:  48 and 48 + 0
# number of successes = 48, number of trials = 48, p-value = 1.192e-14
# alternative hypothesis: true probability of success is not equal to 0.5116919
# 95 percent confidence interval:
#   0.9260272 1.0000000
# sample estimates:
#   probability of success 
# 1 

# INDEL
binom.test(10,10+0,p=exon_noTR_ORF_length/(exon_noTR_ORF_length+exon_noTR_noORF_length))

# Exact binomial test
# 
# data:  10 and 10 + 0
# number of successes = 10, number of trials = 10, p-value = 0.002001
# alternative hypothesis: true probability of success is not equal to 0.5116919
# 95 percent confidence interval:
#   0.6915029 1.0000000
# sample estimates:
#   probability of success 
# 1 

## Adenomas and Adenocarcinomas
# SNV
binom.test(25,25+0,p=exon_noTR_ORF_length/(exon_noTR_ORF_length+exon_noTR_noORF_length))

# Exact binomial test
# 
# data:  25 and 25 + 0
# number of successes = 25, number of trials = 25, p-value = 6.961e-08
# alternative hypothesis: true probability of success is not equal to 0.5116919
# 95 percent confidence interval:
#   0.8628148 1.0000000
# sample estimates:
#   probability of success 
# 1 

# INDEL
binom.test(9,9+0,p=exon_noTR_ORF_length/(exon_noTR_ORF_length+exon_noTR_noORF_length))

# Exact binomial test
# 
# data:  9 and 9 + 0
# number of successes = 9, number of trials = 9, p-value = 0.003983
# alternative hypothesis: true probability of success is not equal to 0.5116919
# 95 percent confidence interval:
#   0.6637329 1.0000000
# sample estimates:
#   probability of success 
# 1 


### HMF
# SNV
binom.test(23,23+27,p=exon_noTR_ORF_length/(exon_noTR_ORF_length+exon_noTR_noORF_length))

# Exact binomial test
# 
# data:  23 and 23 + 27
# number of successes = 23, number of trials = 50, p-value = 0.4826
# alternative hypothesis: true probability of success is not equal to 0.5116919
# 95 percent confidence interval:
#   0.3181492 0.6067580
# sample estimates:
#   probability of success 
# 0.46 

# INDEL
binom.test(5,5+0,p=exon_noTR_ORF_length/(exon_noTR_ORF_length+exon_noTR_noORF_length))

# Exact binomial test
# 
# data:  5 and 5 + 0
# number of successes = 5, number of trials = 5, p-value = 0.06284
# alternative hypothesis: true probability of success is not equal to 0.5116919
# 95 percent confidence interval:
#   0.4781762 1.0000000
# sample estimates:
#   probability of success 
# 1 