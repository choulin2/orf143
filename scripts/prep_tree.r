#prep_tree.r
# Lin Chou
#2023.4.21
### aims:
#(1) trim to species tree from Ensembl





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#####1: prep a guide tree #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


## download the species tree from ENSEMBL
#https://github.com/Ensembl/ensembl-compara/blob/release/104/conf/vertebrates/species_tree.branch_len.nw
# save as: ensembl_vert_species.nw


### library
library(ape)
library(ggplot2)


### read tree
vert_species_tree <- ape::read.tree("ensembl_vert_species.nw")
plot(vert_species_tree)

## number of nodes
vert_species_tree$Nnode
# [1] 226





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#####2: trim the tree to relevant species#####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


### select the 26 species

## use the curated alignment file
# /Users/choulin/scratch/hmp01/hmp01.04/realingn/1kflank.nt.fasta
oneK_species <- c("dasypus_novemcinctus",
                  "jaculus_jaculus",
                  "myotis_lucifugus",
                  "chinchilla_lanigera",
                  "saimiri_boliviensis_boliviensis",
                  "octodon_degus",
                  "dipodomys_ordii",
                  "heterocephalus_glaber_female",
                  "carlito_syrichta",
                  "ictidomys_tridecemlineatus",
                  "urocitellus_parryii",
                  "callithrix_jacchus",
                  "aotus_nancymaae",
                  "gorilla_gorilla",
                  "rhinopithecus_bieti",
                  "rhinopithecus_roxellana",
                  "macaca_nemestrina",
                  "cercocebus_atys",
                  "papio_anubis",
                  "pongo_abelii",
                  "homo_sapiens",
                  "pan_paniscus",
                  "pan_troglodytes",
                  "otolemur_garnettii",
                  "prolemur_simus",
                  "propithecus_coquereli")

## keep tips of interest
plot(keep.tip(vert_species_tree, oneK_species))
write.tree(keep.tip(vert_species_tree, oneK_species))
# (dasypus_novemcinctus:0.1043134,((((((chinchilla_lanigera:0.07290589,heterocephalus_glaber_female:0.0728056):0.00217621,octodon_degus:0.0761923):0.0218854,(ictidomys_tridecemlineatus:0.0138367,urocitellus_parryii:0.0127433):0.07663759):0.00726934,(dipodomys_ordii:0.100795,jaculus_jaculus:0.10254806):0.0001):0.00481115,(((prolemur_simus:0.0371715,propithecus_coquereli:0.0381785):0.03257458,otolemur_garnettii:0.0749096):0.00718259,(carlito_syrichta:0.0797382,(((saimiri_boliviensis_boliviensis:0.02906148,callithrix_jacchus:0.0300164):0.0001,aotus_nancymaae:0.0252751):0.0236248,(((((pan_paniscus:0.00332162,pan_troglodytes:0.00221838):0.00431454,homo_sapiens:0.00659546):0.00185346,gorilla_gorilla:0.00857503):0.00839492,pongo_abelii:0.0171801):0.01370794,((rhinopithecus_roxellana:0.00209318,rhinopithecus_bieti:0.00301682):0.0143382,(macaca_nemestrina:0.00875003,(papio_anubis:0.00663729,cercocebus_atys:0.006350422):0.00113664):0.00871929):0.013672):0.0176123):0.0277667):3.95857e-05):0.021424):0.0001,myotis_lucifugus:0.10948728):0.0001);


### select the 5 ape species for ASR

## use the curated alignment file
# /Users/choulin/scratch/hmp01/hmp01.04/realingn/1kflank.nt.fasta
oneK_species <- c("pongo_abelii",
                  "homo_sapiens",
                  "pan_paniscus",
                  "pan_troglodytes",
                  "gorilla_gorilla")

## keep tips of interest
plot(keep.tip(vert_species_tree, oneK_species))
write.tree(keep.tip(vert_species_tree, oneK_species))
# [1] "((((pan_paniscus:0.00332162,pan_troglodytes:0.00221838):0.00431454,homo_sapiens:0.00659546):0.00185346,gorilla_gorilla:0.00857503):0.00839492,pongo_abelii:0.0171801);"