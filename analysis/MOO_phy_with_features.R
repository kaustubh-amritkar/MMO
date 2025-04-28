if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

library(Biostrings)
library(ggtree)
library(ggimage)
library(ape)
library(ggnewscale)
library(svglite)
library(viridis)
library(dplyr)

phy <- read.tree(file = '/Volumes/bkacar/Kaustubh/MOO/data/foldseek_search/hit_AF_structures_database/TM_score_above_5_structural_comparison/IQ-Tree/MOO_filtered_foldmason_MSA_aa.fa.treefile')

df = data.frame(read.csv("/Volumes/bkacar/Kaustubh/MOO/data/foldseek_search/hit_AF_structures_database/TM_score_above_5_structural_comparison/filtered_hits_all_info_df.csv"))

annot_df <- df %>% select(target, L2, major_group, superkingdom, phylum, class) %>% rename(label = target)

p <- ggtree(phy, layout = "unrooted") %<+% annot_df + geom_treescale()
  
p + geom_tippoint(aes(color=major_group), size = 2)
