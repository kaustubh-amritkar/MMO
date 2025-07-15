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

phy <- read.tree(file = "/Volumes/bkacar/Kaustubh/MOO/data/foldseek_search/hit_AF_structures_database/TM_score_above_5_structural_comparison/IQ-Tree/different_rootings/MOO_filtered_foldmason_MSA_aa.fa.MAD_rooted.treefile")

df = data.frame(read.delim("/Volumes/bkacar/Kaustubh/MOO/data/foldseek_search/hit_AF_structures_database/TM_score_above_5_structural_comparison/filtered_hits_all_info_df.tsv"))

annot_df <- df %>% select(target, L2, major_group, superkingdom, phylum, class) %>% rename(label = target)

p <- ggtree(phy) %<+% annot_df + geom_treescale()
  
p1 = p + geom_tippoint(aes(color=superkingdom), size = 3) 
gheatmap(p1, annot_df[,c("major_group")], offset = 0.5, 
         width = 0.3,
         colnames_position = "top",
         colnames_angle = 90)
