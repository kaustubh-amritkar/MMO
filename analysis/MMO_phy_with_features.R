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

phy <- read.tree(file = "/Volumes/bkacar/Kaustubh/MMO/data/foldseek_search/hit_AF_structures_database/TM_score_above_5_structural_comparison/IQ-Tree/different_rootings/MMO_filtered_foldmason_MSA_aa.fa.MAD_rooted.treefile")

df = data.frame(read.delim("/Volumes/bkacar/Kaustubh/MMO/data/foldseek_search/hit_AF_structures_database/TM_score_above_5_structural_comparison/filtered_hits_all_info_df.tsv"))

annot_df <- df %>% select(target, L2, major_group, superkingdom, phylum, class) %>% rename(label = target)

p <- ggtree(phy) %<+% annot_df + geom_treescale()
  
p1 = p + geom_tippoint(aes(color=superkingdom), size = 3) 
gheatmap(p1, annot_df[,c("major_group")], offset = 0.5, 
         width = 0.3,
         colnames_position = "top",
         colnames_angle = 90)


seq_phy <- read.tree(file = "/Volumes/bkacar/Kaustubh/MMO/data/blast_search/IQTree/different_rootings/MMO_blast_output_seqs.c0.9.mafft.fasta.MAD_rooted.treefile")
blastp_result_df = data.frame(read.delim("/Volumes/bkacar/Kaustubh/MMO/data/blast_search/MMO_blastp_output_after_clustering_information.tsv"))
blastp_annot_df <- blastp_result_df %>% select(id, major_group) %>% rename(label = id)

bp <- ggtree(seq_phy) %<+% blastp_annot_df + geom_treescale()
bp1 = bp + geom_tippoint(aes(color=major_group), size = 3)
