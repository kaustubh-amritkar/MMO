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
library(tidyr)
library(scales)
library(ggnewscale)

### Read structural phylogeny + Map the taxa and/or major_group info on it ###
phy <- read.tree(file = "/Volumes/bkacar/Kaustubh/MMO/data/foldseek_search/hit_AF_structures_database/TM_score_above_5_structural_comparison/IQ-Tree/different_rootings/MMO_filtered_foldmason_MSA_aa.fa.MAD_rooted.treefile")

df = data.frame(read.delim("/Volumes/bkacar/Kaustubh/MMO/data/foldseek_search/hit_AF_structures_database/TM_score_above_5_structural_comparison/filtered_hits_all_info_df.tsv"))

annot_df <- df %>% select(target, L2, major_group, superkingdom, phylum, class) %>% rename(label = target)

p <- ggtree(phy) %<+% annot_df + geom_treescale()
  
#p1 = p + geom_tree(aes(color = major_group))
p1 <- p + geom_tippoint(aes(color=major_group), size = 3) 

### Map the length of the hits as a heatmap ###
l2_df <- annot_df %>% select(label, L2)
rownames(l2_df) <- l2_df$label
l2_df$label <- NULL

p2 <- gheatmap(p1, l2_df, offset=0.5, width = 0.05, colnames_position = "top") + scale_fill_gradient(low = "white", high = "blue", limits=c(200,600), oob = squish, name = "L2")

### Read fasta for active site residues and map them on the phylogeny as a heatmap ###
MMO_active_site_foldmason_MSA <- readAAStringSet("/Volumes/bkacar/Kaustubh/MMO/data/foldseek_search/hit_AF_structures_database/TM_score_above_5_structural_comparison/Active_site_analysis/MMO_filtered_active_site_residue_MSA_from_foldmason.fasta")
MMO_active_site_foldmason_MSA_df <- as.data.frame(as.matrix(MMO_active_site_foldmason_MSA))
MMO_active_site_foldmason_MSA_df$label <- rownames(MMO_active_site_foldmason_MSA_df)
rownames(MMO_active_site_foldmason_MSA_df) <- MMO_active_site_foldmason_MSA_df$label
MMO_active_site_foldmason_MSA_df$label <- NULL

p2 <- p2 + new_scale_fill()
p3 <- gheatmap(p2, MMO_active_site_foldmason_MSA_df, offset = 2.5, width = 0.2, colnames_position = "top") + scale_fill_manual(
  values = c(
    "-" = "gray90", 
    "E" = "darkred", "D" = "#f25c54", 
    "H" = "darkblue",
    "W" = "purple", 
    "Q" = "#b43e8f"
  ),
  na.value = "gray60",
  breaks = c("-", "E", "H", "W", "Q", "D")
)

### Save the figure ### 
ggsave("/Volumes/bkacar/Kaustubh/MMO/Figures/MMO_X_foldseek_phylogeny_with_length_residue_conservation_type_hg.svg", plot = p3)

heatmap_df <- annot_df[, c("label", "major_group")]
rownames(heatmap_df) <- heatmap_df$label
heatmap_df$label <- NULL
gheatmap(p2, heatmap_df, offset = 0.5, 
         width = 0.3,
         colnames_position = "top",
         colnames_angle = 90)


seq_phy <- read.tree(file = "/Volumes/bkacar/Kaustubh/MMO/data/blast_search/IQTree/different_rootings/MMO_blast_output_seqs.c0.9.mafft.fasta.MAD_rooted.treefile")
blastp_result_df = data.frame(read.delim("/Volumes/bkacar/Kaustubh/MMO/data/blast_search/MMO_blastp_output_after_clustering_information.tsv"))
blastp_annot_df <- blastp_result_df %>% select(id, major_group) %>% rename(label = id)

bp <- ggtree(seq_phy) %<+% blastp_annot_df + geom_treescale()
bp1 = bp + geom_tippoint(aes(color=major_group), size = 3)
