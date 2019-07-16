
source("lib.R")

###

PREFIX <- 'all_proteobacteria_M'
COB_CHEL_DESCR <- c(
  'cobaltochelatase subunit CobT',
  'cobaltochelatase',
  'cobalt chelatase')
COB_LIGAND_DESCR <- c('nicotinate-nucleotide--dimethylbenzimidazole phosphoribosyltransferase')
COB_OTHER_DESCR <- c('DUF3520 domain-containing protein')
FS_COLORS <- c('-1' = 'blue', '+1' = 'red', 'None' = 'white')

###

anno_df <- read.delim(paste0(DATA_DIR, PREFIX, '.info.txt'), header = TRUE, as.is = TRUE) %>%
  mutate(fs_len = ifelse(is.na(fs_len), 'None', sprintf('%+d', fs_len))) %>%
  mutate(fs_len = factor(fs_len, levels = names(FS_COLORS))) %>%
  # There shouldn't be any shapes in genes w/o fs
  mutate(shape_alpha = ifelse(fs_len == 'None', 0, 1)) %>%
  mutate(id = as.character(id))
rownames(anno_df) <- anno_df$id
head(anno_df)

tree <- read.tree(paste0(DATA_DIR, PREFIX, '.tree'))

# Verify
length(tree$tip.label)  ==  nrow(anno_df)

# Annotate branches
chl_tips <- anno_df %>% filter(gene %in% c('chlD', 'bchD'), phylum == 'Cyanobacteria') %>% pull(id)
bch_tips <- anno_df %>% filter(gene == 'bchD', phylum == 'Proteobacteria') %>% pull(id)
cob_chel_tips <- anno_df %>% filter(descr %in% COB_CHEL_DESCR) %>% pull(id)
cob_ligand_tips <- anno_df %>% filter(descr %in% COB_LIGAND_DESCR) %>% pull(id)


get_colored_ggtree <- function(my_tree)
{
  ggtree_obj <- ggtree(my_tree) %<+% anno_df +
    # Add colored shapes to tips (according to fs_len)
    geom_tippoint(aes(fill=fs_len, alpha=shape_alpha), shape=21, stroke=0, size=2) +
    scale_fill_manual(values = FS_COLORS, name = 'Frameshift') +
    # Hide some legends: https://stackoverflow.com/a/14604540/310453
    scale_alpha_continuous(guide = FALSE)

  # Highlight available clades
  if(all(chl_tips %in% my_tree$tip.label))
  {
    ggtree_obj <- ggtree_obj +
      geom_cladelabel(node=findMRCA(my_tree, chl_tips), label="chlD (Cyanobacteria)", color='darkgreen', offset.text=.1) +
      geom_hilight(node=findMRCA(my_tree, chl_tips), fill="green", alpha=.3)
  }
  
  if(all(bch_tips %in% my_tree$tip.label))
  {
    ggtree_obj <- ggtree_obj +
      geom_cladelabel(node=findMRCA(my_tree, bch_tips), label="bchD (Proteobacteria)", color='red', offset.text=.1) +
      geom_hilight(node=findMRCA(my_tree, bch_tips), fill="red", alpha=.3)
  }
  
  if(all(cob_chel_tips %in% my_tree$tip.label))
  {
    ggtree_obj <- ggtree_obj +
      geom_cladelabel(node=findMRCA(my_tree, cob_chel_tips), label="cobT\n(chelatase)", color='blue', angle=270, hjust='center', offset.text=.25) +
      geom_hilight(node=findMRCA(my_tree, cob_chel_tips), fill="blue", alpha=.6)
  }
  
  if(all(cob_ligand_tips %in% my_tree$tip.label))
  {
    ggtree_obj <- ggtree_obj +
      geom_cladelabel(node=findMRCA(my_tree, cob_ligand_tips), label="cobT\n(ligand)", color='black', angle=270, hjust='center', offset.text=.25) +
      geom_hilight(node=findMRCA(my_tree, cob_ligand_tips), fill="black", alpha=.3)
  }
  
  ggtree_obj + theme_tree(legend.position="left")
}
#get_colored_ggtree(tree)

# I will remove the whole cob clade except the Co-chelatase subunit branch
global_cob_node <- anno_df %>%
  filter(descr %in% c(COB_CHEL_DESCR, COB_LIGAND_DESCR, COB_OTHER_DESCR)) %>%
  pull(id) %>%
  findMRCA(tree, .)

# Simulate filtering: Get the non-chelatase cob-related genes
cob_all_ids <- extract.clade(tree, global_cob_node)$tip.label
cob_chel_ids <- extract.clade(tree, findMRCA(tree, cob_chel_tips))$tip.label
cob_non_chel_ids <- setdiff(cob_all_ids, cob_chel_ids)

# Verify
length(cob_all_ids) == length(cob_chel_ids) + length(cob_non_chel_ids)

# Identify the branches of the proteins that are not chelatase-subunits
title <- sprintf('%d proteins from %d species', nrow(anno_df), length(unique(anno_df$org_name)))
subT <- sprintf('%d / %d genes with -1 / +1 frameshifts',
                nrow(filter(anno_df, fs_len == '-1')),
                nrow(filter(anno_df, fs_len == '+1')))
global_cob_title <- paste(length(cob_all_ids), "various cobT-related proteins")
get_colored_ggtree(tree) +
  geom_hilight(node=global_cob_node, fill="black", alpha=.3) +
  geom_cladelabel(node=global_cob_node, label=global_cob_title, color='black', angle=270, hjust='center', offset.text=.05) +
  ggtitle(title, subtitle = subT)
ggsave('tree_filtering_M.pdf', path = OUT_DIR)


# Simulate filtering and print all the remaining ones
chel_anno_df <- anno_df %>% filter(!id %in% cob_non_chel_ids)
sprintf('After filtering: %d chelatase subunits from %d species',
        nrow(chel_anno_df), length(unique(chel_anno_df$org_name)))

