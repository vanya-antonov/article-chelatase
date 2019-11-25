
source('lib.R')

###

tree <- read_tree_from_fn('ref_tree_and_all_fs_chlD.tree')
anno_df <- read_anno_for_tree("genes_chel.txt", tree)
length(tree$tip.label) == nrow(anno_df)

# Re-root the tree
outgroup_genes <- anno_df %>% filter(kingdom == 'Archaea' & fs_len == 'None') %>% pull(name)
tree <- root(tree, outgroup = outgroup_genes)
# ggtree(tree)

# Extract the REF_TREE only
ref_genes <- anno_df %>% filter(!is.na(gene)) %>% pull(name)
ref_tree <- keep.tip(tree, ref_genes)
# ggtree(ref_tree)

# Annotate branches
chl_tips <- anno_df %>% filter(gene %in% c('chlD', 'bchD'), phylum == 'Cyanobacteria') %>% pull(name)
bch_tips <- anno_df %>% filter(gene == 'bchD', phylum == 'Proteobacteria') %>% pull(name)
cob_chel_tips <- anno_df %>% filter(gene == 'cobT') %>% pull(name)

my_tree <- ref_tree
title <- sprintf('Total number of genes in the tree = %s', length(my_tree$tip.label))
ggtree(my_tree) %<+% anno_df +
  geom_cladelabel(node=findMRCA(my_tree, chl_tips), label="chlD", color='darkgreen', offset.text=.1) +
  geom_hilight(node=findMRCA(my_tree, chl_tips), fill="green", alpha=.3) +
  geom_cladelabel(node=findMRCA(my_tree, bch_tips), label="bchD", color='red', offset.text=.1) +
  geom_hilight(node=findMRCA(my_tree, bch_tips), fill="red", alpha=.3) +
  geom_cladelabel(node=findMRCA(my_tree, cob_chel_tips), label="cobT", color='blue', angle=270, hjust='center', offset.text=.1) +
  geom_hilight(node=findMRCA(my_tree, cob_chel_tips), fill="blue", alpha=.6) +
  ggtitle(title)
ggsave('ref_tree.pdf', path = paste0(OUT_DIR, 'ref_tree'), width = 4, height = 6)
