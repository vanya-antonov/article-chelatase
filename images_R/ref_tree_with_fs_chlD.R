
source('lib.R')

###

# Genes that do not have good position in the REF_TREE
BAD_GENES <- c('AMED_7054', 'CAGG_RS15520', 'RCAS_RS13095', 'DSHI_RS17730', 'RSP_0274', 'RSPH17025_RS05185', 'BLV13_RS00120', 'SADFL11_RS19460')
FS_COLORS <- c('-1' = 'blue', '+1' = 'red', 'None' = 'white')

###

tree <- read.tree(paste0(DATA_DIR, '190821.ref_tree_with_fs_chlD.tree'))

# Remove some genes
# all(BAD_GENES %in% tree$tip.label)
tree <- drop.tip(tree, BAD_GENES)

anno_df <- read.delim(paste0(DATA_DIR, "genes_chel.txt"), as.is=TRUE) %>%
  filter(name %in% tree$tip.label) %>%
  mutate(fs_len = ifelse(is.na(fs_len), 'None', sprintf('%+d', fs_len))) %>%
  mutate(fs_len = factor(fs_len, levels = names(FS_COLORS))) %>%
  # There shouldn't be any shapes in genes w/o fs
  mutate(shape_alpha = ifelse(fs_len == 'None', 0, 1))
rownames(anno_df) <- anno_df$name
length(tree$tip.label) == nrow(anno_df)
#head(anno_df)

anno_df %>% filter(kingdom == 'Archaea')



# Re-root the tree
outgroup_genes <- anno_df %>% filter(kingdom == 'Archaea' & fs_len == 'None') %>% pull(name)
tree <- root(tree, outgroup = outgroup_genes)
# ggtree(tree)

# Annotate branches
chl_tips <- anno_df %>% filter(gene %in% c('chlD', 'bchD'), phylum == 'Cyanobacteria') %>% pull(name)
bch_tips <- anno_df %>% filter(gene == 'bchD', phylum == 'Proteobacteria') %>% pull(name)
cob_chel_tips <- anno_df %>% filter(gene == 'cobT') %>% pull(name)

title <- sprintf('Total number of genes in the tree = %s', nrow(anno_df))
subT <- sprintf('Number of fs-chlD genes = %s', anno_df %>% filter(fs_len != 'None') %>% nrow())
my_tree <- tree
ggtree(my_tree) %<+% anno_df +
#  geom_tiplab(size = 1) +
  # Add colored shapes to tips (according to fs_len)
  geom_tippoint(aes(fill=fs_len, alpha=shape_alpha), shape=21, stroke=0, size=3) +
  scale_fill_manual(values = FS_COLORS, name = 'Frameshift') +
  # Hide some legends: https://stackoverflow.com/a/14604540/310453
  scale_alpha_continuous(guide = FALSE) +
  geom_cladelabel(node=findMRCA(my_tree, chl_tips), label="chlD", color='darkgreen', offset.text=.1) +
  geom_hilight(node=findMRCA(my_tree, chl_tips), fill="green", alpha=.3) +
  geom_cladelabel(node=findMRCA(my_tree, bch_tips), label="bchD", color='red', offset.text=.1) +
  geom_hilight(node=findMRCA(my_tree, bch_tips), fill="red", alpha=.3) +
  geom_cladelabel(node=findMRCA(my_tree, cob_chel_tips), label="cobT", color='blue', angle=270, hjust='center', offset.text=.1) +
  geom_hilight(node=findMRCA(my_tree, cob_chel_tips), fill="blue", alpha=.6) +
  ggtitle(title, subtitle = subT) +
  theme_tree(legend.position="bottom")
ggsave('ref_tree_with_fs_chlD.pdf', path = OUT_DIR, width = 6, height = 6)

