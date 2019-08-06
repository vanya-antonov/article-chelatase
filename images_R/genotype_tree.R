
source("lib.R")

###

# https://bioc.ism.ac.jp/packages/3.3/bioc/vignettes/ggtree/inst/doc/advanceTreeAnnotation.html
# the new ggtree-way: https://guangchuangyu.github.io/2016/09/ggtree-for-microbiome-data/

# TO BE UPDATED!!!
tree <- read.tree(paste0(DATA_DIR, "orgs_chel.tree"))
length(tree$tip.label)

orgs_df <- read.delim(paste0(DATA_DIR, "org_genotypes.txt"), as.is=TRUE)
rownames(orgs_df) <- orgs_df$dir_name
head(orgs_df)

# Sync the data -- some genomes do not have annotated rRNAs
common_orgs <- intersect(rownames(orgs_df), tree$tip.label)
orgs_df <- orgs_df[common_orgs,]


tar_phylum <- 'Proteobacteria'
tar_orgs_df <- orgs_df %>% filter(phylum == tar_phylum, num_M == 1, num_L > 0)
tar_tree <- keep.tip(tree, tar_orgs_df$dir_name)

small_genotype <- tar_orgs_df %>%
  mutate('L' = ifelse(num_L > 0, 'None', ''),
         'M' = ifelse(num_M_plus > 0, '+1', ifelse(num_M_minus > 0, '-1', 'None')),
         'S' = ifelse(num_S > 0, 'None', '')
  ) %>%
  column_to_rownames('dir_name') %>%
  select('L', 'M', 'S')

Pseudomonas_node <- tar_orgs_df %>%
  filter(genus %in% c('Pseudomonas')) %>%
  pull(dir_name) %>%
  findMRCA(tar_tree, .)

Burkholderia_node <- tar_orgs_df %>%
  filter(genus %in% c('Burkholderia', 'Paraburkholderia')) %>%
  pull(dir_name) %>%
  findMRCA(tar_tree, .)

title <- paste(nrow(tar_orgs_df), tar_phylum, 'species with exactly 1 chlD gene')

ggtree(tar_tree) %>%
  gheatmap(small_genotype, offset=0.3, width=0.5, colnames_offset_y = -1.5) +
  #    geom_tiplab(size=2, align=TRUE, linesize=.5) +
  scale_fill_manual(values=c("None"="black", "-1"="blue", "+1"="red")) +
  geom_cladelabel(node=Pseudomonas_node, label="Pseudomonas") +
  geom_cladelabel(node=Burkholderia_node, label='Burkholderia &\nParaburkholderia') +
  ggtitle(title)
ggsave('genotype_tree_Proteobacteria.pdf', path = OUT_DIR,
       width = 7, height = 5)


# Statistics
# Frameshift in chlD   |   chlI gene is absent
# Y | Y = 65
YY <- tar_orgs_df %>% filter(num_M_fs > 0 & num_S == 0) %>% nrow()
# Y | N = 0
YN <- tar_orgs_df %>% filter(num_M_fs > 0 & num_S > 0) %>% nrow()
# N | Y = 38
NY <- tar_orgs_df %>% filter(num_M_fs == 0 & num_S == 0) %>% nrow()
# N | N = 14
NN <- tar_orgs_df %>% filter(num_M_fs == 0 & num_S > 0) %>% nrow()


mtx <- matrix(c(YY, YN,
                NY, NN),
              ncol=2, byrow=TRUE)
fisher.test(mtx)
