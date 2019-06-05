
source("lib.R")

library(ape)       # read.tree()
#packageVersion("ape")
library(ggtree)

###

# 2016-phyloseq way: https://joey711.github.io/phyloseq/plot_tree-examples.html
# the new ggtree-way: https://guangchuangyu.github.io/2016/09/ggtree-for-microbiome-data/
# https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/treeAnnotation.html

all_data <- dbGetQuery(con,
                       "
                       SELECT DISTINCT genus, tax2, phylum, kingdom, genotype_LDI AS genotype,
                       num_L, num_cobN, num_bchlD_zero, num_bchlD_minus, num_bchlD_plus, num_chlI, num_bchI,
                       (select SUM(pathogen) from chel_orgs_v t where t.genus=o.genus) AS pathogen,
                       (select COUNT(*) from chel_orgs_v t where t.genus=o.genus and t.genotype_LDI = o.genotype_LDI) AS num_orgs,
                       (select COUNT(*) from chel_orgs_v t where t.genus=o.genus) AS total_orgs
                       FROM chel_orgs_v o;
                       ")
all_data <- all_data %>%
  # Select the most frequent genotype for each genus
  arrange(genus, -num_orgs) %>%
  distinct(genus, .keep_all = TRUE) %>%
  # Add new columns
  mutate(num_bchlI = num_bchI + num_chlI,
         percent_orgs = num_orgs/total_orgs,
         descr = sprintf("%s (%.0f%% of %d genome(s))", genus, 100*percent_orgs, total_orgs),
#         pathogen_str = ifelse(pathogen > 0, descr, ''),
         pathogen_str = ifelse(pathogen > 0, genus, '')) %>%
  mutate(TMP = genus) %>% column_to_rownames('TMP')

# Filter by common genera
#tree <- read.tree("../1010.Ba.rRNA_tree/data/_all_genera.rRNA.tree")
tree <- read.tree(paste0(DATA_DIR, "all_genera_rRNA.tree"))
# "624782022.Actinobacteria.Tomitella"  => "Tomitella"
tree$tip.label <- gsub("^.+\\.(\\w+)$", "\\1", tree$tip.label)
common_genera <- intersect(tree$tip.label, all_data$genus)
all_data <- all_data[common_genera,]
tree <- keep.tip(tree, common_genera)

str(all_data)

#tar_phylum <- 'Proteobacteria'
get_ggtree_for_phylum <- function(tar_phylum)
{
  small_data <- all_data %>% filter(phylum == tar_phylum)
  small_tree <- keep.tip(tree, small_data$genus)
  small_genotype <- small_data %>%
    mutate('chlD' = ifelse(num_bchlD_plus > 0, '+1', ifelse(num_bchlD_minus > 0, '-1', 'None')),
           'L' = ifelse(num_L > 0, 'None', ''),
           'chlI' = ifelse(num_bchlI > 0, 'None', '')
    ) %>%
    select(genus, 'chlD', 'L', 'chlI') %>%
    column_to_rownames('genus')
  
  # See the "Modify (tip) labels" sections from: https://guangchuangyu.github.io/software/ggtree/faq/
  labels_df <- small_data %>% select(label=genus, label2=pathogen_str)
  
  # ?geom_tiplab
  ggt <- ggtree(small_tree) %<+% labels_df +
    geom_tiplab(aes(label=label2), size=4, align=TRUE, linesize=.5) +
    theme_tree2()
  ggt <- ggt %>%
    gheatmap(small_genotype, offset=0.4, width=0.4, colnames = FALSE) +
    # ggtitle(paste(nrow(small_data), tar_phylum, 'genera')) +
    scale_fill_manual(values=c("None"="black", "-1"="blue", "+1"="red"))
  ggt %>%
    scale_x_ggtree()  # Add labels to X-axis
}

proteo_gg <- get_ggtree_for_phylum('Proteobacteria')
#ggsave('Proteobacteria.pdf', path = OUT_DIR, bg = "transparent")

actino_gg <- get_ggtree_for_phylum('Actinobacteria')
#ggsave('Actinobacteria.pdf', path = OUT_DIR)

plot_grid(proteo_gg, actino_gg, nrow = 1, labels = 'AUTO')
ggsave('genera_trees.pdf', path=OUT_DIR, width = 10, height = 7)
