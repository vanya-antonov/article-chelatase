
# install.packages('cowplot')
# install.packages('RMySQL')

library(ggplot2)
library(cowplot)   # To have ggplots side-by-side: plot_grid()

library(dplyr)
library(tidyr)     # separate(), gather() and spread() functions
library(tibble)    # for rownames_to_column() and column_to_rownames()

library(RMySQL)

library(circlize)   # colorRamp2()

# library(ComplexHeatmap)
# library(circlize)   # colorRamp2()
# 
library(ape)        # read.tree()
library(phytools)   # findMRCA()

# https://bioc.ism.ac.jp/packages/3.3/bioc/vignettes/ggtree/inst/doc/treeAnnotation.html
library(ggtree)


###

DATA_DIR <- "../data/"
OUT_DIR <- "../images/"

FS_COLORS <- c('Yes' = 'black', 'No' = 'white')
EVALUE_COLORS <- colorRamp2(c(0, 6, 100), c("white", "gray", "black"))

theme_set(theme_bw(base_size = 19))  # increase the font size: https://stackoverflow.com/a/11955412/310453


###

make_genotype_column <- function(all_orgs_df)
{
  # all_orgs_df  is the content of the 'orgs_chel.txt'
  all_orgs_df %>%
    mutate(num_bchlH = num_chlH + num_bchH,
           num_bchlD = num_chlD + num_bchD,
           num_bchlD_fs = num_chlD_fs + num_bchD_fs, 
           num_bchlD_zero = num_bchlD - num_bchlD_fs,
           num_bchlI = num_chlI + num_bchI) %>%
    mutate(cobN_str = ifelse(num_cobN == 0, NA, paste0(num_cobN, 'xcobN')),
           cobT_str = ifelse(num_cobT == 0, NA, paste0(num_cobT, 'xcobT')),
           cobS_str = ifelse(num_cobS == 0, NA, paste0(num_cobS, 'xcobS')),
           
           bchlH_str = ifelse(num_bchlH == 0, NA, paste0(num_bchlH, 'xchlH')),
           bchlD_0_str = ifelse(num_bchlD_zero == 0, NA, paste0(num_bchlD_zero, 'xchlD')),
           bchlD_fs_str = ifelse(num_bchlD_fs == 0, NA, paste0(num_bchlD_fs, 'xfs-chlD')),
           bchlI_str = ifelse(num_bchlI == 0, NA, paste0(num_bchlI, 'xchlI'))) %>%
    unite(cobN_str, cobT_str, cobS_str, bchlH_str, bchlD_0_str, bchlI_str, bchlD_fs_str,
          sep = ', ', col="genotype") %>%
    # Remove NA from genotypes using regexp:
    mutate(genotype = gsub('1x', '', genotype)) %>%
    # Remove NA from genotypes using regexp:
    mutate(genotype = gsub('(, NA)+', '', genotype)) %>%
    mutate(genotype = gsub('^NA, ', '', genotype))
}

read_tree_from_fn <- function(fn)
{
  # Genes that do not have good position in the REF_TREE
  BAD_GENES <- c('AMED_7054', 'CAGG_RS15520', 'RCAS_RS13095', 'DSHI_RS17730', 'RSP_0274', 'RSPH17025_RS05185', 'BLV13_RS00120', 'SADFL11_RS19460')
  tree <- read.tree(paste0(DATA_DIR, fn))
  tree <- drop.tip(tree, BAD_GENES)
  return(tree)
}

read_anno_for_tree <- function(anno_fn, tree)
{
  anno_df <- read.delim(paste0(DATA_DIR, anno_fn), as.is=TRUE) %>%
    filter(name %in% tree$tip.label) %>%
    mutate(fs_len = ifelse(is.na(fs_len), 'None', sprintf('%+d', fs_len))) %>%
    mutate(fs_len = factor(fs_len, levels = names(FS_COLORS))) %>%
    # There shouldn't be any shapes in genes w/o fs
    mutate(shape_alpha = ifelse(fs_len == 'None', 0, 1))
  rownames(anno_df) <- anno_df$name
  return(anno_df)
}
