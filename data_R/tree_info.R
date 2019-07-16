
# Get the con object (DBI::dbConnect)
source("private.R")

# Load the DATA_DIR
source("../images_R/lib.R")

###

all_genes_df <- dbGetQuery(con, "select id, name, descr, gene, fs_len, child_fs_len, org_name, kingdom, phylum, genus from chel_feats_v")
rownames(all_genes_df) <- as.character(all_genes_df$id)
head(all_genes_df)

filter_and_save <- function(all_genes_df, prefix)
{
  #prefix <- 'all_proteobacteria_M'
  tree <- read.tree(paste0(DATA_DIR, prefix, '.tree'))
  
  filtered_genes_df <- all_genes_df[tree$tip.label,]
  write.table(filtered_genes_df, file = paste0(DATA_DIR, prefix, ".info.txt"),
              row.names = FALSE, quote = FALSE, sep = "\t")
  
  # Make sure we have info for all feats
  all(tree$tip.label %in% rownames(filtered_genes_df))
}

filter_and_save(all_genes_df, 'all_proteobacteria_M')
filter_and_save(all_genes_df, 'all_proteobacteria_S')
