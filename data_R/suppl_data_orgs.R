
source("../images_R/lib.R")


###

all_orgs_df <- read.delim(paste0(DATA_DIR, "orgs_chel.txt"), as.is=TRUE) %>%
  # Generate the genotype column
  make_genotype_column()
head(all_orgs_df)

all_genes_df <- read.delim(paste0(DATA_DIR, "genes_chel.txt"), as.is=TRUE) %>%
  # Use 'chlIDH' for 'bchIDH'
  mutate(chel_gene = case_when(chel_gene == 'bchI' ~ 'chlI',
                               chel_gene == 'bchD' ~ 'chlD',
                               chel_gene == 'bchH' ~ 'chlH',
                               TRUE ~ chel_gene))
unique(all_genes_df$chel_gene)
head(all_genes_df)

# Concatenate gene IDs and merge with orgs
orgs_df <- all_genes_df %>%
  group_by(org_name, chel_gene) %>%
  summarise(gene_id = paste(name, collapse=", ")) %>%
  spread(chel_gene, gene_id, fill = '-') %>%
  inner_join(all_orgs_df, by= c('org_name' = 'name'))

orgs_df %>%
  select(org_name, kingdom, phylum, tax2, num_M_fs, genotype, cobN, cobT, cobS, chlH, chlD, chlI) %>%
  rename(num_fs_chlD = num_M_fs) %>%
  arrange(org_name) %>%
  write.table(file = paste0(DATA_DIR, "suppl_data_orgs.txt"),
              row.names = FALSE, quote = FALSE, sep = "\t")
