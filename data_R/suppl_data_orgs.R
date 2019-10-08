
source("../images_R/lib.R")


###


all_orgs_df <- read.delim(paste0(DATA_DIR, "orgs_chel.txt"), as.is=TRUE)
head(all_orgs_df)

# Generate the genotype column
orgs_df <- make_genotype_column(all_orgs_df)
head(orgs_df)

orgs_df %>%
  select(name, kingdom, phylum, tax2, num_M_fs, genotype) %>%
  rename(num_fs_chlD = num_M_fs) %>%
  arrange(name) %>%
  write.table(file = paste0(DATA_DIR, "suppl_data_orgs.txt"),
              row.names = FALSE, quote = FALSE, sep = "\t")
