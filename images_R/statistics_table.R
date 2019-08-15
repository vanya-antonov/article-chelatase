
source("lib.R")

###

orgs_df <- read.delim(paste0(DATA_DIR, "orgs_all.txt"), as.is=TRUE)
head(orgs_df, n = 4)

# Bacteria phyla with fs-chlD
orgs_df %>%
  filter(kingdom == 'Bacteria', num_bchlD_fs > 0) %>%
  pull(phylum) %>%
  unique()
tar_phyla <- c('Euryarchaeota', 'Crenarchaeota', 'Thaumarchaeota',
               "Proteobacteria", "Actinobacteria", "Chloroflexi", "Spirochaetes", "Firmicutes", "Bacteroidetes", "Cyanobacteria",
               'Other')

orgs_df %>%
  mutate(phylum2 = ifelse(phylum %in% tar_phyla, phylum, 'Other')) %>%
  # Define the order of phyla
  mutate(phylum2 = factor(phylum2, levels = tar_phyla)) %>%
  group_by(kingdom, phylum2) %>%
  summarise(total_orgs = n(),
            orgs_with_M = sum(ifelse(num_M > 0, 1, 0)),
            num_M_genes = sum(num_M),
            orgs_with_bchlD_fs = sum(num_bchlD_fs > 0),
            num_bchlD_fs = sum(num_bchlD_fs)) %>%
  write.table(file = paste0(OUT_DIR, "statistics_table.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
