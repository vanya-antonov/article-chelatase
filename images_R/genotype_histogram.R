
source("lib.R")

###

# TAR_PHYLA <- c('Proteobacteria', 'Actinobacteria', 'Chloroflexi',  'Euryarchaeota')
# TAR_PHYLA <- c('Proteobacteria', 'Actinobacteria', 'Euryarchaeota')
TAR_PHYLA <- c('Proteobacteria', 'Actinobacteria', 'Euryarchaeota', 'Cyanobacteria', 'Other phyla')

EXPECTED_GENOTYPES <- c(
  "1xcobN, 1xcobT, 1xcobS, 1x(b)chlH, 1x(b)chlD, 1x(b)chlI",
  "1xcobN, 1xfs-chlD",
  "1x(b)chlH, 1x(b)chlD, 1x(b)chlI",
  "1xcobN, 1x(b)chlD, 1x(b)chlI")

###

all_orgs_df <- read.delim(paste0(DATA_DIR, "orgs_chel.txt"), as.is=TRUE)
head(all_orgs_df)

# Generate the genotype column
orgs_df <- make_genotype_column(all_orgs_df) %>%
  mutate(phylum = ifelse(phylum %in% TAR_PHYLA, phylum, 'Other phyla')) %>%
  # Add * to the expected genotypes
  mutate(genotype = ifelse(genotype %in% EXPECTED_GENOTYPES, paste0(genotype, '*'), genotype))
head(orgs_df)
  
TOP_GENOTYPES <- orgs_df %>%
  group_by(genotype) %>%
  summarise(n = n()) %>%
  arrange(-n) %>%
  pull(genotype) %>%
  head(n = 8)
TOP_GENOTYPES <- c(TOP_GENOTYPES, 'Other')
  
orgs_df %>%
  mutate(genotype = ifelse(genotype %in% TOP_GENOTYPES, genotype, 'Other')) %>%
  # Define the order
  mutate(genotype = factor(genotype, levels = rev(TOP_GENOTYPES))) %>%
  mutate(phylum = factor(phylum, levels = TAR_PHYLA)) %>%
  ggplot() +
  aes(x = genotype, fill = phylum) +
  geom_bar(col = 'black') +
  ylab('Number of genomes') +
  xlab('Chelatase genotype') +
#  scale_fill_manual(
#    name = 'Genotype contains fs-chlD gene: ',
#    values = c('No' = 'gray', 'Yes' = 'red')) +
#  theme(legend.position="top") +
  coord_flip()
ggsave(paste0('genotype_histogram.pdf'), path = OUT_DIR, w=13, h=6)


# Statistics ----
orgs_df %>%
  group_by(genotype) %>%
  summarise(n = n()) %>%
  arrange(n) %>%
  as.data.frame()
  

# Number of Proteobacteria with expected genotypes: 175 (54.5%)
n_proteo <- orgs_df %>% filter(phylum == 'Proteobacteria') %>% nrow()
n_good_proteo <- orgs_df %>% filter(phylum == 'Proteobacteria', genotype %in% EXPECTED_G) %>% nrow()
sprintf("%d (%.1f%%)", n_good_proteo, 100*n_good_proteo/n_proteo)

# Expected genotypes with fs-chlD: 59
orgs_df %>% filter(phylum == 'Proteobacteria', genotype %in% EXPECTED_G, num_M_fs > 0) %>% nrow()

# Expected genotypes w/o fs-chlD: 116
orgs_df %>% filter(phylum == 'Proteobacteria', genotype %in% EXPECTED_G, num_M_fs == 0) %>% nrow()

