
source("lib.R")

###

# TAR_PHYLA <- c('Proteobacteria', 'Actinobacteria', 'Chloroflexi',  'Euryarchaeota')
# TAR_PHYLA <- c('Proteobacteria', 'Actinobacteria', 'Euryarchaeota')
TAR_PHYLA <- c('Proteobacteria', 'Actinobacteria', 'Euryarchaeota', 'Cyanobacteria', 'Other phyla')

###

all_orgs_df <- read.delim(paste0(DATA_DIR, "orgs_chel.txt"), as.is=TRUE)
head(all_orgs_df)

# Generate the genotype column
orgs_df <- make_genotype_column(all_orgs_df) %>%
  mutate(phylum = ifelse(phylum %in% TAR_PHYLA, phylum, 'Other phyla'))
head(orgs_df)
  
TOP_GENOTYPES <- orgs_df %>%
  group_by(genotype) %>%
  summarise(n = n()) %>%
  arrange(-n) %>%
  pull(genotype) %>%
  head(n = 9)
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

