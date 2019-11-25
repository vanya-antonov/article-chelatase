
source("lib.R")

theme_set(theme_bw(base_size = 25))  # increase the font size: https://stackoverflow.com/a/11955412/310453


###

# TAR_PHYLA <- c('Proteobacteria', 'Actinobacteria', 'Chloroflexi',  'Euryarchaeota')
# TAR_PHYLA <- c('Proteobacteria', 'Actinobacteria', 'Euryarchaeota')
TAR_PHYLA <- c('Proteobacteria', 'Actinobacteria', 'Euryarchaeota', 'Cyanobacteria', 'Other phyla')

###

all_orgs_df <- read.delim(paste0(DATA_DIR, "orgs_chel.txt"), as.is=TRUE)
head(all_orgs_df)

# Generate the genotype column
orgs_df <- make_genotype_column(all_orgs_df) %>%
  filter(num_M_fs > 0) %>%
  mutate(phylum = ifelse(phylum %in% TAR_PHYLA, phylum, 'Other phyla')) %>%
  mutate(genotype = ifelse(num_L > 0, paste0(genotype, '*'), genotype))
head(orgs_df)

TOP_GENOTYPES <- orgs_df %>%
  group_by(genotype) %>%
  summarise(n = n()) %>%
  arrange(-n) %>%
  pull(genotype)

title <- sprintf('Total number of genomes = %s', nrow(orgs_df))
subT <- sprintf('Number of fs-chlD genes = %s', sum(orgs_df$num_M_fs))
orgs_df %>%
#  mutate(genotype = ifelse(genotype %in% TOP_GENOTYPES, genotype, 'Other')) %>%
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
  coord_flip() +
  ggtitle(title, subtitle = subT) +
  theme(legend.position="bottom") +
  guides(fill = guide_legend(ncol = 1))
ggsave(paste0('genotype_fs_chlD_histogram.pdf'), path = OUT_DIR, w=13, h=15)

