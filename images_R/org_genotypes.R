
source("lib.R")

###

# TAR_PHYLA <- c('Proteobacteria', 'Actinobacteria', 'Chloroflexi',  'Euryarchaeota')
TAR_PHYLA <- c('Proteobacteria', 'Actinobacteria', 'Euryarchaeota')

###

all_orgs_df <- read.delim(paste0(DATA_DIR, "orgs_chel.txt"), as.is=TRUE)
head(orgs_df)

# Generate the genotype column
orgs_df <- all_orgs_df %>%
  filter(phylum %in% TAR_PHYLA) %>%
  mutate(num_M_zero = num_M - num_M_fs) %>%
  mutate(L_str = ifelse(num_L == 0, NA, paste0(num_L, 'xL')),
         M_0_str = ifelse(num_M_zero == 0, NA, paste0(num_M_zero, 'xM')),
         M_fs_str = ifelse(num_M_fs == 0, NA, paste0(num_M_fs, 'xfs_chlD')),
         S_str = ifelse(num_S == 0, NA, paste0(num_S, 'xS'))) %>%
  unite(L_str, M_0_str, S_str, M_fs_str, sep = ',', col="genotype") %>%
  # Remove NA from genotypes using regexp:
  # '1xL,1xM_0,NA,NA'   =>  '1xL,1xM_0'
  # '3xL,1xM_0,NA,1xS'  =>  '3xL,1xM_0,1xS'
  mutate(genotype = gsub('(,NA)+', '', genotype))
head(orgs_df)

orgs_df %>%
  group_by(phylum, genotype) %>%
  summarise(n = n()) %>%
  arrange(phylum, -n) %>%
  as.data.frame()

# Expected
EXPECTED_G <- c(
  '1xL,1xM,1xS',    # Expected
  '2xL,2xM,2xS',    # Expected
  '1xL,2xM,2xS',            # Proteobacteria
  
  # with fs_chlD
  '1xL,1xfs_chlD',  # Expected
  '2xL,1xfs_chlD',   # Euryarchaeota (M.f.), Proteobacteria
  '3xL,1xfs_chlD',   # Euryarchaeota (M.s.)
  '2xL,1xM,1xS,1xfs_chlD'  # Proteobacteria
)

NOT_EXPECTED_G <- c(
  # without fs_chlD
  '1xL,1xM',   # Actinobacteria, Euryarchaeota, Proteobacteria
  '1xL,2xM',   # Actinobacteria
  '1xL,1xM,1xS,1xfs_chlD',  # Proteobacteria (D.a.)
  'Other'
)

TOP_GENOTYPES <- c(EXPECTED_G, NOT_EXPECTED_G)

orgs_df %>%
  mutate(genotype = ifelse(genotype %in% TOP_GENOTYPES, genotype, 'Other')) %>%
  # Define the order
  mutate(genotype = factor(genotype, levels = rev(TOP_GENOTYPES))) %>%
  mutate(phylum = factor(phylum, levels = TAR_PHYLA)) %>%
  ggplot() +
  aes(x = genotype, fill = ifelse(num_M_fs > 0, 'Yes', 'No')) +
  geom_bar(col = 'black') +
  ylab('Number of genomes') +
  xlab('Chelatase genotype') +
#  ggtitle(title, subtitle = subT) +
  # scale_fill_manual(
  #   name = 'fs-chlD gene',
  #   values = c('red', 'gray')) +
  scale_fill_manual(
    name = 'Genotype contains fs-chlD gene: ',
    values = c('No' = 'gray', 'Yes' = 'red')) +
  theme(legend.position="top") +
  coord_flip() +
  facet_wrap(.~phylum, nrow = 1, scales = 'free_x')
ggsave(paste0('org_genotypes.pdf'), path = OUT_DIR, w=10, h=5)
