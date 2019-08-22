
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
orgs_df <- all_orgs_df %>%
  mutate(phylum = ifelse(phylum %in% TAR_PHYLA, phylum, 'Other phyla')) %>%
  mutate(num_bchlH = num_chlH + num_bchH,
         num_bchlD = num_chlD + num_bchD,
         num_bchlD_fs = num_chlD_fs + num_bchD_fs, 
         num_bchlD_zero = num_bchlD - num_bchlD_fs,
         num_bchlI = num_chlI + num_bchI) %>%
  mutate(cobN_str = ifelse(num_cobN == 0, NA, paste0(num_cobN, 'xcobN')),
         cobT_str = ifelse(num_cobT == 0, NA, paste0(num_cobT, 'xcobT')),
         cobS_str = ifelse(num_cobS == 0, NA, paste0(num_cobS, 'xcobS')),
         
         bchlH_str = ifelse(num_bchlH == 0, NA, paste0(num_bchlH, 'x(b)chlH')),
         bchlD_0_str = ifelse(num_bchlD_zero == 0, NA, paste0(num_bchlD_zero, 'x(b)chlD')),
         bchlD_fs_str = ifelse(num_bchlD_fs == 0, NA, paste0(num_bchlD_fs, 'xfs-chlD')),
         bchlI_str = ifelse(num_bchlI == 0, NA, paste0(num_bchlI, 'x(b)chlI'))) %>%
  unite(cobN_str, cobT_str, cobS_str, bchlH_str, bchlD_0_str, bchlI_str, bchlD_fs_str,
        sep = ', ', col="genotype") %>%
  # Remove NA from genotypes using regexp:
  mutate(genotype = gsub('(, NA)+', '', genotype)) %>%
  mutate(genotype = gsub('^NA, ', '', genotype)) %>%
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
