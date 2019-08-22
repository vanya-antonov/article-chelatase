
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
orgs_df <- all_orgs_df %>%
  filter(num_M_fs > 0) %>%
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


# Statistics ----

# "the small chelatase subunit gene was absent in 140 out of the 150 fs-chlD containing genomes
# that had at least one large chelatase subunit gene"
n_orgs_with_L <- orgs_df %>% filter(num_L > 0) %>% nrow()
n_orgs_with_L_wo_S <- orgs_df %>% filter(num_L > 0 & num_S == 0) %>% nrow()
sprintf("the small chelatase subunit gene was absent in %s out of the %s fs-chlD containing genomes (%.0f%%)...",
        n_orgs_with_L_wo_S, n_orgs_with_L, 100*n_orgs_with_L_wo_S/n_orgs_with_L)


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

