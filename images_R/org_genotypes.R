
source("lib.R")

###

orgs_df <- read.delim(paste0(DATA_DIR, "org_genotypes.txt"), as.is=TRUE)
head(orgs_df)

df <- orgs_df %>%
  mutate(num_M_zero = num_M_tree - num_M_fs) %>%
  mutate(L_str = ifelse(num_L == 0, NA, paste0(num_L, 'xL')),
         M_0_str = ifelse(num_M_zero == 0, NA, paste0(num_M_zero, 'xM_0')),
         M_fs_str = ifelse(num_M_fs == 0, NA, paste0(num_M_fs, 'xM_fs')),
         S_str = ifelse(num_S_tree_wo_kids == 0, NA, paste0(num_S_tree_wo_kids, 'xS'))) %>%
  unite(L_str, M_0_str, M_fs_str, S_str, sep = ',', col="genotype") %>%
  # Remove NA from genotypes using regexp:
  # '1xL,1xM_0,NA,NA'   =>  '1xL,1xM_0'
  # '3xL,1xM_0,NA,1xS'  =>  '3xL,1xM_0,1xS'
  mutate(genotype = gsub('(,NA)+', '', genotype))

df %>%
#  filter(num_L == 1) %>%
  group_by(genotype) %>%
  summarise(n_orgs = n()) %>%
  arrange(-n_orgs) %>%
  filter(n_orgs > 10) %>%
  pull(genotype)

TOP_GENOTYPES_ALL <- c("2xL,1xM_0,1xM_fs,1xS", "1xL,1xM_fs", "1xL,1xM_0",
                       "1xL,1xM_0,1xM_fs,1xS", "2xL,2xM_0,1xM_fs,1xS",
                       "1xL,1xM_0,1xS", 'Other')

NUM_L_COLORS <- c(
  '0' = 'white',
  '1' = 'lightgreen',
  '2' = 'orange',
  '>3' = 'red')

PHYLUM = "Proteobacteria"

gg_df <- filter(df, phylum == PHYLUM)

title <- sprintf('Total number of genomes = %d (%s)', nrow(gg_df), PHYLUM)
subT <- sprintf('# chlD genes = %d (with fs) + %d (without fs) = %d',
                sum(gg_df$num_M_fs), sum(gg_df$num_M_zero), sum(gg_df$num_M_tree))

gg_df %>%
  mutate(genotype = ifelse(genotype %in% TOP_GENOTYPES_ALL, genotype, 'Other')) %>%
  # Define the order
  mutate(genotype = factor(genotype, levels = rev(TOP_GENOTYPES_ALL))) %>%
  mutate(num_L_str = ifelse(num_L >= 3, '>3', as.character(num_L))) %>%
  # Define the order
  mutate(num_L_str = factor(num_L_str, levels = c('0', '1', '2', '>3'))) %>%
  ggplot() +
  aes(x = genotype, fill = num_L_str) +
  geom_bar(col = 'black') +
  ylab('Number of genomes') +
  xlab('') +
  ggtitle(title, subtitle = subT) +
  scale_fill_manual(
    name = '# L-genes\nin genotype',
    values = NUM_L_COLORS) +
    # breaks = names(ORG_TYPE_COLORS),
    # labels = ORG_TYPE_LABELS) +
  # # Legend items in one column
  # guides(fill=guide_legend(ncol=1)) +
  # theme(legend.position="bottom") +
  coord_flip()
ggsave('org_genotypes_proteo.pdf', path = OUT_DIR)


# Candidates to search for in-frame frameshifting
# df %>%
#   filter(genotype == '1xL,1xM_0') %>%
#   select(1:6, genotype) %>%
#   arrange(tax2, genus)
