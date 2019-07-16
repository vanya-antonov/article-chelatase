
source("lib.R")

###
# Histogram by phylum ----

orgs_df <- read.delim(paste0(DATA_DIR, "orgs_all.txt"), as.is=TRUE) %>%
  mutate(org_type = case_when(num_bchlD_fs > 0 ~ 'bchlD_fs',
                              num_bchlD > 0 ~ 'bchlD_wo_fs',
                              num_bchlD == 0 ~ 'no_bchlD'))
head(orgs_df, n = 20)

most_frequent <- c("Proteobacteria", "Actinobacteria", "Archaea", "Chloroflexi", "Firmicutes", "Cyanobacteria", 'Other')
#most_frequent <- c("Proteobacteria", "Actinobacteria", "Archaea", "Firmicutes", "Cyanobacteria", 'Other')

ORG_TYPE_COLORS <- c(
  'bchlD_fs' = 'red',
  'bchlD_wo_fs' = 'lightgreen',
  'no_bchlD' = 'grey')
ORG_TYPE_LABELS <- c(
  'bchlD_fs' = 'Frameshifted chlD gene(s)',
  'bchlD_wo_fs' = 'Normal chlD gene(s)',
  'no_bchlD' = 'No chlD gene(s)')

orgs_df %>%
  mutate(taxa = case_when(kingdom == 'Archaea' ~ 'Archaea',
                          phylum %in% most_frequent ~ phylum,
                          TRUE ~ 'Other')) %>%
  # Define the order of phyla
  mutate(taxa = factor(taxa, levels = rev(most_frequent))) %>%
  # Define the order of org_types
  mutate(org_type = factor(org_type, levels = rev(names(ORG_TYPE_COLORS)))) %>%
  ggplot() +
  aes(x = taxa, fill = org_type) +
  geom_bar(col = 'black') +
  ylab('Number of species') +
  xlab('') +
  ggtitle(paste(nrow(orgs_df), 'prokaryotic species'),
          subtitle = 'NCBI Reference and Representative genomes') +
  scale_fill_manual(
    name = NULL,
    values = ORG_TYPE_COLORS,
    breaks = names(ORG_TYPE_COLORS),
    labels = ORG_TYPE_LABELS) +
  # Legend items in one column
#  guides(fill=guide_legend(ncol=1)) +
  coord_flip() +
  theme(legend.position="bottom")
ggsave('all_orgs_hist.pdf', path = OUT_DIR)
#ggsave('all_orgs_hist.pdf', path = OUT_DIR, width = 7)

#?scale_fill_manual
