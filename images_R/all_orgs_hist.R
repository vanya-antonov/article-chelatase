
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
  'bchlD_fs' = 'Genomes with frameshifted chlD gene(s)',
  'bchlD_wo_fs' = 'Genomes with normal medium subunit gene(s)',
  'no_bchlD' = 'Genomes without medium subunit gene(s)')

num_fs_genes <- sum(orgs_df$num_bchlD_fs)
num_fs_orgs <- orgs_df %>% filter(num_bchlD_fs > 0) %>% nrow
subT <- sprintf('Number of chlD genes with frameshifts = %d (from %d genomes)',
                num_fs_genes, num_fs_orgs)
title <- paste('Total number of prokaryotic genomes = ', nrow(orgs_df))

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
  ylab('Number of genomes') +
  xlab('') +
  ggtitle(title, subtitle = subT) +
  scale_fill_manual(
    name = NULL,
    values = ORG_TYPE_COLORS,
    breaks = names(ORG_TYPE_COLORS),
    labels = ORG_TYPE_LABELS) +
  # Legend items in one column
  guides(fill=guide_legend(ncol=1)) +
  coord_flip() +
  theme(legend.position="bottom")
ggsave('all_orgs_hist.pdf', path = OUT_DIR)
