
source('lib.R')

###

CHEL_GENE_NAMES <- c(
  'chlI', 'chlD', 'chlH',
  'bchI', 'bchD', 'bchH',
  'cobS', 'cobT', 'cobN')

###

all_chel_genes <- read.delim(paste0(DATA_DIR, "genes_chel.txt"), as.is=TRUE) %>%
  filter(gene %in% CHEL_GENE_NAMES)

title <- sprintf('All annotated genes = %s', nrow(all_chel_genes))
all_chel_genes %>%
  ggplot() +
  aes(x = prot_len, fill = gene) +
  geom_density(alpha = 0.4) +
  scale_x_continuous(breaks = seq(0, 1500, 100), limits = c(200, 1500)) +
  ggtitle(title) +
  theme_bw() +
  xlab('Protein length (aa)') +
  facet_grid(gene ~ ., scale = 'free_y')
ggsave('known_chel_lengths.pdf', path=OUT_DIR)
