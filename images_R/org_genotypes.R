
source("lib.R")

###

NUM_L_COLORS <- c(
  '0' = 'white',
  '1' = 'lightgreen',
  '2' = 'orange',
  '>3' = 'red')

###

orgs_df <- read.delim(paste0(DATA_DIR, "org_genotypes.txt"), as.is=TRUE)
head(orgs_df)

# Generate the genotype column
orgs_df <- orgs_df %>%
  mutate(num_M_zero = num_M - num_M_fs) %>%
  mutate(L_str = ifelse(num_L == 0, NA, paste0(num_L, 'xL')),
         M_0_str = ifelse(num_M_zero == 0, NA, paste0(num_M_zero, 'xM_0')),
         M_fs_str = ifelse(num_M_fs == 0, NA, paste0(num_M_fs, 'xM_fs')),
         S_str = ifelse(num_S == 0, NA, paste0(num_S, 'xS'))) %>%
  unite(L_str, M_0_str, M_fs_str, S_str, sep = ',', col="genotype") %>%
  # Remove NA from genotypes using regexp:
  # '1xL,1xM_0,NA,NA'   =>  '1xL,1xM_0'
  # '3xL,1xM_0,NA,1xS'  =>  '3xL,1xM_0,1xS'
  mutate(genotype = gsub('(,NA)+', '', genotype))


genotype_histogram_gg <- function(gg_df, top_genotypes)
{
  top_genotypes <- c(top_genotypes, 'Other')
  title <- sprintf('Total number of genomes = %d', nrow(gg_df))
  subT <- sprintf('# chlD genes = %d (with fs) + %d (without fs) = %d',
                  sum(gg_df$num_M_fs), sum(gg_df$num_M_zero), sum(gg_df$num_M))
  gg_df %>%
    mutate(genotype = ifelse(genotype %in% top_genotypes, genotype, 'Other')) %>%
    # Define the order
    mutate(genotype = factor(genotype, levels = rev(top_genotypes))) %>%
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
    coord_flip()
}

PHYLUM = "Proteobacteria"
TOP_GENOTYPES_PROTEO <- c("2xL,1xM_0,1xM_fs,1xS", "1xL,1xM_fs", "1xL,1xM_0",
                          "1xL,1xM_0,1xM_fs,1xS", "2xL,2xM_0,1xM_fs,1xS","1xL,1xM_0,1xS")

orgs_df %>%
  filter(phylum == PHYLUM) %>%
  genotype_histogram_gg(TOP_GENOTYPES_PROTEO)
ggsave(paste0('org_genotypes_', PHYLUM, '.pdf'), path = OUT_DIR)


# Most frequent genotype
# orgs_df %>%
#   filter(phylum == PHYLUM) %>%
#   group_by(genotype) %>%
#   summarise(n_orgs = n()) %>%
#   arrange(-n_orgs) %>%
#   head(n = 6) %>%
#   pull(genotype)
