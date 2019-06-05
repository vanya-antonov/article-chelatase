
source("lib.R")

###

# Histogram by phylum ----

# df <- dbGetQuery(con, "select phylum, kingdom, num_bchlD_zero, num_bchlD_minus, num_bchlD_plus from chel_orgs_v")
df <- dbGetQuery(con, "select phylum, kingdom, num_M_zero, num_M_minus, num_M_plus from chel_orgs_v")

most_frequent <- c("Actinobacteria", "Proteobacteria", "Archaea", "Cyanobacteria", "Firmicutes", "Chloroflexi", 'Other')
df <- df %>%
  mutate(phylum = ifelse(kingdom == 'Archaea', 'Archaea', phylum)) %>%
  mutate(phylum = ifelse(phylum %in% most_frequent, phylum, 'Other'))

gg_df <- df %>%
  group_by(phylum) %>%
  summarise('None' = sum(num_M_zero),
            '+1' = sum(num_M_plus),
            '-1' = sum(num_M_minus)) %>%
  # summarise('None' = sum(num_bchlD_zero),
  #           '+1' = sum(num_bchlD_plus),
  #           '-1' = sum(num_bchlD_minus)) %>%
  gather(Frameshift, num, -phylum)

#RColorBrewer::display.brewer.all()
ggplot(gg_df) + 
  aes(x = factor(phylum, levels = rev(most_frequent)), y = num, fill = Frameshift) +
  geom_bar(stat = 'identity', col = 'black') + 
  ylab('Number of medium subunit genes') +
  xlab('Phylum') +
  #  scale_fill_brewer(palette = "Set1") +
  coord_flip()
ggsave('all_orgs_hist.pdf', path = OUT_DIR, width = 7)
