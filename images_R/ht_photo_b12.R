
source("lib.R")

library(ComplexHeatmap)
library(circlize)   # colorRamp2()

library(ape)        # read.tree()
#library(phytools)   # findMRCA()

###

#DATA_DIR <- '~/Projects/2018/Baranov/DataLog/'
# photo_b12_df <- read.table(paste0(DATA_DIR, "/1016.Ba.pathogens/Results/1024.org_photo_b12.txt"),
# tree <- read.tree(paste0(DATA_DIR, "/1010.Ba.rRNA_tree/data/_all_species.rRNA.tree"))

# Heatmap ----
all_data <- dbGetQuery(con,
                       "select id AS org_id, name, phylum, kingdom, pathogen,
                       num_bchlD_minus, num_bchlD_plus,
                       evalue_cobN, evalue_chlH, evalue_bchH,
                       evalue_cobT, evalue_chlD, evalue_bchD,
                       evalue_cobS, evalue_chlI, evalue_bchI
                       from chel_orgs_v")
rownames(all_data) <- all_data$org_id

# Add the 'Photosynthetic' and 'B12' columns
photo_b12_df <- read.table(paste0(DATA_DIR, "org_photo_b12.txt"), header=TRUE, sep = "\t", as.is = TRUE, row.names = 1)
all_data <- merge(all_data, photo_b12_df, by = 'row.names')
rownames(all_data) <- all_data$org_id

# Get orgs sorting from the tree
tree <- read.tree(paste0(DATA_DIR, "all_species_rRNA.tree"))
# "624782022.Actinobacteria.Tomitella"  => "624782022"
tree$tip.label <- gsub("^(\\d+).+", "\\1", tree$tip.label)
sorted_orgs <- tree$tip.label[tree$tip.label %in% rownames(all_data)]
all_data <- all_data[sorted_orgs,]
#tail(all_data)

# evalue => log10_evalue
evalue_cols <- grep('^evalue_', colnames(all_data))
get_log10 <- function(x) ifelse(x < 1e-100, 100, -log10(x))
all_data <- all_data %>%
  mutate_at(evalue_cols, replace_na, 1) %>%
  mutate_at(evalue_cols, funs(get_log10)) %>%
  rename_at(evalue_cols, funs(sub("evalue_", "", .)))

# Make the bchlI, bchlD, bchlH
all_data <- all_data %>%
  mutate(bchlI = pmax(bchI, chlI),
         bchlD = pmax(bchD, chlD),
         bchlH = pmax(bchH, chlH))

# Generate 'taxa' column
taxa_cols <- c("Actinobacteria" = "cyan3", "Archaea" = 'black', "Chloroflexi" = 'blue',
               "Cyanobacteria" = "darkgreen", "Firmicutes" = 'orange', "Other" = "gray", "Proteobacteria" = "plum")
most_frequent <- names(taxa_cols)
all_data <- all_data %>%
  mutate(taxa = ifelse(kingdom == 'Archaea', 'Archaea', phylum)) %>%
  mutate(taxa = ifelse(taxa %in% most_frequent, taxa, 'Other'))
rownames(all_data) <- all_data$org_id

# Generate 'frameshift' column
all_data <- all_data %>%
  mutate(frameshift = ifelse(num_bchlD_minus > 0, '-1', 'None')) %>%
  mutate(frameshift = ifelse(num_bchlD_plus > 0, '+1', frameshift))

str(all_data)

# Statistics
cat(100*round(1 - sum(all_data$Photosynthetic == "Yes") / nrow(all_data), 2), "% of the identified species are not photosynthetic\n",
    100*round(sum(all_data$B12 == "Yes") / nrow(all_data), 2), "% of species may be able to synthesize cobalamin (Vitamin B12)\n",
    sep = '')


# Heatmap
save_ht_pdf <- function(ht_list, base_fn, w = 6, h = 6)
{
  pdf(paste0(OUT_DIR, base_fn, '.pdf'), width = w, height = h)
  draw(ht_list,
       show_heatmap_legend = FALSE,
       row_title = sprintf('%d prokaryotic genomes*', nrow(all_data)))
  dev.off()
}
evalue_colors <- colorRamp2(c(0, 6, 100), c("white", "yellow", "red"))

# taxa_ht
taxa_ht <- Heatmap(data.frame(Taxonomy = all_data$taxa),
                   name = "taxa",
                   col = taxa_cols,
                   width = unit(3, "mm"),
                   cluster_columns = FALSE,
                   cluster_rows = FALSE,
                   show_row_names = FALSE,
                   show_heatmap_legend = FALSE,
                   split = all_data$taxa,
                   row_title_rot = 0,
                   row_title_gp = gpar(col = taxa_cols, font = 2),
                   gap = unit(3, "mm"))
#?Heatmap
#quartz()
#taxa_ht

# fs_ha
fs_ha <- rowAnnotation(data.frame(Frameshift = all_data$frameshift),
                       col = list(Frameshift = c('-1' = 'black', '+1' = 'red', 'None' = 'white')),
                       annotation_legend_param = list(Frameshift = list(title = "Frameshift\nin chlD gene")),
                       show_annotation_name = TRUE,
                       width = unit(1, "cm"))
# taxa_ht + fs_ha

# chlD_ht
chlD_ht <- Heatmap(data.frame('chlD' = all_data$bchlD),
                   #all_data[, 'bchlD', drop = FALSE],
                   col = evalue_colors,
                   cluster_columns = FALSE,
                   heatmap_legend_param = list(title = "BLAST\n-log10(E-value)"),
                   width = unit(1, "cm"))
# taxa_ht + fs_ha + chlD_ht

# patho_ha: http://www.bioconductor.org/packages/release/bioc/vignettes/ComplexHeatmap/inst/doc/s4.heatmap_annotation.html#toc_14
patho_rows <- which(all_data$pathogen == 1)
patho_names <- all_data[patho_rows, 'name']
patho_ha <- rowAnnotation(link = row_anno_link(at = patho_rows, labels = patho_names),
                          width = unit(1, "cm") + max_text_width(patho_names))
# save_ht_pdf(taxa_ht + chlD_ht + fs_ha + patho_ha,
#             base_fn = 'ht_1.chlD')

# photo_ha
photo_ha <- rowAnnotation(data.frame(Chlorophyll = all_data$Photosynthetic),
                          col = list(Chlorophyll = c('Yes' = 'green', 'No' = 'white')),
                          width = unit(1, "cm"),
                          show_annotation_name = TRUE,
                          show_legend = FALSE)

# chlIDH_ht
chlIDH_ht <- Heatmap(data.frame('chlH' = all_data$bchlH,
                                'chlD' = all_data$bchlD,
                                'chlI' = all_data$bchlI),
                     col = evalue_colors,
                     cluster_columns = FALSE,
                     width = unit(3, "cm"))
# save_ht_pdf(taxa_ht + photo_ha + chlIDH_ht + fs_ha + patho_ha, w = 8,
#             base_fn = 'ht_2.chlD_photo')

# b12_ha
b12_ha <- rowAnnotation(data.frame(Cobalamin = all_data$B12),
                        col = list(Cobalamin = c('Yes' = 'blue', 'No' = 'white')),
                        width = unit(1, "cm"),
                        show_annotation_name = TRUE,
                        show_legend = FALSE)

# cobNST_ht
cobNST_ht <- Heatmap(all_data[, c('cobN', 'cobT', 'cobS')],
                     col = evalue_colors,
                     cluster_columns = FALSE,
                     width = unit(3, "cm"))
save_ht_pdf(taxa_ht + photo_ha + chlIDH_ht + b12_ha + cobNST_ht + fs_ha + patho_ha, w = 10,
            base_fn = 'ht_photo_b12')

