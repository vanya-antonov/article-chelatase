
source('lib.R')

library(ComplexHeatmap)
library(circlize)   # colorRamp2()

###

MG_CHEL_NAMES  <- c('chlH_bchH', 'chlD_bchD', 'chlI_bchI')
COB_CHEL_NAMES <- c('cobN', 'cobT', 'cobS')

CHL_PATH_GENES <- c("bchE", "chlB_bchB",  "chlG_bchG", "chlL_bchL", "chlM_bchM", "chlN_bchN")
B12_PATH_GENES <- c("cobD_cobC", "cobO", "cobP_cobU", "cobQ", "cobV_cobS", "cysG_cobA")

EVALUE_COLORS <- colorRamp2(c(0, 6, 100), c("white", "yellow", "red"))

TAXA_COLS <- c("Proteobacteria" = "plum", "Actinobacteria" = "cyan3", "Archaea" = 'black', "Chloroflexi" = 'blue',
               "Firmicutes" = 'orange', "Cyanobacteria" = "darkgreen", "Other" = "gray")

HT_COL_W = 0.8

# Global heatmap parameters
# https://jokergoo.github.io/ComplexHeatmap-reference/book/a-list-of-heatmaps.html#change-parameters-globally
ht_opt(
#  heatmap_column_title_gp = gpar(fontsize = 10),
  legend_border = "black",
  heatmap_border = TRUE)

###

all_data <- read.delim(paste0(DATA_DIR, "orgs_chel.txt"), header=TRUE, as.is = TRUE)
rownames(all_data) <- all_data$dir_name
tail(all_data)

# Add evalues
all_evalue <- read.delim(paste0(DATA_DIR, "orgs_evalue.txt"), header=TRUE, as.is = TRUE, row.names=1)
tail(all_evalue)

all_data <- merge(all_data, all_evalue, by = 'row.names') %>%
  column_to_rownames('Row.names')
tail(all_data)

# evalue => log10_evalue
evalue_cols <- grep('^evalue_', colnames(all_data))
get_log10 <- function(x) ifelse(x < 1e-100, 100, -log10(x))
all_data <- all_data %>%
  mutate_at(evalue_cols, replace_na, 1) %>%
  mutate_at(evalue_cols, get_log10) %>%
  rename_at(evalue_cols, funs(sub("evalue_", "", .)))
rownames(all_data) <- all_data$dir_name
head(all_data)

# Add the 'Photosynthetic' and 'B12' columns
# photo_b12_df <- read.table(paste0(DATA_DIR, "181024.org_photo_b12.txt"), header=TRUE, sep = "\t", as.is = TRUE, row.names = 1)
# all_data <- merge(all_data, photo_b12_df, by = 'row.names')
# rownames(all_data) <- all_data$org_id
# tail(all_data)

# Get orgs sorting from the tree
#tree <- read.tree(paste0(DATA_DIR, "181010_all_species.rRNA.tree"))
tree <- read.tree(paste0(DATA_DIR, "orgs_chel.tree"))
# "624782022.Actinobacteria.Tomitella"  => "624782022"
#tree$tip.label <- gsub("^(\\d+).+", "\\1", tree$tip.label)
common_orgs <- intersect(tree$tip.label, rownames(all_data))
sorted_orgs <- tree$tip.label[tree$tip.label %in% common_orgs]
all_data <- all_data[sorted_orgs,]
head(all_data)

# evalue => log10_evalue
# evalue_cols <- grep('^evalue_', colnames(all_data))
# get_log10 <- function(x) ifelse(x < 1e-100, 100, -log10(x))
# all_data <- all_data %>%
#   mutate_at(evalue_cols, replace_na, 1) %>%
#   mutate_at(evalue_cols, get_log10) %>%
#   rename_at(evalue_cols, funs(sub("evalue_", "", .)))
# head(all_data)

# Make the bchlI, bchlD, bchlH
# all_data <- all_data %>%
#   mutate(bchlI = pmax(bchI, chlI),
#          bchlD = pmax(bchD, chlD),
#          bchlH = pmax(bchH, chlH))

# Generate 'taxa' column
all_data <- all_data %>%
  mutate(taxa = ifelse(kingdom == 'Archaea', 'Archaea', phylum)) %>%
  mutate(taxa = ifelse(taxa %in% names(TAXA_COLS), taxa, 'Other')) %>%
  # Make it factor for proper sorting
  mutate(taxa = factor(taxa, levels = names(TAXA_COLS)))
rownames(all_data) <- all_data$dir_name

# Generate 'frameshift' column
all_data <- all_data %>%
  mutate(frameshift = case_when(
    num_M_minus > 0 ~ '-1',
    num_M_plus > 0 ~ '+1',
    TRUE ~ 'None'))

str(all_data)

###
# Heatmap

# taxa_ht
taxa_ht <- Heatmap(data.frame(Taxonomy = all_data$taxa),
                   name = "taxa",
                   col = TAXA_COLS,
                   width = unit(3, "mm"),
                   cluster_columns = FALSE,
                   cluster_rows = FALSE,
                   show_row_names = FALSE,
                   show_heatmap_legend = FALSE,
                   split = all_data$taxa,
                   row_title_rot = 0,
                   gap = unit(0, "mm"),
                   row_title_gp = gpar(col = TAXA_COLS, font = 2))
#taxa_ht

# fs_ha
fs_ha <- rowAnnotation(Frameshift = all_data$frameshift,
                       col = list(Frameshift = c('-1' = 'black', '+1' = 'red', 'None' = 'white')),
                       annotation_legend_param = list(Frameshift = list(title = "Frameshift\nin chlD gene")),
                       show_annotation_name = TRUE,
                       width = unit(1, "cm"))
# taxa_ht + fs_ha

# patho_ha: http://www.bioconductor.org/packages/release/bioc/vignettes/ComplexHeatmap/inst/doc/s4.heatmap_annotation.html#toc_14
# patho_rows <- which(all_data$pathogen == 1)
# patho_names <- all_data[patho_rows, 'name']
# patho_ha <- rowAnnotation(link = anno_mark(at = patho_rows, labels = patho_names),
#                           width = unit(1, "cm") + max_text_width(patho_names))
# taxa_ht + fs_ha + chlD_ht + patho_ha

# photo_ha
# photo_ha <- rowAnnotation(Chlorophyll = all_data$Photosynthetic,
#                           col = list(Chlorophyll = c('Yes' = 'green4', 'No' = 'white')),
#                           width = unit(1, "cm"),
#                           show_annotation_name = TRUE,
#                           show_legend = FALSE)

# chlIDH_ht
chlIDH_ht <- Heatmap(as.matrix(all_data[, MG_CHEL_NAMES]),
                     column_title = "Magnesium\nchelatase\ngenes",
                     col = EVALUE_COLORS,
                     heatmap_legend_param = list('title' = 'BLAST\n-log10(E-value)'),
                     cluster_columns = FALSE,
                     width = unit(3*HT_COL_W, "cm"))
# taxa_ht + fs_ha + chlIDH_ht

# chl_path_ht
chl_path_ht <- Heatmap(as.matrix(all_data[, CHL_PATH_GENES]),
                       column_title = "Chlorophyll\nbiosynthesis\ngenes",
                       col = EVALUE_COLORS,
                       width = unit(6*HT_COL_W, "cm"),
                       cluster_columns = FALSE,
                       cluster_rows = FALSE,
                       show_heatmap_legend = FALSE,
                       show_row_names = FALSE)
# taxa_ht + fs_ha + chlIDH_ht + chl_path_ht

# b12_ha
# b12_ha <- rowAnnotation(Cobalamin = all_data$B12,
#                         col = list(Cobalamin = c('Yes' = 'blue', 'No' = 'white')),
#                         width = unit(1, "cm"),
#                         show_annotation_name = TRUE,
#                         show_legend = FALSE)

# cobNST_ht
cobNST_ht <- Heatmap(as.matrix(all_data[, COB_CHEL_NAMES]),
                     column_title = "Cobalt\nchelatase\ngenes",
                     col = EVALUE_COLORS,
                     cluster_columns = FALSE,
                     show_heatmap_legend = FALSE,
                     width = unit(3*HT_COL_W, "cm"))
# taxa_ht + fs_ha + chlIDH_ht + chl_path_ht + cobNST_ht


# b12_path_ht
b12_path_ht <- Heatmap(as.matrix(all_data[, B12_PATH_GENES]),
                       column_title = "Vitamin B12\nbiosynthesis\ngenes",
                       col = EVALUE_COLORS,
                       width = unit(6*HT_COL_W, "cm"),
                       cluster_columns = FALSE,
                       cluster_rows = FALSE,
                       show_heatmap_legend = FALSE,
                       show_row_names = FALSE)
# taxa_ht + fs_ha + chlIDH_ht + chl_path_ht + cobNST_ht + b12_path_ht

###
# Save heatmaps ----
pdf(paste0(OUT_DIR, 'heatmap_full.pdf'), width = 10, height = 8)
#draw(taxa_ht + fs_ha + chlIDH_ht + chl_path_ht + photo_ha + cobNST_ht + b12_path_ht + b12_ha,
draw(taxa_ht + fs_ha + chlIDH_ht + chl_path_ht + cobNST_ht + b12_path_ht,
     #show_heatmap_legend = FALSE,
     row_title = sprintf('%d prokaryotic genomes', nrow(all_data)))
dev.off()
