
source('lib.R')

library(ComplexHeatmap)
library(circlize)   # colorRamp2()

###

TAR_ORGS <- c(
  'Delftia acidovorans SPH-1',
  'Methanocaldococcus fervens AG86',
  'Methanocaldococcus sp. FS406-22',
  'Rhodobacter sphaeroides 2.4.1',
  'Mycobacterium tuberculosis H37Rv',
  'Nostoc punctiforme PCC 73102')

MG_CHEL_NAMES  <- c('chlH_bchH', 'chlD_bchD', 'chlI_bchI')
COB_CHEL_NAMES <- c('cobN', 'cobT', 'cobS')

CHL_PATH_GENES <- c("bchE", "chlB_bchB",  "chlG_bchG", "chlL_bchL", "chlM_bchM", "chlN_bchN")
B12_PATH_GENES <- c("cobD_cobC", "cobO", "cobP_cobU", "cobQ", "cobV_cobS", "cysG_cobA")

EVALUE_COLORS <- colorRamp2(c(0, 6, 100), c("white", "yellow", "red"))

HT_COL_W = 0.8
HT_RECT_GP = gpar(col = "black", lwd = 2)

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
head(all_data)

# Make the bchlI, bchlD, bchlH
# all_data <- all_data %>%
#   mutate(bchlI = chlI_bchI,
#          bchlD = chlD_bchD,
#          bchlH = chlH_bchH)

# Generate 'taxa' column
# most_frequent <- names(TAXA_COLS)
# all_data <- all_data %>%
#   mutate(taxa = ifelse(kingdom == 'Archaea', 'Archaea', phylum)) %>%
#   mutate(taxa = ifelse(taxa %in% most_frequent, taxa, 'Other')) %>%
#   # Make it factor for proper sorting
#   mutate(taxa = factor(taxa, levels = most_frequent))
# rownames(all_data) <- all_data$org_id

# Generate 'frameshift' column
all_data <- all_data %>%
  mutate(frameshift = ifelse(num_M_minus > 0, '-1', 'None')) %>%
  mutate(frameshift = ifelse(num_M_plus > 0, '+1', frameshift))

all_data <- filter(all_data, name %in% TAR_ORGS) %>%
  column_to_rownames('dir_name')

str(all_data)

###
# Heatmap

# chlIDH_ht
chlIDH_ht <- Heatmap(as.matrix(all_data[, MG_CHEL_NAMES]),
                     column_title = "Magnesium\nchelatase\ngenes",
                     col = EVALUE_COLORS,
                     heatmap_legend_param = list(
                       'title' = 'tBLASTn -log10(E-value)',
                       direction = "horizontal"),
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
                     row_split = -all_data$num_M_fs,
                     row_title = NULL,   # Hide cluster names
                     rect_gp = HT_RECT_GP,
                     row_names_side = "left",
                     row_labels = all_data$name,
                     height = unit(nrow(all_data)*HT_COL_W, "cm"),
                     width = unit(3*HT_COL_W, "cm"))
#draw(chlIDH_ht, heatmap_legend_side = "bottom")

# fs_ha
fs_ha <- rowAnnotation(Frameshift = all_data$frameshift,
                       col = list(Frameshift = c('-1' = 'black', '+1' = 'red', 'None' = 'white')),
                       annotation_legend_param = list(Frameshift = list(title = "Frameshift\nin chlD gene")),
                       show_annotation_name = TRUE,
                       width = unit(1, "cm"))

# chl_path_ht
chl_path_ht <- Heatmap(as.matrix(all_data[, CHL_PATH_GENES]),
                       column_title = "Chlorophyll\nbiosynthesis\ngenes",
                       col = EVALUE_COLORS,
                       cluster_columns = FALSE,
                       cluster_rows = FALSE,
                       show_heatmap_legend = FALSE,
                       show_row_names = FALSE,
                       rect_gp = HT_RECT_GP,
                       width = unit(6*HT_COL_W, "cm"))
# chlIDH_ht + chl_path_ht

# cobNST_ht
cobNST_ht <- Heatmap(as.matrix(all_data[, COB_CHEL_NAMES]),
                     column_title = "Cobalt\nchelatase\ngenes",
                     col = EVALUE_COLORS,
                     cluster_columns = FALSE,
                     show_heatmap_legend = FALSE,
                     show_row_names = FALSE,
                     rect_gp = HT_RECT_GP,
                     width = unit(3*HT_COL_W, "cm"))
# chlIDH_ht + chl_path_ht + cobNST_ht


# b12_path_ht
b12_path_ht <- Heatmap(as.matrix(all_data[, B12_PATH_GENES]),
                       column_title = "Vitamin B12\nbiosynthesis\ngenes",
                       col = EVALUE_COLORS,
                       cluster_columns = FALSE,
                       cluster_rows = FALSE,
                       show_heatmap_legend = FALSE,
                       show_row_names = FALSE,
                       rect_gp = HT_RECT_GP,
                       width = unit(6*HT_COL_W, "cm"))
# chlIDH_ht + chl_path_ht + cobNST_ht + b12_path_ht

###
# Save heatmaps ----
pdf(paste0(OUT_DIR, 'heatmap_Da_Mf_Ms.pdf'), width = 9.5, height = 4.5)
draw(chlIDH_ht + chl_path_ht + cobNST_ht + b12_path_ht + fs_ha,
     heatmap_legend_side = "bottom")
dev.off()
