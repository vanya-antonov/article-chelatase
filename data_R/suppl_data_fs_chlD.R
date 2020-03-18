
# Get the con object (DBI::dbConnect)
source("../data_R/private.R")

# Load the DATA_DIR
#source("../images_R/lib.R")
library(dplyr)

DATA_DIR <- "../data/"
OUT_DIR <- "../images/"

###


SQL_fs_chlD <- "
select id, parent_id, org_name, seq_id, start, end, strand, name, descr, fs_len, translation AS long_product, seq_nt AS long_cds
from chel_feats_v
where fs_len in (-1, 1);"
data_df <- dbGetQuery(con, SQL_fs_chlD)
head(data_df)

SQL_short_prot <- "select id AS parent_id, translation AS short_product, seq_nt AS short_cds from chel_feats_v where num_kids > 0;"
short_prot_df <- dbGetQuery(con, SQL_short_prot)
head(short_prot_df)

right_join(data_df, short_prot_df, by = 'parent_id') %>%
  select(-id, -parent_id) %>%
  arrange(org_name) %>%
  write.table(file = paste0(DATA_DIR, "suppl_data_fs_chlD.txt"),
              row.names = FALSE, quote = FALSE, sep = "\t")
