
# Get the con object (DBI::dbConnect)
source("../data_R/private.R")

# Load the DATA_DIR
source("../images_R/lib.R")

###


SQL <- "
select name, descr, gene, prot_len, fs_len, org_name, kingdom, phylum, genus
from chel_feats_v
where chel_subunit is not NULL"
data_df <- dbGetQuery(con, SQL)
write.table(data_df, file = paste0(DATA_DIR, "genes_chel.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

