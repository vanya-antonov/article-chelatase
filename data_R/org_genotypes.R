
# Get the con object (DBI::dbConnect)
source("../data_R/private.R")

# Load the DATA_DIR
source("../images_R/lib.R")

###

# num_cobN, num_chlH, num_bchH,
# num_cobT, num_chlD, num_bchD,
# num_cobS, num_chlI, num_bchI

SQL <- "
select o.id, o.dir_name, o.name, o.genus, o.phylum, o.kingdom, tax2,
  num_L, num_M_tree AS num_M, num_S_tree AS num_S,
  num_M_fs, num_M_plus, num_M_minus
from chel_orgs_v o
where num_M_tree > 0;
"
orgs_all_df <- dbGetQuery(con, SQL)
write.table(orgs_all_df, file = paste0(DATA_DIR, "org_genotypes.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

