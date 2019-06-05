
# Get the con object (DBI::dbConnect)
source("private.R")

# Load the DATA_DIR
source("../images_R/lib.R")

###

orgs_all_df <- dbGetQuery(con, "select id, kingdom, phylum, genus, name from orgs order by name")
write.table(orgs_all_df, file = paste0(DATA_DIR, "orgs_all.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

orgs_chel_df <- dbGetQuery(con, "select id, kingdom, phylum, genus, name, num_M_zero, num_M_minus, num_M_plus from chel_orgs_v order by name")
write.table(orgs_chel_df, file = paste0(DATA_DIR, "orgs_chel.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
