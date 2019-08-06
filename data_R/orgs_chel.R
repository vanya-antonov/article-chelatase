
# Get the con object (DBI::dbConnect)
source("../data_R/private.R")

# Load the DATA_DIR
source("../images_R/lib.R")

###

SQL <- "
select id, dir_name, name, kingdom, phylum, genus,
  num_M_zero, num_M_minus, num_M_plus
from chel_orgs_v order by name;
"
orgs_chel_df <- dbGetQuery(con, SQL)
write.table(orgs_chel_df, file = paste0(DATA_DIR, "orgs_chel.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
