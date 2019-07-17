
# Get the con object (DBI::dbConnect)
source("../data_R/private.R")

# Load the DATA_DIR
source("../images_R/lib.R")

###

# num_cobN, num_chlH, num_bchH,
# num_cobT, num_chlD, num_bchD,
# num_cobS, num_chlI, num_bchI

SQL <- "
select o.id, o.name, o.genus, o.phylum, o.kingdom, tax2,
  num_L, num_M_tree, num_M_fs,
  (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_subunit_tree='S' and t.num_kids=0) AS num_S_tree_wo_kids
from chel_orgs_v o
where num_M_tree > 0;
"
orgs_all_df <- dbGetQuery(con, SQL)
write.table(orgs_all_df, file = paste0(DATA_DIR, "org_genotypes.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

