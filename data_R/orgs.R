
# Get the con object (DBI::dbConnect)
source("private.R")

# Load the DATA_DIR
source("../images_R/lib.R")

###

SQL <- "
select o.id, o.name, o.genus, o.phylum, o.kingdom,
    (select t.value from org_params t where t.parent_id=o.id and t.name='taxonomy' and t.num=2) AS tax2,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_subunit='M') AS num_M,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_gene in ('chlD', 'bchD')) AS num_bchlD,
    (select count(distinct t.id) from chel_feats_v t where t.org_id=o.id and t.chel_gene in ('chlD', 'bchD') and t.fs_len in (-1, 1)) AS num_bchlD_fs
from orgs o;
"
orgs_all_df <- dbGetQuery(con, SQL)
write.table(orgs_all_df, file = paste0(DATA_DIR, "orgs_all.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

orgs_chel_df <- dbGetQuery(con, "select id, kingdom, phylum, genus, name, num_M_zero, num_M_minus, num_M_plus from chel_orgs_v order by name")
write.table(orgs_chel_df, file = paste0(DATA_DIR, "orgs_chel.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
