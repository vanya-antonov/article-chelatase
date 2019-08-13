
# Get the con object (DBI::dbConnect)
source("../data_R/private.R")

# Load the DATA_DIR
source("../images_R/lib.R")

###

genes_df <- read.delim(
  url('https://raw.githubusercontent.com/vanya-antonov/django_gtdb2/master/chelatase_db/data/pathway_genes.txt'),
  comment.char = '#', header = TRUE, as.is = TRUE)
#cat(unique(sort(genes_df$chel_gene_group)), sep = '\n')

sql_str <- "
select 
  (select t.value from org_params t where t.parent_id=o.id and t.name='dir_name') AS dir_name,

  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='chlI_bchI' and num_kids=0) AS evalue_chlI_bchI,
  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='chlD_bchD') AS evalue_chlD_bchD,
  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='chlH_bchH') AS evalue_chlH_bchH,

  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='cobN') AS evalue_cobN,
  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='cobS') AS evalue_cobS,
  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='cobT') AS evalue_cobT,

  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='bchE') AS evalue_bchE,
  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='chlB_bchB') AS evalue_chlB_bchB,
  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='chlG_bchG') AS evalue_chlG_bchG,
  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='chlL_bchL') AS evalue_chlL_bchL,
  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='chlM_bchM') AS evalue_chlM_bchM,
  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='chlN_bchN') AS evalue_chlN_bchN,
  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='cobD_cobC') AS evalue_cobD_cobC,
  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='cobO') AS evalue_cobO,
  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='cobP_cobU') AS evalue_cobP_cobU,
  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='cobQ') AS evalue_cobQ,
  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='cobV_cobS') AS evalue_cobV_cobS,
  (select min(t.chel_evalue) from chel_feats_v t where t.org_id=o.id and t.chel_gene_group='cysG_cobA') AS evalue_cysG_cobA
from chel_orgs_v o;
"

evalue_df <- dbGetQuery(con, sql_str)
evalue_df[is.na(evalue_df)] <- 1
head(evalue_df)
write.table(evalue_df, file = paste0(DATA_DIR, "orgs_evalue.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

