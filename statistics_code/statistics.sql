
###
# General statistics

# "In total, there were 135 and 36 chlD genes with predicted -1 and +1 frameshifts, respectively"
select fs_len, count(*) from chel_feats_v where fs_len in (1, -1) group by fs_len;


# "Among the 38 fs-chlD genes from Actinobacteria, 27 were present in the genomes with the "1 x cobN, 1 x fs-chlD" genotype"
select * from chel_feats_v where phylum = 'Actinobacteria' and fs_len in (1, -1);
select * from chel_orgs_v where phylum = 'Actinobacteria' and genotype = '1xcobN, 1xfs-chlD';

# "Namely, there were 29 and 9 fs-chlD genes with +1 and -1 frameshifts, respectively"
select fs_len, count(*) from chel_feats_v where phylum = 'Actinobacteria' and fs_len in (1, -1) group by fs_len;

# "there were 139 "1 x cobN, 1 x chlD" genomes"
select count(*) from chel_orgs_v where phylum = 'Actinobacteria' and genotype = '1xcobN, 1xchlD';

# "More than half of the actinobacterial genomes with this genotype were from the Streptomyces genera."
select count(*) from chel_orgs_v where genus = 'Streptomyces' and genotype = '1xcobN, 1xchlD';

# "out of the 25 archaeal genomes with frameshifted chlD genes, 23 had the "N x L, 1 x fs-chlD" genotype"
select * from chel_orgs_v where kingdom = 'Archaea' and num_M_fs > 0;
select * from chel_orgs_v where kingdom = 'Archaea' and num_M_fs > 0 and num_L > 0 and num_S=0;


###
# Cyanobacteria

# "In 75 out of the 76 analyzed genomes the subunits of the magnesium chelatase were encoded by the full set of genes"
select * from chel_orgs_v where phylum = 'Cyanobacteria';

select * from chel_orgs_v where phylum = 'Cyanobacteria'
and num_chlH+num_bchH > 0 and num_chlD+num_bchD > 0 and num_chlI+num_bchI > 0;

# "Interestingly, more than half of the cyanobacterial genomes (43 out of 76) contained two chlH genes"
select * from chel_orgs_v where phylum = 'Cyanobacteria'
and num_chlH+num_bchH = 2 and num_chlD+num_bchD = 1 and num_chlI+num_bchI = 1;


###
# Archaea

# "Namely, out of the 25 archaeal genomes with frameshifted chlD genes, 22 had the "N x L, 1 x fs-chlD" genotype where L is the large subunit gene"
select * from chel_orgs_v where kingdom = 'Archaea' and num_M_fs > 0;
select * from chel_orgs_v where kingdom = 'Archaea' and num_M_fs > 0 and num_M_zero=0 and num_L > 0 and num_S = 0;


select * from chel_orgs_v where kingdom = 'Archaea' and num_M_fs = 0;
select * from chel_orgs_v where kingdom = 'Archaea' and num_M_fs = 0 and num_L > 0 and num_S > 0;


select genotype, count(*) as N from chel_orgs_v where kingdom = 'Archaea' and num_M_fs > 0
group by genotype order by N desc;

select genotype, count(*) as N from chel_orgs_v where kingdom = 'Archaea' and num_M_fs > 0
group by genotype order by N desc;

