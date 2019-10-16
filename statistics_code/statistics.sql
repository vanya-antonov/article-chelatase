
# "Among the 38 fs-chlD genes from Actinobacteria"
select * from chel_feats_v where phylum = 'Actinobacteria' and fs_len in (1, -1);

# "there were 117 "1 x cobN, 1 x chlD" genomes"
select * from chel_orgs_v where phylum = 'Actinobacteria' and genotype = '1xcobN, 1xchlD';


# "More than half of the actinobacterial genomes with this genotype were from the Streptomyces genera."
select * from chel_orgs_v where genus = 'Streptomyces' and genotype = '1xcobN, 1xchlD';

# "out of the 25 archaeal genomes with frameshifted chlD genes, 23 had the "N x L, 1 x fs-chlD" genotype"
select * from chel_orgs_v where kingdom = 'Archaea' and num_M_fs > 0;
select * from chel_orgs_v where kingdom = 'Archaea' and num_M_fs > 0 and num_L > 0 and num_S=0;

