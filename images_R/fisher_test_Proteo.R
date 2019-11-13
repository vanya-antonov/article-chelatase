
source("lib.R")

###

orgs_df <- read.delim(paste0(DATA_DIR, "orgs_chel.txt"), as.is=TRUE)
head(orgs_df)

###
# Test presense of a frameshift mutation

# Genomes with exactly 1 M-subunit and at least 1 L-subunit
total <- filter(orgs_df, phylum == 'Proteobacteria', num_M == 1, num_L > 0) %>% nrow()

# fs-chlD = Y; small absent = Y
YY <- filter(orgs_df, phylum == 'Proteobacteria', num_M == 1, num_L > 0, num_M_fs >  0, num_S == 0) %>% nrow()
YN <- filter(orgs_df, phylum == 'Proteobacteria', num_M == 1, num_L > 0, num_M_fs >  0, num_S >  0) %>% nrow()
NY <- filter(orgs_df, phylum == 'Proteobacteria', num_M == 1, num_L > 0, num_M_fs == 0, num_S == 0) %>% nrow()
NN <- filter(orgs_df, phylum == 'Proteobacteria', num_M == 1, num_L > 0, num_M_fs == 0, num_S >  0) %>% nrow()

mtx <- matrix(c(YY, YN,
                NY, NN),
              ncol=2, byrow=TRUE)
fisher.test(mtx)


###
# Test the presense of a slippery sites (from 1113.Ba.fisher_test_Proteo)

# chlD contains slippery site = Y; reduced genotype = Y
site_YY <- 60
site_NY <- 101 - 60
site_YN <- 2
site_NN <- 112 - 2

site_mtx <- matrix(c(site_YY, site_YN,
                     site_NY, site_NN),
                   ncol=2, byrow=TRUE)
fisher.test(site_mtx)$p.value

