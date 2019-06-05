
# install.packages('cowplot')
# install.packages('RMySQL')

library(ggplot2)
library(cowplot)   # To have ggplots side-by-side: plot_grid()

library(dplyr)
library(tidyr)     # gather() and spread() functions
library(tibble)    # for rownames_to_column() and column_to_rownames()

library(RMySQL)

# library(ComplexHeatmap)
# library(circlize)   # colorRamp2()
# 
# library(ape)        # read.tree()
# library(phytools)   # findMRCA()


###

DATA_DIR <- "../data/"
OUT_DIR <- "../images/"

con <- DBI::dbConnect(RMySQL::MySQL(), 
                      host = "10.0.1.254",
                      # dbname = 'gtdb2',
                      dbname = 'chelatase_db',
                      user = "ivan",
                      password = "123")

theme_set(theme_bw(base_size = 19))  # increase the font size: https://stackoverflow.com/a/11955412/310453

