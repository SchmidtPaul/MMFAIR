setwd("D:/Coding/MMFAIR")
library("rmarkdown")

rmarkdown::clean_site()  # delete old files
rmarkdown::render_site(encoding="UTF-8") # render all files new; UTF-8 for ä, ö, ü, ß

# get r files from rmd files
knitr::purl("heterogeneous_error_variance.Rmd",documentation=0)
