#!/usr/bin/Rscript
library(httr)
set_config(config(ssl_verifypeer = 0L))

### Install BBSclim package
#install.packages("devtools")

#devtools::install_github('crushing05/BBS.tenstop')
#devtools::install_github('crushing05/BBS.fiftystop')
#devtools::install_github('richfitz/datastorr')
#devtools::install_github('crushing05/rBBS')
#devtools::install_github('richfitz/remake')
#remake::install_missing_packages("remake.yml")
#devtools::install_github('krlmlr/here')

#devtools::install_github('SMBC-NZP/BBSclim')


library(BBSclim)

model_opts(tenstops = TRUE, limit.cores = 50, psi.test = FALSE, gam.test = FALSE)

spp <- read.csv("inst/spp_list.csv")
spp_names <- as.character(spp$spp)
vals <- list(spp_names = whisker::iteratelist(spp_names, value="spp_name"))

str <- whisker::whisker.render(readLines("config/remake.yml.whisker"), vals)
writeLines(str, "remake.yml")

remake::make()

