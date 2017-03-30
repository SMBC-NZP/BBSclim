# Set global options
model_opts(psi.test = TRUE, gam.test = TRUE)

## Format makefile given current species list
library(whisker)

spp <- read.csv("inst/spp_list.csv")
spp_names <- as.character(spp$spp)
vals <- list(spp_names = iteratelist(spp_names, value="spp_name"))

str <- whisker.render(readLines("config/remake.yml.whisker"), vals)
writeLines(str, "remake.yml")

remake::make()


