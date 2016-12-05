
code_lookup <- read.csv("data-raw/code_lookup.csv")

code_lookup <- select(code_lookup, -group)

devtools::use_data(code_lookup, internal = TRUE)
