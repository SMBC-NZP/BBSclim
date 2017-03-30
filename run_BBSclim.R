
httr::set_config(httr::config(ssl_verifypeer = 0L))

### Install BBSclim package
install.packages("devtools")
devtools::install_github('SMBC-NZP/BBSclim')

### Run analysis
BBSclim::run_BBSclim()



