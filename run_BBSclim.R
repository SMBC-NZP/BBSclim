
httr::set_config(config(ssl_verifypeer = 0L))

### Install BBSclim package
install.packages("devtools")
devtools::install_github('SMBC-NZP/BBSclim')


BBSclim::run_BBSclim()



