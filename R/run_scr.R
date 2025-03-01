library(rstan)
library(jsonlite)
config <- fromJSON("./config.json")

options(mc.cores = 4)
rstan_options(auto_write = TRUE)

load("./processed_data/scr_stan_data.rda")
m_init <- stan_model("./stan/secr.stan")