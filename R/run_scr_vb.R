## if rerun pre processing ###
#source("./R/gis.R")
#source("./R/handling_traps.R")
#source("./R/get_detection_hist.R")
#source("./R/form_grid.R")
#source("./R/prepare_scr.R")

library(rstan)
library(jsonlite)
config <- fromJSON("./config.json")
load("./processed_data/grid_obj.rda")
load("./processed_data/jaguar_trap_mats_js.rda")

options(mc.cores = 10)
rstan_options(auto_write = TRUE)

load("./processed_data/scr_stan_data.rda")
m_init <- stan_model("./stan/secr_simpler.stan")
stan_data$everdetected <- as.integer(stan_data$everdetected)

set.seed(12345)
m_fit <- vb(m_init,  data = stan_data, 
                    iter = config$stan_iters, 
                  init = function() list(log_psi = log(c(0.15, 0.15)), 
                                         log_p0 = log(c(.2, 0.2)), 
                                         sigma = log(c(1.,1.)),
                                         beta = rep(0,stan_data$n_env)), 
                  sample_file = "./res/scr_stan_fit12345_vb.csv")

save(m_fit, jaguar_trap_mats, stan_data, config,
     file="./res/scr_stan_fit12345_vb.rda")