load("./processed_data/grid_obj.rda")
load("./processed_data/jaguar_trap_mats_js.rda")

cap_mat = jaguar_trap_mats$cap_mat
trap_mat = jaguar_trap_mats$trap_mat
ind_ids = jaguar_trap_mats$ids$ind_ids


N = dim(cap_mat)[1]
n_trap = dim(cap_mat)[2]

# sum over all occasions
y_red = apply(cap_mat, c(1,2), sum)
ever_detected = rowSums(y_red) > 0

augment_sex = rep(0, N - nrow(ind_ids))
augment_sex[1:floor(length(augment_sex)/2)] = 1
sex = c(ind_ids$Sexidx, augment_sex) # the rest being half half
deployred = rowSums(trap_mat)

stan_data = grid_obj
stan_data$y_red = y_red
stan_data$everdetected = ever_detected
stan_data$sex = sex
stan_data$deployred = deployred
stan_data$N = N
stan_data$n_trap = n_trap

save(stan_data, file = "./processed_data/scr_stan_data.rda")
