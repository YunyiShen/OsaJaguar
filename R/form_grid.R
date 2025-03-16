library(jsonlite)
#load("./processed_data/traps_loc.rda")
load("./processed_data/jaguar_trap_mats_js.rda")

config <- jsonlite::fromJSON("config.json")
scaling <- config$scaling

grid = read.csv("./processed_data/grid_points_secr.csv")
grid_pts = as.matrix(grid[,c("x","y")]/scaling)
trap_X = jaguar_trap_mats$ids$trap_ids[,c("x","y")]/scaling

env_covar = as.matrix(grid[,c(3,4,5,7,8,9,10)])
env_covar[,2:4] = scale(env_covar[,2:4]) # scale continuous covariates
#env_covar = scale(env_covar) # scale all
print(nrow(grid_pts))
distsqr = matrix(0, nrow = nrow(grid_pts), ncol = nrow(trap_X))
for (i in 1:nrow(grid_pts)){
    print(i/nrow(grid_pts))
    for (j in 1:nrow(trap_X)){
        distsqr[i,j] = sum((grid_pts[i,] - trap_X[j,])^2)
    }
}

grid_obj = list(X = trap_X,
                n_grid = nrow(grid_pts),
                grid_pts = grid_pts, 
                distsqr = distsqr,
                n_env = ncol(env_covar), 
                envX = env_covar)

save(grid_obj, file = "./processed_data/grid_obj.rda")
