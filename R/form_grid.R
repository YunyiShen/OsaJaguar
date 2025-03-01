library(sf)
library(terra)
library(jsonlite)
load("./processed_data/traps_loc.rda")

config <- jsonlite::fromJSON("config.json")
scaling <- config$scaling

grid = read.csv("./processed_data/grid_points_secr.csv")
grid_pts = as.matrix(grid[,c("x","y")]/scaling)
trap_X = st_coordinates(traps_sf)/scaling

env_covar = as.matrix(grid[,c(3,4,5,7,8,9,10)])
env_covar = scale(env_covar)

distsqr = matrix(0, nrow = nrow(grid_pts), ncol = nrow(trap_X))
for (i in 1:nrow(grid_pts)){
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
