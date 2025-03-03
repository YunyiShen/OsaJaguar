library(rstan)
library(sf)
library(ggplot2)
library(reshape)
library(fields)
source("./R/util.R")

load("./res/scr_stan_fit12345_vb.rda")

z <- rstan::extract(m_fit, c("z"))$z
NN = rowSums(z)
jpeg("./res/Figs/popsize_est.jpg", 
        width = 6, height = 4, units = "in",
        res = 500)
par(mar = c(2.5,2.5,1.7,.5), mgp = c(1.5, 0.5, 0))
hist(NN, freq = TRUE, xlab = "Population size", main = "")
dev.off()

beta_env = rstan::extract(m_fit, c("betaenv"))$betaenv
jpeg("./res/Figs/env_beta.jpg", 
        width = 8 * 1.2, 
        height = 4 * 1.2, units = "in",res = 500)
par(mfrow = c(2,4),mar = c(2.5,2.5,1.7,.5), mgp = c(1.5, 0.5, 0))

plot(density(beta_env[,1]), main = "Conservation status", xlab = "",font.main = 1)
polygon(density(beta_env[,1]), col = "#9b9b9b")
abline(v=0, lwd = 2)

plot(density(beta_env[,2]), main = "Elevation", xlab = "",font.main = 1)
polygon(density(beta_env[,2]), col = "#9b9b9b")
abline(v=0, lwd = 2)

plot(density(beta_env[,3]), main = "Ruggness", xlab = "",font.main = 1)
polygon(density(beta_env[,3]), col = "#9b9b9b")
abline(v=0, lwd = 2)

plot(density(beta_env[,4]), main = "Human development", xlab = "",font.main = 1)
polygon(density(beta_env[,4]), col = "#9b9b9b")
abline(v=0, lwd = 2)

plot(density(beta_env[,5]), main = "Economical forest", xlab = "",font.main = 1)
polygon(density(beta_env[,5]), col = "#9b9b9b")
abline(v=0, lwd = 2)

plot(density(beta_env[,6]), main = "Wetland, mangrove", xlab = "",font.main = 1)
polygon(density(beta_env[,6]), col = "#9b9b9b")
abline(v=0, lwd = 2)

plot(density(beta_env[,7]), main = "Grassland", xlab = "",font.main = 1)
polygon(density(beta_env[,7]), col = "#9b9b9b")
abline(v=0, lwd = 2)

dev.off()

jpeg("./res/Figs/density_est.jpg", 
        width = 6, height = 6, units = "in",
        res = 500)
par(mar = c(.5,.5,1.7,.5), mgp = c(1.5, 0.5, 0))
range_pts = apply(stan_data$grid_pts, 2, max) - apply(stan_data$grid_pts, 2, min)
s = rstan::extract(m_fit, c("s"))$s
density_est <- SCRdensity(s,z,stan_data$grid_pts,
                    nx = ceiling(range_pts[1] * 10/3) + 0, 
                    ny = ceiling(range_pts[2] * 10/3) + 0, 
                    main = "density estimate", 
                    Xl = min(stan_data$grid_pts[,1])-.3, 
                    Xu = max(stan_data$grid_pts[,1])+.3,
                    Yl = min(stan_data$grid_pts[,2])-.3, 
                    Yu = max(stan_data$grid_pts[,2])+.3,
                    plotit = FALSE
                    )

image(density_est$grid$xg, 
        density_est$grid$yg, 
        density_est$Dn, 
        zlim = c(0.,45), 
        xlab = "", ylab = "", 
        col = gray.colors(45, start = 0., 
                          end = 0.9, gamma = .4, rev = TRUE), 
        xaxt='n', yaxt='n',font.main = 1)

points(stan_data$X[stan_data$deployred>0,], pch = 1)
for(j in 1:14){ # 13 seen individuals
    points(stan_data$X[stan_data$yred[j,]>0,], pch = 19)
}

points(stan_data$grid_pts, pch = 20, 
         col = adjustcolor("red", alpha.f = 0.8), cex = 0.3)
dev.off()
#points(grid_objs$traplocs[rowSums(JS_stan_data$deploy[,,i])>0,], pch = 1)