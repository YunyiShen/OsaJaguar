library(rstan)
library(sf)
library(ggplot2)
library(reshape)
source("./R/util.R")

postfix = "_withTico_2km"
load( paste0("./res/scr_stan_fit12345_2024",postfix,".rda"))

z <- rstan::extract(m_fit, c("z"))$z
NN = rowSums(z)
beta_env = rstan::extract(m_fit, c("betaenv"))$betaenv
sex = stan_data$sex

females = z[,sex == 0]
males = z[,sex == 1]
NN_female = rowSums(females)
NN_male = rowSums(males)

detection_rate = exp(rstan::extract(m_fit, "log_p0")$log_p0)
decay = exp(-rstan::extract(m_fit, "log_sigma")$log_sigma)



###### population size ######
jpeg( paste0("./res/Figs/popsize_est", postfix, ".jpg"), 
        width = 12/2, height = 4/2, units = "in",
        res = 500)
par(mar = c(2.5,2.5,1.7,.5), mgp = c(1.5, 0.5, 0), mfrow = c(1,3))
hist(NN, breaks = seq(min(NN) - 0.5, max(NN) + 0.5, by = 1),
        freq = FALSE, xlab = "Population size", main = "", ylab = "posterior")
hist(NN_female, breaks = seq(min(NN_female) - 0.5, max(NN_female) + 0.5, by = 1),
        freq = FALSE, xlab = "Female", main = "", ylab = "")
hist(NN_male, breaks = seq(min(NN_male) - 0.5, max(NN_male) + 0.5, by = 1),
        freq = FALSE, xlab = "Male", main = "", ylab = "")
dev.off()


###### environemnt #####
jpeg(paste0("./res/Figs/env_beta", postfix, ".jpg"), 
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


####### density #######
jpeg(paste0("./res/Figs/density_est", postfix, ".jpg"), 
        width = 8, height = 6, units = "in",
        res = 500)
par(mar = c(.5,.5,1.7,.5), mgp = c(1.5, 0.5, 0))
#range_pts = apply(stan_data$grid_pts, 2, max) - apply(stan_data$grid_pts, 2, min)
s = rstan::extract(m_fit, c("s"))$s
density_est = SCRdensity_tiff(s, z, stan_data$grid_pts)
writeRaster(density_est, paste0("./res/density_est", postfix, ".tiff"),
                overwrite=TRUE)


plot(density_est, col = gray.colors(30, start = 0., 
                          end = 0.9, gamma = 1., rev = TRUE),
                cex.axis = 1.5
                          )
mtext(expression("Density (100 km"^2*")"), side = 4, line = 3, cex = 1.2)

points(stan_data$X[stan_data$deployred>0,] * config$scaling, pch = 1, cex = 0.5)
for(j in 1:14){ # 13 seen individuals
    points(stan_data$X[stan_data$yred[j,]>0,] * config$scaling, pch = 19, cex = 0.5)
}

#points(stan_data$grid_pts * config$scaling, pch = 20, 
#         col = adjustcolor("red", alpha.f = 0.8), cex = 0.3)


legend("bottomright", 
       #inset=c(1.2, 0),
       legend = c("Camera traps",
                 "Jaguars detected"
                  ),
       pch = c(1,19), cex = c(1,1), 
       bty = "n",
       col = c("black","black"),
       y.intersp = 2, 
       xjust = 0
       )


dev.off()

##### other vital rates #####
jpeg(paste0("./res/Figs/vital_rates", postfix, ".jpg"), 
        width = 6 * 1.2, 
        height = 3 * 1.2, units = "in",res = 500)
par(mfrow = c(2,2),mar = c(2.5,2.5,1.7,.5), mgp = c(1.5, 0.5, 0))

plot(density(detection_rate[,1]), 
        main = "Detection rate", xlab = "",
        font.main = 1, ylab = "Female posterior",
        xlim = c(min(detection_rate), max(detection_rate))
        )
polygon(density(detection_rate[,1]), col = "#9b9b9b")

plot(density(decay[,1]), main = "Detection length scale", 
        xlab = "",font.main = 1, ylab = "",
        xlim = c(min(decay), max(decay))
        )
polygon(density(decay[,1]), col = "#9b9b9b")

plot(density(detection_rate[,2]), main = "", xlab = "",
        font.main = 1, ylab = "Male posterior",
        xlim = c(min(detection_rate), max(detection_rate))
        )
polygon(density(detection_rate[,2]), col = "#9b9b9b")

plot(density(decay[,2]), main = "", xlab = "",
                font.main = 1, ylab = "",
        xlim = c(min(decay), max(decay))
                )
polygon(density(decay[,2]), col = "#9b9b9b")
dev.off()