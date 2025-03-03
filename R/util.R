SCRdensity <- function(s, z, pts, plotit = TRUE,
                      nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL, 
                      Yu = NULL, scalein = 100, scaleout = 100, 
                      col="gray",ncolors = 10,whichguy=NULL,plot_scale = TRUE,calculate_sd = FALSE,...){

    Sxout <- pts[s,1] |> 
      matrix(nrow = nrow(s))

    Syout <- pts[s,2] |> 
      matrix(nrow = nrow(s))
  
  
  niter <- nrow(Sxout)
  if (is.null(Xl)) {
    Xl <- min(Sxout) * 0.999
    Xu <- max(Sxout) * 1.001
    Yl <- min(Syout) * 0.999
    Yu <- max(Syout) * 1.001
  }
  xg <- seq(Xl, Xu, , nx)
  yg <- seq(Yl, Yu, , ny)
  guy<-col(Sxout)
  
  
  #Sxout <- cut(Sxout, breaks = xg)
  #Syout <- cut(Syout, breaks = yg)
  if(is.null(whichguy)){
    area <- (yg[2] - yg[1]) * (xg[2] - xg[1]) * scalein
    if(calculate_sd){
      Dn_sd <- sapply(1:nrow(Sxout), function(i, Sxout, Syout, z, xg,yg, area, scalein){
        table(cut(Sxout[i,z[i,]==1], breaks = xg), 
            cut(Syout[i,z[i,]==1], breaks = yg))/area * scalein
      }, Sxout, Syout, z, xg,yg, area, scalein) |>
        apply(1,sd) |> 
        matrix(nrow = nx-1)
    }
    else Dn_sd <- NULL
    
    Sxout <- cut(Sxout[z == 1], breaks = xg)
    Syout <- cut(Syout[z == 1], breaks = yg)
    Dn <- table(Sxout, Syout)/niter
    
    Dn <- (Dn/area) * scaleout
  }
  else{
    Sxout <- cut(Sxout[z == 1], breaks = xg)
    Syout <- cut(Syout[z == 1], breaks = yg)
    Dn<-table(Sxout[guy==whichguy],Syout[guy==whichguy] )/niter
  }
  
  cat("mean: ", mean(Dn), fill = TRUE)
  if(plotit){
    par(mar = c(3, 3, 3, 6))
    if (col == "gray") {
      cc <- seq(3, 17, , 10)/20
      cc <- rev(gray(cc))
    }
    else cc <- terrain.colors(ncolors)
  
  
    image(xg, yg, Dn, col = cc,...)
    if(plot_scale) image.scale(Dn, col = cc)
    box()
  }

  return(list(grid = list(xg = xg, yg = yg), Dn = Dn, Dn_sd = Dn_sd))
}