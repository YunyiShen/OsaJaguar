library(terra)
library(sf)

tri_function <- function(x) {
  center <- x[ceiling(length(x) / 2)]  # Get the center cell
  sum(abs(center - x), na.rm = TRUE)  # Sum of absolute differences
}


#### deal with human influence index ####
ele = rast("./rawraster/SRTM_30m_31971_AmistOsa.tif")
ghm = rast("./rawraster/gHM.tif")

mask_for_ghm = project(ele, crs(ghm))
ghm = resample(ghm, mask_for_ghm)
ghm = project(ghm, crs(ele))

writeRaster(ghm, "./processed_raster/gHM.tif")

#### deal with elevation, ruggness and lclu ####
lclu = rast("./rawraster/lulc_2017.tif")
ele = rast("./rawraster/SRTM_30m_31971_AmistOsa.tif")
rugg = focal(ele, w = matrix(1, 3, 3), fun = tri_function, na.policy = "omit")

lclu = project(lclu, crs(ele), method = "near", res = c(1000, 1000))
ele = resample(ele, lclu)
rugg = resample(rugg, lclu)

writeRaster(ele, "./processed_raster/elevation.tif")
writeRaster(rugg, "./processed_raster/ruggedness.tif")
writeRaster(lclu, "./processed_raster/lulc.tif")

#### get grid points #### 

cell_indices <- cells(lclu)

# Get the x and y coordinates of pixel centers
pixel_centers <- xyFromCell(lclu, cell_indices) |> as.data.frame()
names(pixel_centers) <- c("x","y")
pixel_centers_sf = st_as_sf(pixel_centers, coords = c("x", "y"), crs = crs(ele))
write.csv(pixel_centers, "./gridpoints.csv")



########## connect envvirables to pixel points ########

pixel_centers = read.csv("./gridpoints.csv", row.names = 1)
ele = rast("./processed_raster/elevation.tif")
rugg = rast("./processed_raster/ruggedness.tif")
lclu = rast("./processed_raster/lulc.tif")
# LCLU code 1-9: palm plantation, Mangrove, Water, 
#   Grassland, Urban, Old growth forest, Secondary forest, 
#   Wetland, Teak
ghm = rast( "./processed_raster/gHM.tif")

full_mask = (ele>-10000) * (lclu >= 0)
plot(full_mask)
writeRaster(full_mask, "./processed_raster/mask.tif")

# reserves ####
reserve = st_read("./rawraster/ManagedLandsCR.shp")
reserve = st_transform(reserve, crs(ele))



inside = st_intersects(pixel_centers_sf, reserve)
inside = sapply(inside, function(w){1*(length(w)>0)})

pixel_centers$conservation = inside
pixel_centers$ele = extract(ele, pixel_centers_sf)[,2]
pixel_centers$rugg = extract(rugg, pixel_centers_sf)[,2]
pixel_centers$lclu = extract(lclu, pixel_centers_sf)[,2]
pixel_centers$ghm = extract(ghm, pixel_centers_sf)[,2]

## group lclu 
## 1+9: economical plants
## 2+8: wetland and mangrove
## 4: grassland
## 5: urban 
## 6+7: forest, use as reference

pixel_centers$econforest = (pixel_centers$lclu == 1 | pixel_centers$lclu == 9) * 1
pixel_centers$wetland_mangrove = (pixel_centers$lclu == 2 | pixel_centers$lclu == 8) * 1
pixel_centers$grassland = (pixel_centers$lclu == 4) * 1
pixel_centers$urban = (pixel_centers$lclu == 5) * 1


pixel_centers_filtered <- as.data.frame(pixel_centers) |> na.omit()
pixel_centers_filtered <- pixel_centers_filtered[pixel_centers_filtered$lclu!=3,] # filter out water

#na.omit(pixel_centers_filtered)
plot(ele)
points(y~x, data = pixel_centers_filtered)

write.csv(pixel_centers_filtered, "./grid_points_secr.csv", row.names = F)

##### filter traps ######




