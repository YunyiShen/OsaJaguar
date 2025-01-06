library(terra)
library(sf)
traps = read.csv("Effort_Locations_Combined.csv")
mask = rast("./processed_raster/mask.tif")

traps_sf = st_as_sf(traps, coords = c("Trap.Station.Longitude",
                                      "Trap.Station.Latitude"), 
                    crs = 4326) |>
  st_transform(crs = crs(mask))




