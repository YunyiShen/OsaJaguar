library(terra)
library(sf)
traps = read.csv("./processed_data/Effort_Locations_Combined.csv")
mask = rast("./processed_data/processed_raster/mask.tif")


#### deal with trap locations ####
traps_sf = st_as_sf(traps, coords = c("Trap.Station.Longitude",
                                      "Trap.Station.Latitude"), 
                    crs = 4326) |>
  st_transform(crs = crs(mask))

traps_sf$x = st_coordinates(traps_sf)[,1]
traps_sf$y = st_coordinates(traps_sf)[,2]

# test if the traps are within the mask
traps_sf$mask = extract(mask, st_coordinates(traps_sf)) |> unlist()
traps_sf = traps_sf[!is.na( traps_sf$mask),]

save(traps_sf, file = "./processed_data/traps_loc.rda")

#### deal with detections ####
detections = read.csv("./rawdata/Detection_History_Combined.csv")
detections_sf = st_as_sf(detections, coords = c("Trap.Station.Longitude",
                                      "Trap.Station.Latitude"), 
                         crs = 4326) |>
  st_transform(crs = crs(mask))

detections_sf$x = st_coordinates(detections_sf)[,1]
detections_sf$y = st_coordinates(detections_sf)[,2]

# test if the detections are within the mask
detections_sf$mask = extract(mask, st_coordinates(detections_sf)) |> unlist()
detections_sf = detections_sf[!is.na( detections_sf$mask),]

unknown_ind = grepl("Unknown", detections_sf$Individual.ID) | grepl("Undetermined", detections_sf$Individual.ID)
detections_sf = detections_sf[!unknown_ind,]

save(detections_sf, file = "./processed_data/detections.rda")