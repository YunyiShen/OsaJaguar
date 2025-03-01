library(sf)
library(jsonlite)

load("./processed_data/traps_loc.rda")
load("./processed_data/detections.rda")
config <- jsonlite::fromJSON("config.json")

date_range <- config$date |> as.Date() # setup a date range
occ_days <- config$occasion # number of days in a secondary occasion
M <- config$max_individual # number of augmented individuals

year_range <- format(date_range,"%Y")

## handle occasions
traps_sf$Setup_date = as.Date(traps_sf$Start_Date, format = "%m/%d/%Y")
traps_sf$Retrieval_date = as.Date(traps_sf$End_Date, format = "%m/%d/%Y")
overall_start_date = min(traps_sf$Setup_date)
traps_sf$Setup_occasion = floor(as.numeric(traps_sf$Setup_date - overall_start_date)/occ_days) + 1
traps_sf$Retrieval_occasion = floor(as.numeric(traps_sf$Retrieval_date - overall_start_date)/occ_days) + 1

## handle detection occasions
detections_sf$Date = as.Date(detections_sf$Date.Time, format = "%m/%d/%y %H:%M")
detections_sf$Occasion = floor(as.numeric(detections_sf$Date - overall_start_date)/occ_days) + 1

## enumeration of traps ## 
trap_ids <- data.frame(unique(as.data.frame(traps_sf)[,c("Station.Name","x","y")]))
trap_ids$id <- 1:nrow(trap_ids)

ind_ids <- data.frame(unique(as.data.frame(detections_sf)[,c("Individual.ID", "Sex")]))
ind_ids$id <- 1:nrow(ind_ids)
ind_ids$Sexidx <- as.numeric(ind_ids$Sex == "male")

K_max = max(traps_sf$Retrieval_occasion)
n_trap = nrow(trap_ids)

# whether trap is active at each occasion
trap_mat = matrix(0, ncol = K_max, nrow = n_trap)

for (i in 1:n_trap){
    station = trap_ids$Station.Name[i]
    tmp = traps_sf[traps_sf$Station.Name == station,]
    for (j in 1:nrow(tmp)){
        trap_mat[i, tmp$Setup_occasion[j]:tmp$Retrieval_occasion[j]] = 1
    }
}

# whether individual is detected at each occasion
cap_mat <- array(0,c(M, n_trap, K_max))

for(i in 1:nrow(detections_sf)){
    ind_id = detections_sf$Individual.ID[i]
    ind_idx = which(ind_ids$Individual.ID == ind_id)
    trap_idx = which(trap_ids$Station.Name == detections_sf$Station.Name[i])
    occ_idx = detections_sf$Occasion[i]
    if (occ_idx > K_max){
        next
    }
    cap_mat[ind_idx, trap_idx, occ_idx] = 1
}

jaguar_trap_mats <- list(cap_mat = cap_mat, trap_mat = trap_mat,
                    ids = list(trap_ids = trap_ids, ind_ids = ind_ids),
                    occa_design = list(start_date = overall_start_date
                                       )
                    )

save(jaguar_trap_mats, file = "./processed_data/jaguar_trap_mats_js.rda")