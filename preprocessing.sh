module load anaconda/2022a
source activate Rstan

Rscript ./R/gis.R
Rscript ./R/handling_traps.R
Rscript ./R/get_detection_hist.R
Rscript ./R/form_grid.R
Rscript ./R/prepare_scr.R