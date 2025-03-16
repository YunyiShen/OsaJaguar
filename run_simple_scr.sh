#!/bin/bash
module load anaconda/2023a
source activate Rstan

Rscript ./R/gis.R
Rscript ./R/handling_traps.R
Rscript ./R/get_detection_hist.R
Rscript ./R/form_grid.R
Rscript ./R/prepare_scr.R
#Rscript ./R/run_scr_vb.R
Rscript ./R/run_scr.R
Rscript ./R/get_density_plot.R