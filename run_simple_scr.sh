#!/bin/bash
module load anaconda/2022a
source activate Rstan

Rscript ./R/run_scr_vb.R
Rscript ./R/run_scr.R