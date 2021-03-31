# megafauna_mpas

## This repository contains the major scripts and functions for the MPA/marine megafauna project.
Note: (March 2021) This repository is active and incomplete but will contain the full set of scripts/functions by the publishing date.

#### 'a1_gridUD_loop'
This is the main script for running the gridded UD analysis on the multi-species dataset. The script loops through each species dataset and calls the following three functions (found in the 'functions' folder): /
1- fn_HomeRangePub.R: calculates home range and core areas
2- fn_add_weights.R: weights location data based on movement patterns of each species
3- fn_coreUD_SP.R: called by fn_HomeRangePub.R if the user selects the "Seaman-Powell" method for calculating core area use
