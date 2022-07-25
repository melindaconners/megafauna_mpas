# megafauna_mpas

#### 'a1_gridUD_loop'
This is the main script for running the gridded UD analysis on the multi-species dataset. The script loops through each species dataset and calls the following three functions (found in the 'functions' folder): 
- fn_HomeRangePub.R: calculates home range and core areas
- fn_add_weights.R: weights location data based on movement patterns of each species
- fn_coreUD_SP.R: called by fn_HomeRangePub.R if the user selects the "Seaman-Powell" method for calculating core area use
