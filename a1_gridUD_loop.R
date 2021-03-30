
rm(list = ls(all = TRUE))

# * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #
# Gridded Utilization Distributions
# Script Processes: this script will iterate over tracking datasets to construct gridded utilization distributions
# Written for Conners et al, XXXX
# contact: melinda.conners@stonybrook.edu
# Last updated: January 10, 2020
# * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #


# * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #
# Set Environment
# * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #

# Load Packages -------------------------------------------------------------------------

library(maptools)
data(wrld_simpl)
library(tidyverse)
library(lubridate)
library(sp)
library(plyr)
library(purrr)
library(adehabitatLT)
library(adehabitatHR)
library(spatstat)
library(pracma)
library(geosphere)
library(raster)


# Load Functions -------------------------------------------------------------------------

fundir <- "~/00 pew functions/"  
source(paste0(fundir,"AddWeightsFunction.R"))
source(paste0(fundir,"HomeRangeFunctionPub.R"))
source(paste0(fundir,"SPCoreFunction.R"))
source(paste0(fundir,"Lon180to360.R"))
source(paste0(fundir,"Lon360to180.R"))
source(paste0(fundir,"km2Degree.R"))

########################################################################################### 
# Begin Analytical Framework
########################################################################################### 

# Set Root Directory: This the main directory hub where files for this project live.
root<- '~/Documents/PewMarMega/analysis/'            

# Set Data Directory: This is where the interpolated (from adehabitatLT inttraj function) animal tracking datasets (in traj class) live:
datadir <- paste0(root,'datadir/')          

# Set Deposit Directory: This is where the UD results and any associated figures will live:
dropdir <- paste0(root,'ud_out/')    


###########################################################################################################################################################   
# PreAmble: Determining Grid Size. If you don't know if your spp is localized, intermediate, large, or vast ranged, use a grid size scaled to median max range
###########################################################################################################################################################   
# Note: The gridded UD analysis requires user-input of grid size. We used a fixed grid size for each home range size class (e.g. localized, intermediate, large, and vast). 
# For our initial run, we didn't know the home range size of our datasets, so we used adjustible grid sizes per dataset, informed by their median maximum range calculated per dataset 
# (See script XXXXX.r to calculate track characteristics, including median max range). 
# We assign a variable 'meta' which contains a datatable that includes median max range per dataset(s) for initial runs and home range size class after dataset is categorized for final run. 

meta<-read.csv(paste0(root,'metadir/meta_medMaxRange.csv')) # meta is a datatable that contains the columns: 1. 'sp_name'  and 2. 'median_max_range_km' and 3. 'sc' home range size class (will be NA until final iteration)

#######################################################################################################################################################
# Proceed with Gridded UDs 
####################################################################################################################################################### 


# Set Central Place Foragers vector -------------------------------
CPF_vec <- c("BFAL_018", "LAAL_002", "MABO_003") # Enter datasetID of central place foragers

# Set Global Parameters  ------------------------------------------
hr_ud       <- 90                         # Home Range UD designation (typically 95 or 90)
core_method <- "SP"                       # Core Method: Can be "50" for using the 50UD or "SP" for using the Seaman-Powell method
area_proj   <- "laea"                     # Specify projection for area calculations
weighted    <- "TRUE"                     # are grids weighted? (TRUE or FALSE) - if TRUE, grids will be weighted differently for central-place-foraging species and for nomadic species
grid_adj    <- "FALSE"                    # "TRUE" for an adjusted grid that relies on and "FALSE" for a fixed grid
meta        <- meta                       # metadata table with range information (relevant for adjustible grid) and home range size class 'sc' (relevant for fixed grid size depending on home range size class of dataset)


# Loop through datasets to construct gridded UDs --------------------

setwd(datadir)

a<- list.files(pattern="*.csv") # list of all interpolated datasets--

for (i in 1:length(a)) { # Loop through individual species/datasets--
  
  int_traj<-read.csv(a[i]) # read in interpolated datafile i
  sp_name <- int_traj$species_ds[1]
  
  spds_m <- AddWeightsFunction(sp_name, int_traj, CPF_vec) # Step 1: Add Weights for Gridded UD Analysis 
  
  spds_m$Lon<-Lon180to360(spds_m$Lon) 
  
  # Set Local Parameters:  --------------------------------------------------------------------------
  lon     <- spds_m$Lon                     # longitude positions in 0-360 degrees
  lat     <- spds_m$Lat                     
  weightv <- spds_m$Weight                  # weights column 
  
  # Run Home Range Function --------------------------------------------------------------------------
  
  ud_list <- HomeRangeFunctionPub(dropdir, sp_name, lon, lat, weightv, hr_ud, core_method, area_proj, weighted, grid_adj, meta)
  
  area_tab <- ud_list$area_tab 
  
  if (i==1) {
    ud_table <- area_tab
  }else{
    ud_table <- rbind(ud_table, area_tab)
  }
  
  save(ud_list,file=paste0(dropdir,sp_name, "_Poly_Grid_List.Rdata"))
  
  rm(area_tab, spds_m, sp_name, lat, lon, weightv, ud_list)   

  }


# Save datatable with speciesID, 90UD_area, core_thresh_UD, coreUD_area, samplesize 
require(dplyr)
tbl_df = mutate_if(ud_table[,c(1:7)], is.numeric, as.integer)
tbl_df <- tbl_df %>%
  mutate(latcenter = ud_table$latcenter,loncenter = ud_table$loncenter)

ndf<-as.data.frame(matrix(NA, nrow=length(a),ncol=2))
colnames(ndf) <- c("ds","n")
for (i in 1:length(a)) {
  m<-read.csv(a[i])
  ndf$ds[i]<-strsplit(a[i],"_")[[1]][1]
  ndf$n[i]<-length(unique(m$animID))
}

tbl_df2 <- tbl_df %>% mutate(n=ndf$n)
tbl_df2 <- tbl_df2 %>% dplyr::select(spp_ds, areaKM90, Seaman_Core_Thresh, areaKM_spCore, n)
write.csv(tbl_df2,file=paste0(dropdir,"df_homerange_areas_KM.csv"))







