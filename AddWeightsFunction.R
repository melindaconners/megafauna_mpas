# ------------------------------------------------------------------------------------------------------------------
# This Script Adds a Column of Weights for central-place foragers and non central-place foragers to be used in the gridded UD Analysis. 
# For CPFs, weights are calculated as the inverse number of trips. 
# For non-CPFs, weights are calculated with trip duration (with the 85% threshold (from Block et al.2011) )
# 
# Written by Melinda Conners, melinda.conners@stonybrook.edu for Conners et al XXXX, PNAS
# Last modified January 10, 2020
# ------------------------------------------------------------------------------------------------------------------

####################################################################################################################
# Input 
####################################################################################################################

# 1. sp_name <- species_name or dataset_id character vector for dataset i (e.g. "BFAL_018")
# 2. int_traj<- csv file for dataset i which is the interpolated datafile (output from the redisltraj function from adehabitatLT package 
#               with three additional columns: 
#               1. 'species_ds' (character: species/dataset ID) '
#               2. 'animID' (character: individual animal ID) 
#               3. 'tripID' (character: tripid (trip number) of animal specified in animID)) 
#               *** all columns in int-traj: x,y,date,dx,dy,dist,dt,R2n,abs.angle,rel.angle,id,burst,pkey,species_ds,animID,tripID
# 3. CPF_vec<- a character vector of species_name or dataset_id that will be used to check if dataset is a CPF or not. This will determine the weighting scheme used. e.g. CPF_vec <- c("BFAL_018", "LAAL_002", "MABO_003") or CPF_vec <- "NA"


####################################################################################################################
# Output 
####################################################################################################################

# A n x 7 matrix ('newmat') with the following columns:
# 1. Species  # species, or, dataset identifier
# 2. Time     # POSIXct (format:                                                                                      ??????????
# 3. Lat      # -90 to 90
# 4. Lon      # 360                                                                                                   ????????
# 5. Weight   # weight to be used in gridded UD analysis
# 6. tripID   # a unique identifier of which trip, for animals that have multiple trips (e.g. many CPFs)
# 7. animID   # unique animal identifier

AddWeightsFunction<-function(sp_name, int_traj, CPF_vec) {

  # 'int_traj' is the interpolated dataframe
  # 'CPF_vec' is a character vector of sp_name(s) names which specifies if the species is a central-place forager or not. 
  
  # Set traj as dataframe for loop
     m_df<-as.data.frame(int_traj)
     m_df$weight <- NA
   
  # Weights for Central Place Foragers:  -----------------------------------------------------------------------------------
    if (length(which(sp_name==CPF_vec))>=1) { # If this dataset is of a CPF, do the following:
      
      nbirds<-unique(m_df$animID) # Loop through each individual animal to count its trips
      
      for (k in 1:length(nbirds)) {  
        sub_anim<-m_df[m_df$animID==nbirds[k],]
        ntrips<-length(unique(sub_anim$tripID))
        m_df$weight[m_df$animID==nbirds[k]]<-1/ntrips
      }
      rm(k)
      
      wt <- m_df$weight # weight vector
     
     
      }else{ # If not a CPF -----------------------------------------------------------------------------------------------
        
        # Set weights based on species-level information of trip durations 
        trips<-unique(m_df$tripID)
        trip_lengths<-as.data.frame(matrix(NA,ncol=2,nrow=length(trips)))
        colnames(trip_lengths) <- c("tripID", "V1")
        
        for (k in 1:length(trip_lengths$V1)) {
          tmp<-m_df[m_df$tripID==trips[k],]
          trip_lengths$tripID[k] <- as.character(trips[k])
          trip_lengths$V1[k] <- length(tmp$x)
        }
        
        raz<-matrix(NA,nrow=max(trip_lengths$V1),ncol=length(trip_lengths$tripID)) # This creates an empty matrix that will be filled with "ones" that represent trip duration
        # Each column in raz is a different trip
        # Each row in raz will become either a 1 or an NA - if there is a 1 it means there is a location at that index
        # Then each row will be summed to give the count of trips at that specific 'trip length' index. 
 
        # Create the value matrix (raz) ---------------------------------------------------------------------
        for (j in 1:length(trip_lengths$tripID)){
             raz[1:trip_lengths$V1[j],j]<-1
        }
        rm(j)

        # Sum rows to get number of individuals at each location index --------------------------------------
        countvec<-apply(raz,1,sum, na.rm=TRUE) # countvec should be the same length as the largest trip duration
        weights<-1/countvec
      
        # 85th-%-ile of  track-lengths for species dataset "i":
        t_85val<- quantile(trip_lengths$x, .85)
    
       # Loop through trips and add weights ---------------------------------------------------------------------------------
       trips<-unique(m_df$id)
      
       for (j in 1:length(trips)) { # --------------------------------------------------------------------------------------
       
           trip_m <- m_df[which(m_df$id==trips[j]),]
           trip_m$index <- seq(1,length(trip_m$y),by=1)
      
           # add vec with threshed weights of length of current trip. 
           threshed_weights=weights[1:length(trip_m$x)] 
           tix<- which(trip_m$index==round(as.numeric(t_85val),0))
            
              if (pracma::isempty(tix)) {
                threshed_weights<-threshed_weights
              }else{
                threshed_weights[tix:length(trip_m$x)]<-threshed_weights[tix]
              }
      
          threshed_weights<-as.data.frame(threshed_weights)
     
          if (j==1) {
            threshed_weights_full <- as.data.frame(threshed_weights)
          }else{
            threshed_weights_full <- rbind(threshed_weights_full, threshed_weights)
          }
     
          rm(threshed_weights)
      
         }
         rm(j) 
 
         wt<-matrix(unlist(threshed_weights_full))  # weight vector
        
      } # End non CPF Part of Else Loop >> end result is a 'wt' vector ----------------------------------------------------------------
        
     
     
     
     
      # set new mat that only contains relevant information from the interpolated mat. 
      newmat<-matrix(NA, nrow=length(m_df$x),ncol=7)
      newmat<-as.data.frame(newmat)
      colnames(newmat) <- c('Species',"Time","Lat","Lon","Weight","tripID","animID")
  
      newmat$Weight <- as.numeric(wt)
      newmat$Species<- sp_name 
      newmat$Time   <- m_df$date
      newmat$Lat    <- m_df$y
      newmat$Lon    <- m_df$x
      newmat$animID <- m_df$animID
      newmat$tripID <- m_df$tripID

      
      return(newmat)
      
        
     } # End
     
     
     
     