HomeRangeFunctionPub <- function(dropdir, sp_name, lon, lat, weightv, hr_ud, core_method, area_proj, weighted, grid_adj, gridd, sc) {
  
  
  #------------------------------------------------------------------------------------------
  # Create Spatial Points from lon lat input
  #------------------------------------------------------------------------------------------
  
  # ------------------------------------
  mat <- as.data.frame(matrix(NA, nrow=length(lat),ncol=3))
  colnames(mat) <- c("lon","lat", "weight")
  mat$lon <- lon
  mat$lat <- lat
  mat$weight <- weightv
  coordinates(mat) = c("lon","lat")
  proj4string(mat) <- CRS("+proj=longlat +ellps=WGS84")
  
  # Find central lon and central lat to be used for Lambert Equal Area Azimuthal projection later in script:
  summ_sppi <- summary(mat@coords)
  lon_centr<-as.numeric(str_split(str_split(summ_sppi[4,],pattern=":")[[1]][2],pattern=" ")[[1]][1])
  lat_centr<-ifelse(is.na(as.numeric(str_split(str_split(summ_sppi[4,],pattern=":")[[2]][2],pattern=" ")[[1]][1])), as.numeric(str_split(str_split(summ_sppi[4,],pattern=":")[[2]][2],pattern=" ")[[1]][3]), as.numeric(str_split(str_split(summ_sppi[4,],pattern=":")[[2]][2],pattern=" ")[[1]][1]))
  
  #-----------------------------------------------------------------------------------------
  # Begin Gridded Utilization Distribution Analysis
  #------------------------------------------------------------------------------------------
  #-----------------------------------------------------------------------------------------
  # CREATE BASE GRID
  bb <- sp::bbox(mat)
  
  #-----------------------------------------------------------------------------------------
  # Establish Grid Size
  #------------------------------------------------------------------------------------------
  
  
  if (grid_adj == TRUE) { # Iniital run so use adjustible grid size scaled to median max range of dataset
    # 5% fixed grid scaling fctr
    fctr <- .05 
    maxRangei<-meta$median_max_range_km[which(meta$spp==as.character(sp_name))]
    gd<-vms::km2Degree(mean(mat$lon),abs(mean(mat$lat)),(maxRangei*fctr)) # grid -> 5% of spp medn max range
    
  }else{
    # Localized :     5 km2
    # Intermediate : 10 km2
    # Large:         25 km2
    # Vast:          50 km2
    sc<- as.character(meta$hr_sc[which(meta$sp_name==as.character(sp_name))]) # home range size class
    
    if (sc=="localized") {
      gd<-km2Degree(mean(mat$lon),mean(mat$lat),5)
    }else if (sc=="intermediate") {
      gd<-km2Degree(mean(mat$lon),mean(mat$lat),10)
    }else if (sc=="large") {
      gd<-km2Degree(mean(mat$lon),mean(mat$lat),25)
    }else if (sc=="vast") {
      gd<-km2Degree(mean(mat$lon),mean(mat$lat),50)
    } else {
      gd<-km2Degree(mean(mat$lon),mean(mat$lat),75) # LEST Pacific only
    } 
  }
  
  
  #-----------------------------------------------------------------------------------------
  # Proceed with Gridded UD
  #------------------------------------------------------------------------------------------
  
  cs <- c(gd,gd) # Grid Cell Size in Degrees
  
  # Establish spatial grid ---------------------------------------------------------------------------------
  cc <- bb[,1]+(cs/2)
  cd <- ceiling(diff(t(bb))/cs)
  sppi_grid <- GridTopology(cellcentre.offset = cc, cellsize = cs, cells.dim = cd) 
  sppi_sp <- SpatialGrid(sppi_grid) # grid of bounding box
  
  #  Project spatial grid to CRS (Coordinate Reference System) ----------------------------------
  proj4string(sppi_sp) <- CRS("+proj=longlat +ellps=WGS84") # Define Projection
  p4s <- CRS(proj4string(sppi_sp)) 
  sppi_sg <- SpatialGrid(sppi_grid, proj4str = p4s) 
  
  # Use below if want to see tracks and base grid
  # plot(sppi_sp) 
  # plot(mat, add=TRUE, pch=19, cex = .05, col="red")
  
  # --------------------------------------------------------------------------------------------------------
  # Convert lat lon animal tracks to spatial raster data. And, project
  # --------------------------------------------------------------------------------------------------------
  sppi_sppts <- SpatialPoints(mat)
  sppi_raster <- raster(sppi_sg) #create raster LAYER of animal locations - but has no values yet
  projection(sppi_raster) <- "+proj=longlat +ellps=WGS84" # project
  
  if (weighted=="TRUE") { 
    wt <- mat$weight  
    grid_ct_sppi<-rasterize(sppi_sppts, sppi_raster, wt, fun=sum, na.rm=TRUE)
  }else{
    grid_ct_sppi <- rasterize(sppi_sppts, sppi_raster,fun='count', na.rm=TRUE)
  }
  
  values(grid_ct_sppi) <- ifelse(is.na(values(grid_ct_sppi)), 0, values(grid_ct_sppi)) #turn NAs into 0
  # plot(grid_ct_sppi) # if you like. this shows the number of hits per cell (weighted).
  
  
  # --------------------------------------------------------------------------------------------------------
  # Create the Gridded Utilization Distribution
  # --------------------------------------------------------------------------------------------------------

    # Normalize:
    s <- sum(values(grid_ct_sppi), na.rm=TRUE)
    values(grid_ct_sppi) <- values(grid_ct_sppi)/s
    df <- data.frame(v=values(grid_ct_sppi), id=seq(1,ncell(grid_ct_sppi)))
    df.2 <- df[rev(order(df$v)),] 
    df.2$cs <- cumsum(df.2[,'v'])
    rownames(df.2)<-NULL
    df.2$cs[df.2$cs>0.9999999] <- NA # Replace 1s with NAs
    
    #CREATE BREAKS
    df.2$c<-cut(df.2$cs,breaks=seq(0,1,by=0.05),labels=seq(100,5,by=-5))
    df.3 <- df.2[order(df.2$id),]
    values(grid_ct_sppi) <- as.numeric(df.3$c)
    plot(grid_ct_sppi)
    
    # --------------------------------------------------------------------------------------------------------
    # Quantifying Gridded UD: Convert Raster to Polygon to Count Cells and Measure Area 
    # --------------------------------------------------------------------------------------------------------
    # Count Number of Cells
    # Project Cells to Geogrpahic Coordinate System
    
    if(min(getValues(grid_ct_sppi),na.rm=TRUE) == max(getValues(grid_ct_sppi),na.rm=TRUE)) {
      polyi.100 <-rasterToPolygons(grid_ct_sppi, fun=function(x){x>0},dissolve=TRUE)
      polyi.95  <- rasterToPolygons(grid_ct_sppi, fun=function(x){x>0},dissolve=TRUE)
      polyi.50  <-  rasterToPolygons(grid_ct_sppi, fun=function(x){x>0},dissolve=TRUE)
    }else{
      polyi.100 <- rasterToPolygons(grid_ct_sppi, fun=function(x){x>0},dissolve=TRUE) # dont use dissolve=TRUE if you want to count the number of cells!
      polyi.95 <-  rasterToPolygons(grid_ct_sppi, fun=function(x){x>0 & x<20},dissolve=TRUE)
      polyi.90 <-  rasterToPolygons(grid_ct_sppi, fun=function(x){x>0 & x<19},dissolve=TRUE)
      if(min(getValues(grid_ct_sppi),na.rm=TRUE) >= 11) {
        polyi.50 <-  rasterToPolygons(grid_ct_sppi, fun=function(x){x>0 & x<min(getValues(grid_ct_sppi),na.rm=TRUE)+1},dissolve=TRUE)
      }else{
        polyi.50 <-  rasterToPolygons(grid_ct_sppi, fun=function(x){x>0 & x<11},dissolve=TRUE)
      }  
    }
    
    SPoutList <- SPCoreFunction(grid_ct_sppi)
    polyi.SPcore <- SPoutList[[1]]
    ud_core_sppi <- SPoutList[[2]]
    
    #plot UD with contour lines around Core  -------------------------------------------------------------------
    qwikcont.core<-rasterToContour(grid_ct_sppi,levels=ud_core_sppi) #quick contour for simple plot
    
    png(filename=paste(dropdir,sp_name,"_","_fullUD",".png",sep=""),
        width= 3.25,
        height    = 3.25,
        units     = "in",
        res       = 1200,
        pointsize = 6)
    plot(grid_ct_sppi,col=alpha(rev(topo.colors(20)),.6))
    # plot(qwikcont95,add=TRUE,lwd=0.2,col="white")
    plot(qwikcont.core,add=TRUE,lwd=1,col="blue")
    plot(wrld_simpl, add=TRUE, col="grey")
    plot(elide(wrld_simpl, shift = c(360, 0)), add = TRUE, col="grey")
    dev.off()
    
    # --------- PROJECT TO LAMBERT AZIMUTHAL EQUAL AREA ----------------------------------------------------
    
    # Centering on midpoint of locations for each dataset 
    poly.p100<-spTransform(polyi.100,  paste("+proj=laea +ellps=WGS84 +lon_0=",lon_centr,"+lat_0=",lat_centr ,"  +x_0=0 +y_0=0 ",sep=""))
    poly.p95<-spTransform(polyi.95,  paste("+proj=laea +ellps=WGS84 +lon_0=",lon_centr,"+lat_0=",lat_centr ,"  +x_0=0 +y_0=0 ",sep=""))
    poly.p90<-spTransform(polyi.90,  paste("+proj=laea +ellps=WGS84 +lon_0=",lon_centr,"+lat_0=",lat_centr ,"  +x_0=0 +y_0=0 ",sep=""))
    poly.p50<-spTransform(polyi.50,  paste("+proj=laea +ellps=WGS84 +lon_0=",lon_centr,"+lat_0=",lat_centr ,"  +x_0=0 +y_0=0 ",sep=""))
    polyp.core<-spTransform(polyi.SPcore,  paste("+proj=laea +ellps=WGS84 +lon_0=",lon_centr,"+lat_0=",lat_centr ,"  +x_0=0 +y_0=0 ",sep=""))
    
    # The coordinate reference system is not longitude/latitude. Use rgeos::gArea instead
    aSP<-round(sum(rgeos::gArea(polyp.core))/1e6,2) #km^2
    a50<-round(sum(rgeos::gArea(poly.p50))/1e6,2) #km^2
    a90<-round(sum(rgeos::gArea(poly.p90))/1e6,2) #km^2
    a95<-round(sum(rgeos::gArea(poly.p95))/1e6,2) #km^2
    a100<-round(sum(rgeos::gArea(poly.p100))/1e6,2) #km^2
    
    #add to table
    tb<-as.data.frame(matrix(NA,nrow=1,ncol=9))
    colnames(tb)<-c("spp_ds","areaKM_100","areaKM_95","areaKM_90","areaKM_50","Seaman_CoreThresh","areaKM_spCore", "loncenter", "latcenter")
    tb$spp_ds<-sp_name
    tb$areaKM_100<- a100
    tb$areaKM_95<-a95
    tb$areaKM_90<-a90
    tb$areaKM_50<-a50
    tb$Seaman_CoreThresh <- ud_core_sppi*5
    tb$areaKM_spCore<-aSP
    tb$loncenter <- lon_centr
    tb$latcenter <- lat_centr
    
    # Make List of Objects to Return
    ud_out_list<-list("area_tab" = tb,"sppi_raster" = grid_ct_sppi, "HomeRangeGrid" = polyi.90 , "HomeRangePoly" = poly.p90, "CoreGrid" = polyi.SPcore, "CorePoly" =polyp.core)
    # ud_out_list<-list("sppi_raster" = grid_ct_sppi)
    
    return(ud_out_list)
    
  
  
}
