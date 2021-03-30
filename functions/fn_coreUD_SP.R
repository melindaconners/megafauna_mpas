

SPCoreFunction<- function(grid_ct_sppi) {
    
    # --------------------------------------------------------------------------------------------------------
    # Quantifying Gridded UD: Using Seaman-Powell Method for Core Definition:
    # --------------------------------------------------------------------------------------------------------
   
    #########################################################
    # ------  Determine the Core Area      ------- #
    #########################################################
   # This step calculates relative area based on number of cells (not actual area)
 
    grid1<-(grid_ct_sppi==1) # 5
    area1<-(sum(getValues(grid1),na.rm=TRUE)*res(grid1)[1]*res(grid1)[2])
    grid2<-(grid_ct_sppi==2) # 10
    area2<-(sum(getValues(grid2),na.rm=TRUE)*res(grid2)[1]*res(grid2)[2])
    grid3<-(grid_ct_sppi==3) # 15
    area3<-(sum(getValues(grid3),na.rm=TRUE)*res(grid3)[1]*res(grid3)[2])
    grid4<-(grid_ct_sppi==4) # 20
    area4<-(sum(getValues(grid4),na.rm=TRUE)*res(grid4)[1]*res(grid4)[2])
    grid5<-(grid_ct_sppi==5) # 25
    area5<-(sum(getValues(grid5),na.rm=TRUE)*res(grid5)[1]*res(grid5)[2])
    grid6<-(grid_ct_sppi==6) # 30
    area6<-(sum(getValues(grid6),na.rm=TRUE)*res(grid6)[1]*res(grid6)[2])
    grid7<-(grid_ct_sppi==7) # 35
    area7<-(sum(getValues(grid7),na.rm=TRUE)*res(grid7)[1]*res(grid7)[2])
    grid8<-(grid_ct_sppi==8) # 40
    area8<-(sum(getValues(grid8),na.rm=TRUE)*res(grid8)[1]*res(grid8)[2])
    grid9<-(grid_ct_sppi==9) # 45
    area9<-(sum(getValues(grid9),na.rm=TRUE)*res(grid9)[1]*res(grid9)[2])
    grid10 <-(grid_ct_sppi==10) # 50
    area10 <-(sum(getValues(grid10), na.rm=TRUE) * res(grid10)[1]*res(grid10)[2])
    grid11<-(grid_ct_sppi==11) # 55
    area11<-(sum(getValues(grid11),na.rm=TRUE)*res(grid11)[1]*res(grid11)[2])
    grid12<-(grid_ct_sppi==12) # 60
    area12<-(sum(getValues(grid12),na.rm=TRUE)*res(grid12)[1]*res(grid12)[2])
    grid13<-(grid_ct_sppi==13) # 65
    area13<-(sum(getValues(grid13),na.rm=TRUE)*res(grid13)[1]*res(grid13)[2])
    grid14<-(grid_ct_sppi==14) # 70
    area14<-(sum(getValues(grid14),na.rm=TRUE)*res(grid14)[1]*res(grid14)[2])
    grid15<-(grid_ct_sppi==15) # 75
    area15<-(sum(getValues(grid15),na.rm=TRUE)*res(grid15)[1]*res(grid15)[2])
    grid16<-(grid_ct_sppi==16) # 80
    area16<-(sum(getValues(grid16),na.rm=TRUE)*res(grid16)[1]*res(grid16)[2])
    grid17<-(grid_ct_sppi==17) # 85
    area17<-(sum(getValues(grid17),na.rm=TRUE)*res(grid17)[1]*res(grid17)[2])
    grid18<-(grid_ct_sppi==18) # 90
    area18<-(sum(getValues(grid18),na.rm=TRUE)*res(grid18)[1]*res(grid18)[2])
    grid19<-(grid_ct_sppi==19) # 95
    area19<-(sum(getValues(grid19),na.rm=TRUE)*res(grid19)[1]*res(grid19)[2])
    grid20 <- (grid_ct_sppi==20) # 100
    area20 <- (sum(getValues(grid20), na.rm=TRUE) * res(grid20)[1]*res(grid20)[2])

    ud_area <- c(area1, area2, area3, area4, area5, area6, area7, area8, area9, area10,area11, area12, area13, area14, area15, area16, area17, area18, area19, area20)
    
    cpclass<-c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
    cpclass<-sort(cpclass,decreasing=TRUE)
    area_sppi<-c(0,ud_area)
    area_sppi<-sort(area_sppi,decreasing=TRUE)
    prop_area<-area_sppi/max(area_sppi) #PCTPROB
    
    area_line<-c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)
    area_line<-sort(area_line,decreasing=TRUE)
    
    # Calculate area of each isopleth
    # This is a rough calculation of area (it's a RELATIVE calculation.)
    # It's basically just counting grids of each isopleth and then multipying it by the area of the grid cells. 
  
    # par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0) #default
    # png(filename=paste(root,"02 analysis/gridded UD/10262018/core_thresh_plots/",sp_name,"_PowellCore_LinePlot",".png",sep=""),
    #     width= 4,
    #     height    = 3.25,
    #     units     = "in",
    #     res       = 1200,
    #     pointsize = 6)
    # plot(cpclass,prop_area,type="o",pch=19,xlab="Probability of Use (UD contour)",ylab="Proportion of Area")
    # lines(cpclass,area_line,lwd=2)
    
    #########################################################
    # Calculate distance, find the maximum distance, and plot
    #########################################################
    
    dist_powell<-abs(area_line-prop_area) 
    core_thresh<-max(dist_powell)
    foo <- which(dist_powell==core_thresh)
    if (length(foo)>1) {
     foo<-foo[which(foo==max(foo))] 
    }
    ud_core_sppi <- (cpclass[foo]/5)-1
    abline(v=ud_core_sppi*5,col="red")

    dev.off()
    
    if(ud_core_sppi<min(getValues(grid_ct_sppi), na.rm=TRUE)) {
      ud_core_sppi <- min(getValues(grid_ct_sppi), na.rm=TRUE)
    }
    
    # Create a raster with just the core area classes
    core_raster <- grid_ct_sppi
    core_raster_dissolve <- grid_ct_sppi
    values(core_raster) <- ifelse(values(grid_ct_sppi) <= ud_core_sppi, values(grid_ct_sppi), NA)
    values(core_raster_dissolve) <- ifelse(values(grid_ct_sppi) <= ud_core_sppi, 1, NA)
  
    # Powell-Seaman Core
    # Convert Raster Values to Polygon Grid and Contour Lines
    polyi.SPcore <-  rasterToPolygons(core_raster_dissolve, fun=function(x){x>0 & x<ud_core_sppi},dissolve=TRUE)
    
    # plot(polyi.SPcore, col=alpha("grey", .5))
    # plot(mat, add=TRUE, cex = .1, pch=16)
    
    return(list(polyi.SPcore, ud_core_sppi))
}




