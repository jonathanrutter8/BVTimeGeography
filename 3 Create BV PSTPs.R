########################################################################################
################  PAIRWISE VOXEL-BASED PSTPs FOR BIRDS AND VESSELS  ####################
########################################################################################

# Takes in CSVs of bird and vessel data
# Start by ensuring you have a folder called “PSTP Outputs FINAL” in your WD with the following empty folders:
  # "B Subset Dataframes"
  # "BV PSTP Bricks"
  # "Interaction Bricks"
  # "Summary Tables"
# Next, adjust section B and function arguments as required
# Next, double check input column names, as these could impact creating regular trajectories
  # See Bdf4 and Vdf4 objects within section C
# This script is set up to run in parallel

####################### CONTENTS ###########################

# A. LOAD PACKAGES
# B. SET GLOBAL DATA AND PARAMETERS
# C. DEFINE FUNCTION TO GET B AND V PSTPS
# 1. Subset B and V trajectories
# 2. Set coordinate reference system
# 3. Set temporal extent and time steps
# 4. Create regular trajectory for bird (B)
# 5. Create regular trajectory for vessel (V)
# 6. Set spatial extent and define space-time cube
# 7. Construct bird PSTP
# 8. Construct vessel PSTP
# 9. Create PSTP Summary Table
# D. GET B AND V PSTPS USING FUNCTION



############# A. LOAD PACKAGES  ###################

library(plyr)
library(tidyverse)
library(lubridate)
library(sp)
library(move)
library(raster) # May need 2022 version of this package. If facing errors, uninstall and reinstall using package "versions"
library(svMisc)
library(ncdf4)
library(parallel)
library(pbapply)




############# B. SET GLOBAL DATA AND PARAMETERS ################

# Set WD
WD <- "_____________________"
setwd(WD)

# Read in CSVs
Bdata <- read.csv("Rainier Clean Clipped/Rainier_Clip_Dissertation.csv")
Vdata <- read.csv("GFW Queried Data/GFW_5daytest_clean.csv")

# Set max velocities (B is a constant of units m/s, V is a unitless multiplier.)
B_max_velocity <- 30 # Based on Merkel et al. 2016
Vmultiplier <- 1.5

# Set telemetry error
sigma2_2 <- 12 # 12m error for Rainier GPS fixes

# Set CRS of input data: WGS 1984
WGS1984 <- CRS("+init=epsg:4326")

# Prep for parallel processing
numCores <- detectCores()
cl <- makeCluster(numCores)
clusterEvalQ(cl, {
  
  library(plyr)
  library(tidyverse)
  library(lubridate)
  library(sp)
  library(move)
  library(raster)
  library(svMisc)
  library(ncdf4)
  library(parallel)
  library(pbapply)
  
  # Set WD
  WD <- "C:/Users/jonat/Desktop/University of Oxford BCM/Oxford BCM 2021-2022/BCM Dissertation/R Scripts and Data/"
  setwd(WD)
  
  # Read in CSVs
  Bdata <- read.csv("Rainier Clean Clipped/Rainier_Clip_Dissertation.csv")
  Vdata <- read.csv("GFW Queried Data/GFW_5daytest_clean.csv")
  
  # Set max velocities (B is a constant of units m/s, V is a unitless multiplier.)
  B_max_velocity <- 30 # Based on Merkel et al. 2016
  Vmultiplier <- 1.5
  
  # Set telemetry error
  sigma2_2 <- 12 # 12m error for Rainier GPS fixes
  
  # Set CRS of input data: WGS 1984
  WGS1984 <- CRS("+init=epsg:4326")
})







############### C. DEFINE FUNCTION TO GET B AND V PSTPS ######################

# This function creates 1 B and 1 V PSTP and saves both as raster bricks in the "BV PSTP Bricks" folder.
# It also creates a B subset dataframe and saves it in the "B Subset Dataframes" folder.
# Function returns the elapsed time.

# Parameters:
# SubTrajNum_V: the vessel trajectory subset (can loop over this if needed)
# dt: voxel temporal resolution in seconds (5 mins = 300 secs, 30 mins = 1800 secs)
# dxy: voxel spatial resolution in metres

# START OF FUNCTION
get.BV.PSTP <- function(SubTrajNum_V, dt, dxy){
  
  ############ 1. SUBSET B AND V TRAJECTORIES #############
  
  # Start the clock
  start_time <- Sys.time()
  
  # Subset vessel dataframe
  Vdf1 <- Vdata %>%
    mutate(birdID = as.factor(birdID), DateTime = ymd_hms(DateTime)) %>%
    filter(SubTraj == SubTrajNum_V)
  
  # Determine first and last DateTimes, birdID, vesselID of vessel dataframe
  first_DT_Vdf1 <- min(Vdf1$DateTime)
  last_DT_Vdf1 <- max(Vdf1$DateTime)
  birdID_Vdf1 <- Vdf1$birdID[1]
  vesselID_Vdf1 <- Vdf1$vesselID[1]
  
  # Determine indices of Bdata that are nearest to first and last DateTimes of Vdf1 (same bird)
  Bdata1 <- Bdata %>% filter(ID == birdID_Vdf1)
  first_index_Bdata1 <- which(abs(ymd_hms(Bdata1$DateTime) - first_DT_Vdf1) == 
                                min(abs(ymd_hms(Bdata1$DateTime) - first_DT_Vdf1)))
  last_index_Bdata1 <- which(abs(ymd_hms(Bdata1$DateTime) - last_DT_Vdf1) == 
                               min(abs(ymd_hms(Bdata1$DateTime) - last_DT_Vdf1)))
  
  # Subset to corresponding B df
  Bdf1 <- Bdata1[first_index_Bdata1:last_index_Bdata1,] %>%
    mutate(ID = as.factor(ID), DateTime = ymd_hms(DateTime))
  
  # Return early (abort the function) if number of B positions is less than 2
  if(nrow(Bdf1) < 2){return(0)}
  
  
  
  
  
  ################ 2. SET COORDINATE REFERENCE SYSTEM ###############################
  
  # Turn trajectories into move objects
  Bdf2 <- move(x = Bdf1$Long,
               y = Bdf1$Lat,
               time = Bdf1$DateTime,
               proj = WGS1984,
               animal = Bdf1$ID,
               data = Bdf1,
               removeDuplicatedTimestamps = TRUE)
  Vdf2 <- move(x = Vdf1$Long,
               y = Vdf1$Lat,
               time = Vdf1$DateTime,
               proj = WGS1984,
               animal = Vdf1$ssvid,
               data = Vdf1,
               removeDuplicatedTimestamps = TRUE)
  
  # Transform coordinates to Azimuthal Equi-distance projection centered at at center of bird track (default)
  Bdf2 <- spTransform(Bdf2, center = TRUE)
  Vdf2 <- spTransform(Vdf2, CRSobj = Bdf2@proj4string)
  
  
  
  
  
  
  
  ############## 3. SET TEMPORAL EXTENT AND TIME STEPS ###################################
  
  # Figure out temporal extent, with whole number of time steps
  # Did not apply a temporal buffer because GFW query was already buffered
  Tmin <- max(c(min(Bdf1$DateTime), min(Vdf1$DateTime))); Tmin # Don't round
  Tmax_init <- min(c(max(Bdf1$DateTime), max(Vdf1$DateTime))); Tmax_init
  Trange_init <- abs(as.numeric(difftime(Tmax_init, Tmin, units = "secs"))); Trange_init
  Trange <- round_any(Trange_init, dt, f = ceiling); Trange # Make Trange a multiple of dt
  Tmax <- Tmin + Trange; Tmax
  
  
  # Determine number of time steps
  n_Tsteps <- (Trange / dt)
  
  # Get vector of time steps (centroids)
  vector_Tsteps <- Tmin
  for(i in 1:(n_Tsteps - 1)){
    vector_Tsteps <- c(vector_Tsteps, (Tmin + (i * dt)))
  }
  
  
  
  
  ############### 4. CREATE REGULAR TRAJECTORY FOR BIRD (B) #############
  
  # Interpolate to times listed in vector_Tsteps
  Bdf3 <- interpolateTime(Bdf2,
                          time = vector_Tsteps,
                          spaceMethod = "euclidean")
  
  # Transform into data.frame and return
  Bdf4 <- as.data.frame(Bdf3)
  
  # Subset to columns: coords.x1, coords.x2, timestamps, ID, Meta, locQual, 
  # dt_int, dist_int, speed, birdhour, sigma1_2, sensor
  Bdf4 <- Bdf4[,c("coords.x1", "coords.x2", "timestamps", "ID", "Meta", "LocQual", 
                  "dt_int", "dist_int", "speed", "birdhour", "sigma1_2", "sensor")]
  
  # If sensor = "interpolateTime", change to "Interpolated", change locQual to "I"
  # Then ensure data types are correct, change names of columns
  Bdf4$LocQual <- as.character(Bdf4$LocQual)
  Bdf4$sensor <- as.character(Bdf4$sensor)
  for(i in 1:nrow(Bdf4)){
    if(Bdf4[i,"sensor"] == "interpolateTime"){
      Bdf4[i,"LocQual"] <- "I"
      Bdf4[i,"sensor"] <- "Interpolated"
    } else {
      Bdf4[i,"sensor"] <- "Real"
    }
  }
  Bdf4$LocQual <- as.factor(Bdf4$LocQual)
  Bdf4$sensor <- as.factor(Bdf4$sensor)
  Bdf4$timestamps <- ymd_hms(Bdf4$timestamps)
  Bdf4 <- Bdf4 %>% rename(Xcoords = coords.x1, Ycoords = coords.x2, DateTime = timestamps, Interpolated = sensor)
  
  # Add in missing real positions: make data.frame of reals that resembles interpolated
  Bdf5 <- as.data.frame(Bdf2)
  Bdf5 <- Bdf5[,c("coords.x1", "coords.x2", "timestamps", "ID", "Meta", "LocQual", 
                  "dt_int", "dist_int", "speed", "birdhour", "sigma1_2")]
  Bdf5 <- Bdf5 %>% 
    mutate(Interpolated = "Real") %>%
    rename(Xcoords = coords.x1, Ycoords = coords.x2, DateTime = timestamps)
  
  # Rbind
  Bdf6 <- rbind(Bdf4, Bdf5)
  
  # Reorder by DateTime, remove duplicate timestamps, fill in speed and sigma1_2 for interpolated positions
  Bdf6 <- Bdf6[order(Bdf6[,"DateTime"]),]
  Bdf6 <- Bdf6 %>% 
    fill(speed, .direction = "up") %>% 
    fill(sigma1_2, .direction = "down") %>% # Fill in variances for upcoming interpolated positions
    fill(sigma1_2, .direction = "up") # Then fill those that were missed at start of trajectory
  Bdf6 <- Bdf6[!duplicated(Bdf6$DateTime),]
  
  
  # If an Interpolated position is within timestep/2 of a Real position, change to Real
  Bdf7 <- Bdf6
  for(i in 1:nrow(Bdf7)){
    if(i == 1){
      if((Bdf7[i, "Interpolated"] == "Interpolated") && 
         (Bdf7[i+1, "Interpolated"] == "Real") && 
         (abs(as.numeric(difftime(Bdf7[i, "DateTime"], Bdf7[i+1, "DateTime"], units = "secs"))) <= (dt/2))){
        Bdf7[i, "Interpolated"] <- "Real"
      }
    } else if(i == nrow(Bdf7)) {
      if((Bdf7[i, "Interpolated"] == "Interpolated") && 
         (Bdf7[i-1, "Interpolated"] == "Real") && 
         (abs(as.numeric(difftime(Bdf7[i, "DateTime"], Bdf7[i-1, "DateTime"], units = "secs"))) <= (dt/2))){
        Bdf7[i, "Interpolated"] <- "Real"
      }
    } else {
      if((Bdf7[i, "Interpolated"] == "Interpolated") && 
         (((Bdf7[i-1, "Interpolated"] == "Real") && 
           (abs(as.numeric(difftime(Bdf7[i, "DateTime"], Bdf7[i-1, "DateTime"], units = "secs"))) <= (dt/2))) || 
          ((Bdf7[i+1, "Interpolated"] == "Real") && 
           (abs(as.numeric(difftime(Bdf7[i, "DateTime"], Bdf7[i+1, "DateTime"], units = "secs"))) <= (dt/2))))){
        Bdf7[i, "Interpolated"] <- "Real"
      }
    }
  }
  
  # Remove Real positions that do not fall directly on a voxel centroid
  Bdf8 <- Bdf7[Bdf7$DateTime %in% vector_Tsteps,]
  
  # Make first and last position on the trajectory Real
  # In some cases may introduce a small amount of error, but this is critical for STP to work
  Bdf8[1, "Interpolated"] <- "Real"
  Bdf8[nrow(Bdf8), "Interpolated"] <- "Real"
  
  # Initialise Index1, t_since (time since previous Real), t_until (time until next Real),
  # Initialise alpha and variance (see Buchin et al 2012)
  Bdf8 <- Bdf8 %>%
    mutate(Index1 = row_number(), t_since = 0, t_until = 0, alpha = 0, variance = 0)
  
  # Calculate t_since and t_until: For each time step...
  for(i in 1:nrow(Bdf8)){
    if(Bdf8[i, "Interpolated"] == "Real"){
      
      # Calculate t_since (units in seconds)
      if(i < nrow(Bdf8)){
        for(j in 1:(nrow(Bdf8) - i)){
          if(Bdf8[i+j, "Interpolated"] == "Interpolated"){
            Bdf8[i+j, "t_since"] <- abs(as.numeric(difftime(Bdf8[i, "DateTime"], Bdf8[i+j, "DateTime"], units = "secs")))
          } else {
            break
          }
        }
      }
      
      # Calculate t_until (units in seconds)
      if(i > 1){
        for(k in 1:(i - 1)){
          if(Bdf8[i-k, "Interpolated"] == "Interpolated"){
            Bdf8[i-k, "t_until"] <- abs(as.numeric(difftime(Bdf8[i, "DateTime"], Bdf8[i-k, "DateTime"], units = "secs")))
          } else {
            break
          }
        }
      }
    }
  }
  
  # Calculate alpha and variance: For each time step...
  for(i in 1:nrow(Bdf8)){
    
    # Calculate alpha (value from 0-1 as agent moves from Real1 to Real2)
    Bdf8[i, "alpha"] <- ifelse(Bdf8[i, "Interpolated"] == "Real", 
                               0, 
                               Bdf8[i, "t_since"] / (Bdf8[i, "t_since"] + Bdf8[i, "t_until"]))
    
    # Calculate variance (increases towards midpoint between Real1 and Real2)
    Bdf8[i, "variance"] <- ((Bdf8[i, "t_since"] + Bdf8[i, "t_until"]) * Bdf8[i, "alpha"] * (1 - Bdf8[i, "alpha"]) * Bdf8[i, "sigma1_2"]) +
      ((Bdf8[i, "alpha"]^2 + (1 - Bdf8[i, "alpha"])^2) * sigma2_2)
  }
  
  # Write CSV of B subset regular trajectory
  setwd(paste(WD, "/PSTP Outputs TEST THIN VOXELS/B Subset Dataframes/", sep = "")) # CHANGE HERE (Outputs FINAL becomes Outputs TEST THIN VOXELS)
  write.csv(Bdf8, paste("B_RegTraj_vox", dxy/1000, "_SubTraj", SubTrajNum_V, ".csv", sep = ""))
  setwd(WD)
  
  
  
  
  
  
  
  ########### 5. CREATE REGULAR TRAJECTORY FOR VESSEL (V) #################
  
  # Similar to #4 above, but with Vdf instead of Bdf, and no dBBMM variance calculation
  # Note that column names are different
  
  # Interpolate to times listed in vector_Tsteps
  Vdf3 <- interpolateTime(Vdf2,
                          time = vector_Tsteps,
                          spaceMethod = "euclidean")
  
  # Transform into data.frame and return
  Vdf4 <- as.data.frame(Vdf3)
  
  # Subset to columns: coords.x1, coords.x2, timestamps, ID, Meta, locQual, 
  # dt_int, dist_int, speed, sensor
  Vdf4 <- Vdf4[,c("coords.x1", "coords.x2", "timestamps", "flag", "ssvid", "fishing", 
                  "dt_int", "dist_int", "speed", "sensor")]
  
  # If sensor = "interpolateTime", change to "Interpolated"
  # Then ensure data types are correct, change names of columns
  Vdf4$sensor <- as.character(Vdf4$sensor)
  for(i in 1:nrow(Vdf4)){
    if(Vdf4[i,"sensor"] == "interpolateTime"){
      Vdf4[i,"sensor"] <- "Interpolated"
    } else {
      Vdf4[i,"sensor"] <- "Real"
    }
  }
  Vdf4$sensor <- as.factor(Vdf4$sensor)
  Vdf4$timestamps <- ymd_hms(Vdf4$timestamps)
  Vdf4 <- Vdf4 %>% rename(Xcoords = coords.x1, Ycoords = coords.x2, DateTime = timestamps, Interpolated = sensor)
  
  # Add in missing real positions: make data.frame of reals that resembles interpolated
  Vdf5 <- as.data.frame(Vdf2)
  Vdf5 <- Vdf5[,c("coords.x1", "coords.x2", "timestamps", "flag", "ssvid", "fishing", 
                  "dt_int", "dist_int", "speed")]
  Vdf5 <- Vdf5 %>% 
    mutate(Interpolated = "Real") %>%
    rename(Xcoords = coords.x1, Ycoords = coords.x2, DateTime = timestamps)
  
  # Rbind
  Vdf6 <- rbind(Vdf4, Vdf5)
  
  # Reorder by DateTime, remove duplicate timestamps, fill in speed for interpolated positions
  Vdf6 <- Vdf6[order(Vdf6[,"DateTime"]),]
  Vdf6 <- Vdf6 %>% fill(speed, .direction = "up")
  Vdf6 <- Vdf6[!duplicated(Vdf6$DateTime),]
  
  
  # If an Interpolated position is within timestep/2 of a Real position, change to Real
  Vdf7 <- Vdf6
  for(i in 1:nrow(Vdf7)){
    if(i == 1){
      if((Vdf7[i, "Interpolated"] == "Interpolated") && 
         (Vdf7[i+1, "Interpolated"] == "Real") && 
         (abs(as.numeric(difftime(Vdf7[i, "DateTime"], Vdf7[i+1, "DateTime"], units = "secs"))) <= (dt/2))){
        Vdf7[i, "Interpolated"] <- "Real"
      }
    } else if(i == nrow(Vdf7)) {
      if((Vdf7[i, "Interpolated"] == "Interpolated") && 
         (Vdf7[i-1, "Interpolated"] == "Real") && 
         (abs(as.numeric(difftime(Vdf7[i, "DateTime"], Vdf7[i-1, "DateTime"], units = "secs"))) <= (dt/2))){
        Vdf7[i, "Interpolated"] <- "Real"
      }
    } else {
      if((Vdf7[i, "Interpolated"] == "Interpolated") && 
         (((Vdf7[i-1, "Interpolated"] == "Real") && 
           (abs(as.numeric(difftime(Vdf7[i, "DateTime"], Vdf7[i-1, "DateTime"], units = "secs"))) <= (dt/2))) || 
          ((Vdf7[i+1, "Interpolated"] == "Real") && 
           (abs(as.numeric(difftime(Vdf7[i, "DateTime"], Vdf7[i+1, "DateTime"], units = "secs"))) <= (dt/2))))){
        Vdf7[i, "Interpolated"] <- "Real"
      }
    }
  }
  
  # Remove Real positions that do not fall directly on a voxel centroid
  Vdf8 <- Vdf7[Vdf7$DateTime %in% vector_Tsteps,]
  
  # Make first and last position on the trajectory Real
  # In some cases may introduce a small amount of error, but this is critical for STP to work
  Vdf8[1, "Interpolated"] <- "Real"
  Vdf8[nrow(Vdf8), "Interpolated"] <- "Real"
  
  # Initialise Index1, t_since (time since previous Real), t_until (time until next Real),
  Vdf8 <- Vdf8 %>%
    mutate(Index1 = row_number(), t_since = 0, t_until = 0)
  
  # Calculate t_since and t_until: For each time step...
  for(i in 1:nrow(Vdf8)){
    if(Vdf8[i, "Interpolated"] == "Real"){
      
      # Calculate t_since (units in seconds)
      if(i < nrow(Vdf8)){
        for(j in 1:(nrow(Vdf8) - i)){
          if(Vdf8[i+j, "Interpolated"] == "Interpolated"){
            Vdf8[i+j, "t_since"] <- abs(as.numeric(difftime(Vdf8[i, "DateTime"], Vdf8[i+j, "DateTime"], units = "secs")))
          } else {
            break
          }
        }
      }
      
      # Calculate t_until (units in seconds)
      if(i > 1){
        for(k in 1:(i - 1)){
          if(Vdf8[i-k, "Interpolated"] == "Interpolated"){
            Vdf8[i-k, "t_until"] <- abs(as.numeric(difftime(Vdf8[i, "DateTime"], Vdf8[i-k, "DateTime"], units = "secs")))
          } else {
            break
          }
        }
      }
    }
  }
  
  
  
  
  
  
  
  ############# 6. SET SPATIAL EXTENT AND DEFINE SPACE-TIME CUBE #################
  
  # Figure out spatial extent, with whole number of dxy steps
  # Round down/up to nearest multiple of 30000m (regardless of voxel size)
  # Add spatial buffer of 30000m to each extent, to encompass wider disks of uncertainty
  # Lat and Long already transformed into Ycoord and Xcoord (centred flat projection)
  Xmin <- round_any(min(c(Bdf8$Xcoords, Vdf8$Xcoords)), 30000, f = floor) - 30000; Xmin
  Xmax <- round_any(max(c(Bdf8$Xcoords, Vdf8$Xcoords)), 30000, f = ceiling) + 30000; Xmax
  Ymin <- round_any(min(c(Bdf8$Ycoords, Vdf8$Ycoords)), 30000, f = floor) - 30000; Ymin
  Ymax <- round_any(max(c(Bdf8$Ycoords, Vdf8$Ycoords)), 30000, f = ceiling) + 30000; Ymax
  Xrange <- Xmax - Xmin
  Yrange <- Ymax - Ymin
  
  # Determine number of grid cells (pixels) and voxels
  # This assumes Trange is in seconds and XYrange are in m
  n_Xsteps <- Xrange/dxy
  n_Ysteps <- Yrange/dxy
  n_pixels <- n_Xsteps * n_Ysteps
  n_voxels <- n_Tsteps * n_pixels
  
  
  # Create a template raster for 1 time step with all values at 0
  # CRS should match that defined above
  template <- raster(ncol = n_Xsteps, 
                     nrow = n_Ysteps, 
                     xmn = Xmin,
                     ymn = Ymin,
                     xmx = Xmax,
                     ymx = Ymax,
                     crs = Bdf2@proj4string)
  
  
  
  
  
  
  
  ################# 7. CONSTRUCT BIRD PSTP #########################
  
  # Initialise a raster stack (starting with 1 layer)
  BSTP <- template
  
  # Set flag for end of trajectory, just in case last fix is not Real
  stop <- FALSE
  
  # OUTER LOOP: For each time step... (i is index of all steps, j and k are indexes of steps between fixes)
  for(i in 1:n_Tsteps){
    
    # If Real fix
    if(Bdf8[i, "Interpolated"] == "Real"){
      
      # Create single spatial point: Real1 (should overwrite from previous iteration)
      # WGS1984 by default, but this gives option to change
      Real1 <- SpatialPoints(coords = Bdf8[i ,c("Xcoords", "Ycoords")],
                             proj4string = Bdf2@proj4string)
      
      # Convert spatial point into a single raster cell (assuming telemetry error is negligible)
      Real1_raster <- rasterize(Real1, template, field = 1, background = 0)
      
      # Put this point into a layer of the PSTP raster stack
      if(i == 1){
        values(BSTP) <- values(Real1_raster)
      } else {
        BSTP <- addLayer(BSTP, Real1_raster)
      }
      
      # If not the last time step...
      if(i < n_Tsteps){
        # Create raster of distance from Real1 (hopefully in metres)
        Real_dist1 <- distanceFromPoints(template, Real1)
        
        # Find next Real fix
        for(j in 1:(n_Tsteps - i)){
          if(Bdf8[i + j, "Interpolated"] == "Real"){
            
            # Create single spatial point: Real2 (should overwrite from previous iteration)
            Real2 <- SpatialPoints(coords = Bdf8[i + j, c("Xcoords", "Ycoords")],
                                   proj4string = Bdf2@proj4string) # Hopefully this works
            
            # Create raster of distance from Real2 (hopefully in metres)
            Real_dist2 <- distanceFromPoints(template, Real2)
            
            # Break this loop
            break
            
            # If no next Real fix left in trajectory, then set stop = TRUE to break outer loop
            # (In theory, this should never happen because last step has been set to Real)
          } else if (j == n_Tsteps - i){
            stop <- TRUE
            break
            
            # If no Real at this step but still one remaining, look at next time step
          } else {
            next
          }
        }
        
        # Break outer loop if necessary
        if(stop == TRUE){break}
        
        # INNER LOOP: For each time step between Real fixes...
        for(k in 1:(n_Tsteps - i)){
          if(Bdf8[i + k, "Interpolated"] == "Interpolated"){
            # Create single spatial point: Interp_CP
            Interp_CP <- SpatialPoints(coords = Bdf8[i + k, c("Xcoords", "Ycoords")],
                                       proj4string = Bdf2@proj4string) # Hopefully this works
            
            # Define accessible distances from Real1 and Real2 (should be in metres)
            access_dist1 <- Bdf8[i + k, "t_since"] * B_max_velocity
            access_dist2 <- Bdf8[i + k, "t_until"] * B_max_velocity
            
            # Create rasters of accessible areas from Real1 and Real2 (0 or 1)
            Real1_access_area <- Real_dist1 <= access_dist1
            Real2_access_area <- Real_dist2 <= access_dist2
            
            # Create raster of boundary (0 or 1)
            boundary_layer <- Real1_access_area * Real2_access_area
            
            # If boundary layer is entirely 0 (can happen if access_dist < dxy)
            if(sum(values(boundary_layer)) == 0){
              
              # Convert spatial point into a single raster cell of value 1
              prob_STD <- rasterize(Interp_CP, boundary_layer, field = 1, background = 0)
              
              # Else do chosen distance weighting function
            } else {
              
              # Create raster of distance from Interp_CP
              Interp_CP_dist <- distanceFromPoints(template, Interp_CP)
              
              # If centre == 0, fill in centre of Interp_CP_dist to avoid doughnut effect
              Interp_CP_dist[] <- ifelse(Interp_CP_dist[] == 0, dxy, Interp_CP_dist[])
              
              # Apply Gaussian distance weighting function to create probability surface (dBBMMM model)
              # Clipped to boundary layer
              # Normalised so that the sum of all probabilities is STD = 1
              Gaussian_decay <- boundary_layer *
                (1/(2*pi*Bdf8[i + k, "variance"])) *
                exp(-(Interp_CP_dist^2)/(2*Bdf8[i + k, "variance"]))
              prob_STD <- Gaussian_decay / cellStats(Gaussian_decay, "sum")
            }
            
            # Put probability surface values into BSTP raster stack
            BSTP <- addLayer(BSTP, prob_STD)
            
          } else {
            break
          }
        }
      } else {
        break
      }
    } else {
      next
    }
  }
  
  # Change to raster brick
  BSTP <- brick(BSTP)
  
  # Write PSTP to file (change file identifier up top)
  setwd(paste(WD, "PSTP Outputs TEST THIN VOXELS/BV PSTP Bricks/", sep = "")) # CHANGE HERE
  writeRaster(BSTP, filename = paste("BSTP_vox", dxy/1000, "_SubTraj", SubTrajNum_V, ".nc", sep = ""), format = "CDF", overwrite = TRUE)
  setwd(WD)
  
  
  
  
  
  
  
  
  ############## 8. CONSTRUCT VESSEL PSTP ############
  
  # Similar to #7 above, but with Vdf instead of Bdf
  # and using velocity multiplier instead of max velocity, to reflect more consistent movement of ships
  
  # Initialise a raster stack (starting with 1 layer)
  VSTP <- template
  
  # Set flag for end of trajectory, just in case last fix is not Real
  stop <- FALSE
  
  # OUTER LOOP: For each time step... (i is index of all steps, j and k are indexes of steps between fixes)
  for(i in 1:n_Tsteps){
    
    # If Real fix
    if(Vdf8[i, "Interpolated"] == "Real"){
      
      # Create single spatial point: Real1 (should overwrite from previous iteration)
      # WGS1984 by default, but this gives option to change
      Real1 <- SpatialPoints(coords = Vdf8[i ,c("Xcoords", "Ycoords")],
                             proj4string = Vdf2@proj4string)
      
      # Convert spatial point into a single raster cell (assuming telemetry error is negligible)
      Real1_raster <- rasterize(Real1, template, field = 1, background = 0)
      
      # Put this point into a layer of the PSTP raster stack
      if(i == 1){
        values(VSTP) <- values(Real1_raster)
      } else {
        VSTP <- addLayer(VSTP, Real1_raster)
      }
      
      # If not the last time step...
      if(i < n_Tsteps){
        # Create raster of distance from Real1 (in metres)
        Real_dist1 <- distanceFromPoints(template, Real1)
        
        # Find next Real fix
        for(j in 1:(n_Tsteps - i)){
          if(Vdf8[i + j, "Interpolated"] == "Real"){
            
            # Create single spatial point: Real2 (should overwrite from previous iteration)
            Real2 <- SpatialPoints(coords = Vdf8[i + j, c("Xcoords", "Ycoords")],
                                   proj4string = Vdf2@proj4string)
            
            # Create raster of distance from Real2 (hopefully in metres)
            Real_dist2 <- distanceFromPoints(template, Real2)
            
            # Break this loop
            break
            
            # If no next Real fix left in trajectory, then set stop = TRUE to break outer loop
            # (In theory, this should never happen because last step has been set to Real)
          } else if (j == n_Tsteps - i){
            stop <- TRUE
            break
            
            # If no Real at this step but still one remaining, look at next time step
          } else {
            next
          }
        }
        
        # Break outer loop if necessary
        if(stop == TRUE){break}
        
        # INNER LOOP: For each time step between Real fixes...
        for(k in 1:(n_Tsteps - i)){
          if(Vdf8[i + k, "Interpolated"] == "Interpolated"){
            # Create single spatial point: Interp_CP
            Interp_CP <- SpatialPoints(coords = Vdf8[i + k, c("Xcoords", "Ycoords")],
                                       proj4string = Vdf2@proj4string)
            
            # Define accessible distances from Real1 and Real2 (should be in metres)
            # Convert speed to m/s and multiply by velocity multiplier
            access_dist1 <- Vdf8[i + k, "t_since"] * (Vdf8[i+k, "speed"] / 3.6) * Vmultiplier
            access_dist2 <- Vdf8[i + k, "t_until"] * (Vdf8[i+k, "speed"] / 3.6) * Vmultiplier
            #access_dist1 <- Vdf8[i + k, "t_since"] * V_max_velocity
            #access_dist2 <- Vdf8[i + k, "t_until"] * V_max_velocity
            
            # Create rasters of accessible areas from Real1 and Real2 (0 or 1)
            Real1_access_area <- Real_dist1 <= access_dist1
            Real2_access_area <- Real_dist2 <= access_dist2
            
            # Create raster of boundary (0 or 1)
            boundary_layer <- Real1_access_area * Real2_access_area
            
            # If boundary layer is entirely 0 (can happen if access_dist < dxy)
            if(sum(values(boundary_layer)) == 0){
              
              # Convert spatial point into a single raster cell of value 1
              prob_STD <- rasterize(Interp_CP, boundary_layer, field = 1, background = 0)
              
              # Else do chosen distance weighting function
            } else {
              
              # Create raster of distance from Interp_CP
              Interp_CP_dist <- distanceFromPoints(template, Interp_CP)
              
              # If centre == 0, fill in centre of Interp_CP_dist to avoid doughnut effect
              Interp_CP_dist[] <- ifelse(Interp_CP_dist[] == 0, dxy, Interp_CP_dist[])
              
              # Simple option: apply inverse distance weighting function to create probability surface
              # Clipped to boundary layer
              # Normalised so that the sum of all probabilities is STD = 1
              # Downs et al 2014b suggest only slight differences between inverse, linear, and Gaussian functions
              IDW_not_normalised <- boundary_layer * (1/Interp_CP_dist)
              prob_STD <- IDW_not_normalised / cellStats(IDW_not_normalised, "sum")
            }
            
            # Put probability surface values into VSTP raster stack
            VSTP <- addLayer(VSTP, prob_STD)
            
          } else {
            break
          }
        }
      } else {
        break
      }
    } else {
      next
    }
  }
  
  # Change to raster brick
  VSTP <- brick(VSTP)
  
  # Write PSTP to file (change file_identifier up top)
  setwd(paste(WD, "PSTP Outputs TEST THIN VOXELS/BV PSTP Bricks/", sep = "")) # CHANGE HERE
  writeRaster(VSTP, filename = paste("VSTP_vox", dxy/1000, "_SubTraj", SubTrajNum_V, ".nc", sep = ""), format = "CDF", overwrite = TRUE)
  setwd(WD)
  
  
  
  
  
  
  ############### 9. CREATE PSTP SUMMARY TABLE ########################
  
  # Calculate elapsed time
  end_time <- Sys.time()
  elapsed_time <- difftime(end_time, start_time, units = "mins")
  
  # Set WD to Summary Tables folder
  setwd(paste(WD, "/PSTP Outputs TEST THIN VOXELS/Summary Tables/", sep = "")) # CHANGE HERE
  
  # If a summary table already exists...
  if(file.exists("PSTP_summary_table.csv")){
    
    # Read CSV
    PSTP_summary <- read.csv("PSTP_summary_table.csv")
    
    # Add row to table
    PSTP_summary <- PSTP_summary %>% 
      mutate(T_start = ymd_hms(T_start), T_end = ymd_hms(T_end)) %>%
      add_row(SubTraj = Vdf1$SubTraj[1], 
              vesselID = Vdf1$vesselID[1], 
              vox_size_xy = dxy/1000,
              vox_size_t = dt/60,
              CP_Lon = as.numeric(gsub(".*lon_0=|\\s.*", "", projection(BSTP))), 
              CP_Lat = as.numeric(gsub(".*lat_0=|\\s.*", "", projection(BSTP))), 
              T_start = Tmin, 
              T_end = Tmax, 
              n_Tsteps = n_Tsteps, 
              n_pixels = n_pixels, 
              n_voxels = n_voxels, 
              CompTime = abs(as.numeric(elapsed_time)))
    
    # Overwrite CSV
    write.csv(PSTP_summary, "PSTP_summary_table.csv", row.names = FALSE)
    
    # If a summary table does not yet exist...  
  } else {
    
    # Create dataframe
    PSTP_summary <- data.frame(SubTraj = integer(), 
                               vesselID = integer(), 
                               vox_size_xy = integer(),
                               vox_size_t = integer(),
                               CP_Lon = double(), 
                               CP_Lat = double(), 
                               T_start = as.Date(x = integer(0), origin = "1970-01-01"), 
                               T_end = as.Date(x = integer(0), origin = "1970-01-01"), 
                               n_Tsteps = integer(), 
                               n_pixels = integer(), 
                               n_voxels = integer(), 
                               CompTime = double())
    
    # Add row to table
    PSTP_summary <- PSTP_summary %>% 
      add_row(SubTraj = Vdf1$SubTraj[1], 
              vesselID = Vdf1$vesselID[1], 
              vox_size_xy = dxy/1000,
              vox_size_t = dt/60,
              CP_Lon = as.numeric(gsub(".*lon_0=|\\s.*", "", projection(BSTP))), 
              CP_Lat = as.numeric(gsub(".*lat_0=|\\s.*", "", projection(BSTP))), 
              T_start = Tmin, 
              T_end = Tmax, 
              n_Tsteps = n_Tsteps, 
              n_pixels = n_pixels, 
              n_voxels = n_voxels, 
              CompTime = abs(as.numeric(elapsed_time)))
    
    # Write to CSV
    write.csv(PSTP_summary, "PSTP_summary_table.csv", row.names = FALSE)
    
  }
  
  
  
  ##########################################
  
  # Message that iteration is complete
  message(paste("Completed SubTraj ", SubTrajNum_V, " voxel size ", 
                dxy/1000, "x", dxy/1000, "x", dt/60, ".\n", sep = ""))
  
  # Reset WD
  setwd(WD)
  
  # Return elapsed time
  return(abs(as.numeric(elapsed_time)))
  
}
# END OF FUNCTION







############### D. GET B AND V PSTPS USING FUNCTION ########################

# Loop through entire dataset in parallel
pbsapply(1:max(Vdata$SubTraj), get.BV.PSTP, dt = 1800, dxy = 30000, cl = cl) # 30x30x30
pbsapply(1:max(Vdata$SubTraj), get.BV.PSTP, dt = 300, dxy = 5000, cl = cl) # 5x5x5
stopCluster(cl)

# Loop through entire dataset (ORIGINAL)
for(subtraj in 1:max(Vdata$SubTraj)){
  
  # 30x30x30 voxel
  get.BV.PSTP(subtraj, 1800, 30000)
  cat("Completed SubTraj ", subtraj, " voxel size 30.", "\n", sep = "")
  
  # 5x5x5 voxel
  get.BV.PSTP(subtraj, 300, 5000)
  cat("Completed SubTraj ", subtraj, " voxel size 5.", "\n", sep = "")
  
}
