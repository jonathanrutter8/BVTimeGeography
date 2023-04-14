########################################################################################
################  PAIRWISE INTERACTIONS OF BIRDS AND VESSELS  ##########################
########################################################################################

# Takes in CSVs of bird and vessel data
# You must have already produced BV PSTPs within your "PSTP Outputs FINAL" folder
# Adjust section B and function arguments as required
# Double check input column names, as these could impact creating regular trajectories
#   See Bdf4 and Vdf4 objects within section D
# This script is set up to run in parallel

####################### CONTENTS ###########################

# A. LOAD PACKAGES
# B. DEFINE GLOBAL PARAMETERS
# C. DEFINE FUNCTION FOR VOXEL-BASED OVERLAP INTERACTIONS
  # 1. Attendance
  # 2. Encounter
  # 3. Summary table
  # D. DEFINE FUNCTION FOR DISTANCE THRESHOLD-BASED OVERLAP INTERACTIONS
# PART I: GET PSTP WITH 1x1x1 VOXELS
  # 1. Subset B and V trajectories
  # 2. Set coordinate reference system
  # 3. Set temporal extent and time steps
  # 4. Create regular trajectory for bird (B)
  # 5. Create regular trajectory for vessel (V)
  # 6. Set spatial extent and define space-time cube
  # 7. Construct bird PSTP
  # 8. Construct vessel PSTP
  # 9. Create PSTP Summary Table
# PART II: GET INTERACTION BRICKS
  # 1. Attendance
  # 2. Encounter
  # 3. Summary table
# E. GET INTERACTION BRICKS USING FUNCTIONS


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


######### B. DEFINE GLOBAL PARAMETERS ###############

# Set WD
WD <- "____________________"
setwd(WD)

# Set distance thresholds (in metres)
# 30km for encounters, 5km for attendance
d0_attendance <- 5000
d0_encounter <- 30000

# Set probability threshold (above which an interaction event is considered to occur)
p0_interaction <- 0.5
p0_lower <- 0.025 # Lower and upper bounds to help quantify uncertainty
p0_upper <- 0.975

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
  
  # Set WD
  WD <- "C:/Users/jonat/Desktop/University of Oxford BCM/Oxford BCM 2021-2022/BCM Dissertation/R Scripts and Data/"
  setwd(WD)
  
  # Set distance thresholds (in metres)
  # 30km for encounters, 5km for attendance
  d0_attendance <- 5000
  d0_encounter <- 30000
  
  # Set probability threshold (above which an interaction event is considered to occur)
  p0_interaction <- 0.5
  p0_lower <- 0.025
  p0_upper <- 0.975
  
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





########### C. DEFINE FUNCTION FOR VOXEL-BASED OVERLAP INTERACTIONS #########

# This function creates 1 attendance and 1 encounter raster brick per SubTraj.
# Each raster brick shows probability of interaction calculated from voxel-based overlap.
# Both netCDF files are saved in the "Interaction Bricks" folder.
# Probability of interaction saved as new column to B Subset CSV.
# Expected duration of interaction saved in summary CSV.
# Function returns the elapsed time.

# Function requires "BV PSTP Bricks" folder containing correctly labelled netCDF files.
# Single parameter: SubTrajNum_V, the vessel trajectory subset (can loop over this if needed)

# START OF FUNCTION
get.vox.interactions <- function(SubTrajNum_V){
  
  # Set WD
  setwd(paste(WD, "PSTP Outputs TEST THIN VOXELS/", sep = "")) # CHANGE HERE
  
  # Start the clock
  start_time <- Sys.time()
  
  ######### 1. ATTENDANCE: OVERLAP OF 5x5x5 VOXELS #########
  
  # Read B and V netCDF files for chosen SubTraj
  BV_filepair <- list.files("BV PSTP Bricks", 
                            pattern = paste("vox5_SubTraj", SubTrajNum_V, ".nc", sep = ""))
  if(length(BV_filepair) < 2){return(0)}
  B_index <- which(grepl("BSTP", BV_filepair))
  V_index <- which(grepl("VSTP", BV_filepair))
  BSTP1 <- brick(paste("BV PSTP Bricks", BV_filepair[B_index], sep = "/"))
  VSTP1 <- brick(paste("BV PSTP Bricks", BV_filepair[V_index], sep = "/"))
  
  # Remake template layer
  template1 <- raster(BSTP1, layer = 1)
  values(template1) <- 0
  
  # Initialise BV overlap raster brick
  prob_BVOL <- template1
  
  # CALCULATE OVERLAP FOR CHOSEN VOXEL SIZE: For each time step...
  for(i in 1:nlayers(BSTP1)){
    
    # Extract raster layers for B and V
    Blayer <- raster(BSTP1, layer = i)
    Vlayer <- raster(VSTP1, layer = i)
    
    # Create raster layer with probability of overlap
    BVOL_layer <- Blayer * Vlayer
    
    # Add to BV overlap raster brick
    if(i == 1){
      values(prob_BVOL) <- values(BVOL_layer)
    } else {
      prob_BVOL <- addLayer(prob_BVOL, BVOL_layer)
    }
  }
  
  # Change to raster brick
  prob_BVOL <- brick(prob_BVOL)
  
  # Save results as netCDF
  setwd(paste(WD, "PSTP Outputs TEST THIN VOXELS/Interaction Bricks", sep = "")) # CHANGE HERE
  writeRaster(prob_BVOL, filename = paste("BVinteraction_vox5_SubTraj", SubTrajNum_V, ".nc", sep = ""), format = "CDF", overwrite = TRUE)
  setwd(paste(WD, "PSTP Outputs TEST THIN VOXELS/", sep = "")) # CHANGE HERE
  
  # Make vector of probabilities
  vector_OL <- c()
  vector_OL[1] <- cellStats(prob_BVOL[[1]], "sum")
  for(step in 2:nlayers(prob_BVOL)){
    vector_OL[step] <- cellStats(prob_BVOL[[step]], "sum")
  }
  
  # Add probability of attendance to B Subset Dataframe
  setwd(paste(WD, "PSTP Outputs TEST THIN VOXELS/B Subset Dataframes", sep = "")) # CHANGE HERE
  Bdf_att <- read.csv(paste("B_RegTraj_vox5_SubTraj", SubTrajNum_V, ".csv", sep = ""))
  Bdf_att <- Bdf_att %>% mutate(prob_att = vector_OL)
  write.csv(Bdf_att, paste("B_RegTraj_vox5_SubTraj", SubTrajNum_V, ".csv", sep = ""))
  setwd(paste(WD, "PSTP Outputs TEST THIN VOXELS/", sep = "")) # CHANGE HERE
  
  # Get expected duration of attendance (in minutes)
  exp_dur_attendance <- 5 * sum(vector_OL >= p0_interaction)
  exp_dur_attendance_lower <- 5 * sum(vector_OL >= p0_lower)
  exp_dur_attendance_upper <- 5 * sum(vector_OL >= p0_upper)
  
  
  
  
  
  
  ######### 2. ENCOUNTER: OVERLAP OF 30x30x30 VOXELS #########
  
  # Read B and V netCDF files for chosen SubTraj
  BV_filepair <- list.files("BV PSTP Bricks", 
                            pattern = paste("vox30_SubTraj", SubTrajNum_V, sep = ""))
  if(length(BV_filepair) < 2){return(0)}
  B_index <- which(grepl("BSTP", BV_filepair))
  V_index <- which(grepl("VSTP", BV_filepair))
  BSTP1 <- brick(paste("BV PSTP Bricks", BV_filepair[B_index], sep = "/"))
  VSTP1 <- brick(paste("BV PSTP Bricks", BV_filepair[V_index], sep = "/"))
  
  # Remake template layer
  template1 <- raster(BSTP1, layer = 1)
  values(template1) <- 0
  
  # Initialise BV overlap raster brick
  prob_BVOL <- template1
  
  # CALCULATE OVERLAP FOR CHOSEN VOXEL SIZE: For each time step...
  for(i in 1:nlayers(BSTP1)){
    
    # Extract raster layers for B and V
    Blayer <- raster(BSTP1, layer = i)
    Vlayer <- raster(VSTP1, layer = i)
    
    # Create raster layer with probability of overlap
    BVOL_layer <- Blayer * Vlayer
    
    # Add to BV overlap raster brick
    if(i == 1){
      values(prob_BVOL) <- values(BVOL_layer)
    } else {
      prob_BVOL <- addLayer(prob_BVOL, BVOL_layer)
    }
  }
  
  # Change to raster brick
  prob_BVOL <- brick(prob_BVOL)
  
  # Save results as netCDF
  setwd(paste(WD, "PSTP Outputs TEST THIN VOXELS/Interaction Bricks", sep = "")) # CHANGE HERE
  writeRaster(prob_BVOL, filename = paste("BVinteraction_vox30_SubTraj", SubTrajNum_V, ".nc", sep = ""), format = "CDF", overwrite = TRUE)
  setwd(paste(WD, "PSTP Outputs TEST THIN VOXELS/", sep = "")) # CHANGE HERE
  
  # Make vector of probabilities
  vector_OL <- c()
  vector_OL[1] <- cellStats(prob_BVOL[[1]], "sum")
  for(step in 2:nlayers(prob_BVOL)){
    vector_OL[step] <- cellStats(prob_BVOL[[step]], "sum")
  }
  
  # Add probability of encounter to B Subset Dataframe
  setwd(paste(WD, "PSTP Outputs TEST THIN VOXELS/B Subset Dataframes", sep = "")) # CHANGE HERE
  Bdf_enc <- read.csv(paste("B_RegTraj_vox30_SubTraj", SubTrajNum_V, ".csv", sep = ""))
  Bdf_enc <- Bdf_enc %>% mutate(prob_enc = vector_OL)
  write.csv(Bdf_enc, paste("B_RegTraj_vox30_SubTraj", SubTrajNum_V, ".csv", sep = ""))
  setwd(paste(WD, "PSTP Outputs TEST THIN VOXELS/", sep = "")) # CHANGE HERE
  
  # Get expected duration of encounter (in minutes)
  exp_dur_encounter <- 30 * sum(vector_OL >= p0_interaction)
  exp_dur_encounter_lower <- 30 * sum(vector_OL >= p0_lower)
  exp_dur_encounter_upper <- 30 * sum(vector_OL >= p0_upper)
  
  
  
  
  
  
  ############### 3. CREATE VOX INTERACTIONS SUMMARY TABLE ########################
  
  # Get elapsed time
  end_time <- Sys.time()
  elapsed_time <- difftime(end_time, start_time, units = "secs")
  
  # Set WD to Summary Tables folder
  setwd(paste(WD, "/PSTP Outputs TEST THIN VOXELS/Summary Tables/", sep = "")) # CHANGE HERE
  
  # If a summary table already exists...
  if(file.exists("Interactions_vox_summary_table.csv")){
    
    # Read CSV
    Interactions_vox_summary <- read.csv("Interactions_vox_summary_table.csv")
    
    # Add row to table
    Interactions_vox_summary <- Interactions_vox_summary %>% 
      add_row(SubTraj = SubTrajNum_V, 
              vesselID = Vdata$vesselID[Vdata$SubTraj == SubTrajNum_V][1], 
              birdID = Vdata$birdID[Vdata$SubTraj == SubTrajNum_V][1],
              Duration_Encounter = exp_dur_encounter,
              Duration_Enc_lower = exp_dur_encounter_lower,
              Duration_Enc_upper = exp_dur_encounter_upper,
              Duration_Attendance = exp_dur_attendance,
              Duration_Att_lower = exp_dur_attendance_lower,
              Duration_Att_upper = exp_dur_attendance_upper,
              Proportion_Attendance = exp_dur_attendance / exp_dur_encounter,
              CompTime = as.numeric(elapsed_time))
    
    # Overwrite CSV
    write.csv(Interactions_vox_summary, "Interactions_vox_summary_table.csv", row.names = FALSE)
    
  } else {
    
    # Create dataframe
    Interactions_vox_summary <- data.frame(SubTraj = integer(), 
                                           vesselID = integer(),
                                           birdID = integer(),
                                           Duration_Encounter = double(),
                                           Duration_Enc_lower = double(),
                                           Duration_Enc_upper = double(),
                                           Duration_Attendance = double(),
                                           Duration_Att_lower = double(),
                                           Duration_Att_upper = double(),
                                           Proportion_Attendance = double(),
                                           CompTime = double())
    
    # Add row to table
    Interactions_vox_summary <- Interactions_vox_summary %>% 
      add_row(SubTraj = SubTrajNum_V, 
              vesselID = Vdata$vesselID[Vdata$SubTraj == SubTrajNum_V][1], 
              birdID = Vdata$birdID[Vdata$SubTraj == SubTrajNum_V][1],
              Duration_Encounter = exp_dur_encounter,
              Duration_Enc_lower = exp_dur_encounter_lower,
              Duration_Enc_upper = exp_dur_encounter_upper,
              Duration_Attendance = exp_dur_attendance,
              Duration_Att_lower = exp_dur_attendance_lower,
              Duration_Att_upper = exp_dur_attendance_upper,
              Proportion_Attendance = exp_dur_attendance / exp_dur_encounter,
              CompTime = as.numeric(elapsed_time))
    
    # Write to CSV
    write.csv(Interactions_vox_summary, "Interactions_vox_summary_table.csv", row.names = FALSE)
    
  }
  
  
  ###################
  
  # Message that iteration is complete
  message(paste("Completed SubTraj ", SubTrajNum_V, " VOB interactions.\n", sep = ""))
  
  # Reset WD
  setwd(WD)  
  
  # Return elapsed time
  return(abs(as.numeric(elapsed_time)))
  
}
# END OF FUNCTION










########## D. DEFINE FUNCTION FOR DISTANCE THRESHOLD-BASED INTERACTIONS #############

# Should loop over SubTrajs with known overlap, based on voxel-based overlap analysis above

# START OF FUNCTION
get.dist.interactions <- function(SubTrajNum_V){
  
  ################ PART I: GET PSTP WITH 1x1x1 VOXELS ######################################
  
  ################### 1. SUBSET B AND V TRAJECTORIES ##########################
  
  # Start the clock
  start_time <- Sys.time()
  
  # Set WD
  setwd(paste(WD, "PSTP Outputs FINAL/", sep = ""))
  
  # Set 1x1x1 voxel size
  dt <- 60
  dxy <- 1000
  
  # Find first/last indices with potential interaction based on vox-based overlap (30km voxel)
  Bdf_OL <- read.csv(paste("B Subset Dataframes/B_RegTraj_vox30_SubTraj", SubTrajNum_V, ".csv", sep = ""))
  first_index_OL <- min(which(Bdf_OL$prob_enc > 0))
  last_index_OL <- max(which(Bdf_OL$prob_enc > 0))
  
  # Buffer the min and max times by 60 mins
  #first_DT_OL <- ymd_hms(Bdf_OL$DateTime[first_index_OL]) - minutes(60)
  #last_DT_OL <- ymd_hms(Bdf_OL$DateTime[last_index_OL]) + minutes(60)
  
  # Narrow down to a portion where SubTraj 9 and 14 are near each other
  first_DT_OL <- ymd_hm("2019/07/10 02:30")
  last_DT_OL <- ymd_hm("2019/07/10 09:30")
  
  # Get portions V trajectory between these times
  Vdf1 <- Vdata %>%
    mutate(birdID = as.factor(birdID), DateTime = ymd_hms(DateTime)) %>%
    filter(SubTraj == SubTrajNum_V, DateTime >= first_DT_OL, DateTime <= last_DT_OL)
  birdID_Vdf1 <- Vdf1$birdID[1]
  vesselID_Vdf1 <- Vdf1$vesselID[1]
  
  # Get portions of B trajectory between these times
  Bdf1 <- Bdata %>%
    filter(ID == birdID_Vdf1) %>%
    mutate(ID = as.factor(ID), DateTime = ymd_hms(DateTime)) %>%
    filter(DateTime >= first_DT_OL, DateTime <= last_DT_OL)
  
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
  setwd(paste(WD, "/PSTP Outputs FINAL/B Subset Dataframes/", sep = ""))
  write.csv(Bdf8, paste("B_RegTraj_vox", dxy/1000, "_SubTraj", SubTrajNum_V, "a.csv", sep = ""))
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
  setwd(paste(WD, "/PSTP Outputs FINAL/BV PSTP Bricks/", sep = ""))
  writeRaster(BSTP, filename = paste("BSTP_vox", dxy/1000, "_SubTraj", SubTrajNum_V, "a.nc", sep = ""), format = "CDF", overwrite = TRUE)
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
  setwd(paste(WD, "/PSTP Outputs FINAL/BV PSTP Bricks/", sep = ""))
  writeRaster(VSTP, filename = paste("VSTP_vox", dxy/1000, "_SubTraj", SubTrajNum_V, "a.nc", sep = ""), format = "CDF", overwrite = TRUE)
  setwd(WD)
  
  # Calculate elapsed time
  end_time <- Sys.time()
  elapsed_time <- difftime(end_time, start_time, units = "secs")
  
  
  
  
  
  ############### 9. CREATE PSTP SUMMARY TABLE ########################
  
  # Set WD to Summary Tables folder
  setwd(paste(WD, "/PSTP Outputs FINAL/Summary Tables/", sep = ""))
  
  # If a summary table already exists...
  if(file.exists("PSTP_summary_table.csv")){
    
    # Read CSV
    PSTP_summary <- read.csv("PSTP_summary_table.csv")
    
    # Add row to table
    PSTP_summary <- PSTP_summary %>% 
      mutate(T_start = ymd_hms(T_start), T_end = ymd_hms(T_end)) %>%
      add_row(SubTraj = Vdf1$SubTraj[1], 
              vesselID = Vdf1$vesselID[1], 
              voxel_size = dxy/1000,
              CP_Lon = as.numeric(gsub(".*lon_0=|\\s.*", "", projection(BSTP))), 
              CP_Lat = as.numeric(gsub(".*lat_0=|\\s.*", "", projection(BSTP))), 
              T_start = Tmin, 
              T_end = Tmax, 
              n_Tsteps = n_Tsteps, 
              n_pixels = n_pixels, 
              n_voxels = n_voxels, 
              CompTime = as.numeric(elapsed_time))
    
    # Overwrite CSV
    write.csv(PSTP_summary, "PSTP_summary_table.csv", row.names = FALSE)
    
  } else {
    
    # Create dataframe
    PSTP_summary <- data.frame(SubTraj = integer(), 
                               vesselID = integer(), 
                               voxel_size = integer(),
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
              voxel_size = dxy/1000,
              CP_Lon = as.numeric(gsub(".*lon_0=|\\s.*", "", projection(BSTP))), 
              CP_Lat = as.numeric(gsub(".*lat_0=|\\s.*", "", projection(BSTP))), 
              T_start = Tmin, 
              T_end = Tmax, 
              n_Tsteps = n_Tsteps, 
              n_pixels = n_pixels, 
              n_voxels = n_voxels, 
              CompTime = as.numeric(elapsed_time))
    
    # Write to CSV
    write.csv(PSTP_summary, "PSTP_summary_table.csv", row.names = FALSE)
    
  }
  
  
  
  
  
  ######### PART II: GET INTERACTION BRICKS ##################
  
  
  
  # Start the clock
  start_time <- Sys.time()
  
  # Initialise BV attendance/encounter raster stacks
  prob_attendance <- template
  prob_encounter <- template
  
  ########### 1. ATTENDANCE: WITHIN d0 = 5km ###############
  
  # For each time step...
  for(i in 1:nlayers(BSTP)){
    
    # Progress bar
    progress(i, progress.bar = TRUE)
    
    # Initialise BV attendance raster disk
    attendance_disk <- template
    
    # Extract raster layers for B and V
    Blayer <- raster(BSTP, layer = i)
    Vlayer <- raster(VSTP, layer = i)
    
    # Get vector of indexes for where P > 0 (within prism boundary)
    boundary <- which(values(Blayer) > 0)
    
    # For each B cell where P > 0...
    for(j in boundary){
      
      # Create raster of V probabilities, clipped to within d0 of B cell
      xy <- xyFromCell(Blayer, j)
      V_within_d0 <- distanceFromPoints(Vlayer, xy) <= d0_attendance # 0 or 1
      Vprob_within_d0 <- Vlayer * V_within_d0
      
      # Multiply P(B) * sum(P(V) within d0, get single P value
      prob_attendance_1cell <- extract(Blayer, j) * cellStats(Vprob_within_d0, "sum")
      
      # Put value into attendance_disk at that cell
      attendance_disk[j] <- prob_attendance_1cell
    }
    
    # Add new layer to BV attendance raster brick
    if(i == 1){
      values(prob_attendance) <- values(attendance_disk)
    } else {
      prob_attendance <- addLayer(prob_attendance, attendance_disk)
    }
  }
  
  # Change to raster brick
  prob_attendance <- brick(prob_attendance)
  
  # Save results as netCDF
  setwd(paste(WD, "PSTP Outputs FINAL/Interaction Bricks", sep = ""))
  writeRaster(prob_attendance, filename = paste("BVinteraction_dist5_SubTraj", SubTrajNum_V, "a.nc", sep = ""), format = "CDF", overwrite = TRUE)
  setwd(paste(WD, "PSTP Outputs FINAL/", sep = ""))
  
  # Make vector of probabilities
  vector_interaction <- c()
  vector_interaction[1] <- cellStats(prob_attendance[[1]], "sum")
  for(step in 2:nlayers(prob_attendance)){
    vector_interaction[step] <- cellStats(prob_attendance[[step]], "sum")
  }
  
  # Add probability of attendance to B Subset Dataframe
  setwd(paste(WD, "PSTP Outputs FINAL/B Subset Dataframes/", sep = ""))
  Bdf_att <- read.csv(paste("B_RegTraj_vox1_SubTraj", SubTrajNum_V, "a.csv", sep = ""))
  Bdf_att <- Bdf_att %>% mutate(prob_att = vector_interaction)
  write.csv(Bdf_att, paste("B_RegTraj_vox1_SubTraj", SubTrajNum_V, "a.csv", sep = ""))
  setwd(paste(WD, "PSTP Outputs FINAL/", sep = ""))
  
  # Get expected duration of encounter
  exp_dur_attendance <- sum(vector_interaction >= p0_interaction)
  exp_dur_attendance_lower <- sum(vector_interaction >= p0_lower)
  exp_dur_attendance_upper <- sum(vector_interaction >= p0_upper)
  
  
  
  
  ############## 2. ENCOUNTER: WITHIN d0 = 30km ##################
  
  # For each time step...
  for(i in 1:nlayers(BSTP)){
    
    # Progress bar
    progress(i, progress.bar = TRUE)
    
    # Initialise BV encounter raster disk
    encounter_disk <- template
    
    # Extract raster layers for B and V
    Blayer <- raster(BSTP, layer = i)
    Vlayer <- raster(VSTP, layer = i)
    
    # Get vector of indexes for where P > 0 (within prism boundary)
    boundary <- which(values(Blayer) > 0)
    
    # For each B cell where P > 0...
    for(j in boundary){
      
      # Create raster of V probabilities, clipped to within d0 of B cell
      xy <- xyFromCell(Blayer, j)
      V_within_d0 <- distanceFromPoints(Vlayer, xy) <= d0_encounter # 0 or 1
      Vprob_within_d0 <- Vlayer * V_within_d0
      
      # Multiply P(B) * sum(P(V) within d0, get single P value
      prob_encounter_1cell <- extract(Blayer, j) * cellStats(Vprob_within_d0, "sum")
      
      # Put value into encounter_disk at that cell
      encounter_disk[j] <- prob_encounter_1cell
    }
    
    # Add new layer to BV encounter raster brick
    if(i == 1){
      values(prob_encounter) <- values(encounter_disk)
    } else {
      prob_encounter <- addLayer(prob_encounter, encounter_disk)
    }
  }
  
  # Change to raster brick
  prob_encounter <- brick(prob_encounter)
  
  # Save results as netCDF
  setwd(paste(WD, "PSTP Outputs FINAL/Interaction Bricks", sep = ""))
  writeRaster(prob_encounter, filename = paste("BVinteraction_dist30_SubTraj", SubTrajNum_V, "a.nc", sep = ""), format = "CDF", overwrite = TRUE)
  setwd(paste(WD, "PSTP Outputs FINAL/", sep = ""))
  
  # Make vector of probabilities
  vector_interaction <- c()
  vector_interaction[1] <- cellStats(prob_encounter[[1]], "sum")
  for(step in 2:nlayers(prob_encounter)){
    vector_interaction[step] <- cellStats(prob_encounter[[step]], "sum")
  }
  
  # Add probability of encounter to B Subset Dataframe
  setwd(paste(WD, "PSTP Outputs FINAL/B Subset Dataframes", sep = ""))
  Bdf_enc <- read.csv(paste("B_RegTraj_vox1_SubTraj", SubTrajNum_V, "a.csv", sep = ""))
  Bdf_enc <- Bdf_enc %>% mutate(prob_enc = vector_interaction)
  write.csv(Bdf_enc, paste("B_RegTraj_vox1_SubTraj", SubTrajNum_V, "a.csv", sep = ""))
  setwd(paste(WD, "PSTP Outputs FINAL/", sep = ""))
  
  # Get expected duration of encounter
  exp_dur_encounter <- sum(vector_interaction >= p0_interaction)
  exp_dur_encounter_lower <- sum(vector_interaction >= p0_lower)
  exp_dur_encounter_upper <- sum(vector_interaction >= p0_upper)
  
  
  
  
  
  ############### 3. CREATE DIST INTERACTIONS SUMMARY TABLE ########################
  
  # Get elapsed time
  end_time <- Sys.time()
  elapsed_time <- difftime(end_time, start_time, units = "mins")
  
  # Set WD to Summary Tables folder
  setwd(paste(WD, "/PSTP Outputs FINAL/Summary Tables/", sep = ""))
  
  # If a summary table already exists...
  if(file.exists("Interactions_dist_summary_table.csv")){
    
    # Read CSV
    Interactions_dist_summary <- read.csv("Interactions_dist_summary_table.csv")
    
    # Add row to table
    Interactions_dist_summary <- Interactions_dist_summary %>% 
      add_row(SubTraj = SubTrajNum_V, 
              vesselID = Vdf1$vesselID[1], 
              birdID = Vdata$birdID[Vdata$SubTraj == SubTrajNum_V][1],
              Duration_Encounter = exp_dur_encounter,
              Duration_Enc_lower = exp_dur_encounter_lower,
              Duration_Enc_upper = exp_dur_encounter_upper,
              Duration_Attendance = exp_dur_attendance,
              Duration_Att_lower = exp_dur_attendance_lower,
              Duration_Att_upper = exp_dur_attendance_upper,
              Proportion_Attendance = exp_dur_attendance / exp_dur_encounter,
              CompTime = abs(as.numeric(elapsed_time)))
    
    # Overwrite CSV
    write.csv(Interactions_dist_summary, "Interactions_dist_summary_table.csv", row.names = FALSE)
    
  } else {
    
    # Create dataframe
    Interactions_dist_summary <- data.frame(SubTraj = integer(), 
                                            vesselID = integer(),
                                            birdID = integer(),
                                            Duration_Encounter = double(),
                                            Duration_Enc_lower = double(),
                                            Duration_Enc_upper = double(),
                                            Duration_Attendance = double(),
                                            Duration_Att_lower = double(),
                                            Duration_Att_upper = double(),
                                            Proportion_Attendance = double(),
                                            CompTime = double())
    
    # Add row to table
    Interactions_dist_summary <- Interactions_dist_summary %>% 
      add_row(SubTraj = SubTrajNum_V, 
              vesselID = Vdf1$vesselID[1], 
              birdID = Vdata$birdID[Vdata$SubTraj == SubTrajNum_V][1],
              Duration_Encounter = exp_dur_encounter,
              Duration_Enc_lower = exp_dur_encounter_lower,
              Duration_Enc_upper = exp_dur_encounter_upper,
              Duration_Attendance = exp_dur_attendance,
              Duration_Att_lower = exp_dur_attendance_lower,
              Duration_Att_upper = exp_dur_attendance_upper,
              Proportion_Attendance = exp_dur_attendance / exp_dur_encounter,
              CompTime = abs(as.numeric(elapsed_time)))
    
    # Write to CSV
    write.csv(Interactions_dist_summary, "Interactions_dist_summary_table.csv", row.names = FALSE)
    
  }
  
  ###################
  
  # Message that iteration is complete
  message(paste("Completed SubTraj ", SubTrajNum_V, " DTB interactions.\n", sep = ""))
  
  # Reset WD
  setwd(WD)  
  
  # Return elapsed time
  return(abs(as.numeric(elapsed_time)))
  
}
# END OF FUNCTION





################### E. GET INTERACTION BRICKS USING FUNCTIONS #####################


########## Loop through dataset (same as below, but with sapply) #############

# Voxel overlap-based interaction (all SubTrajs)
pbsapply(1:max(Vdata$SubTraj), get.vox.interactions, cl = cl)

# Distance threshold-based interaction (only feed in known SubTrajs with overlap)
pbsapply(c(9, 14), get.dist.interactions, cl = cl)

# Close the cluster
stopCluster(cl)



