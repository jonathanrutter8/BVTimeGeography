########################################################################################
#####################  MSC DISSERTATION RESULTS OUTPUTS  ###############################
########################################################################################

# Outputs are in order of manuscript results
# Sections E and F create time series and maps

####################### CONTENTS ###########################

# A. Load packages
# B. Set WD and global parameters
# C. Bird and vessel trajectories
# D. Creation of voxel-based PSTPs
# E. Voxel overlap-based interaction
# F. Distance threshold-based interaction
# G. Duration of interaction


############# A. LOAD PACKAGES #############################

library(tidyverse)
library(lubridate)
library(ggmap)
library(sp)
library(move)
library(raster)
library(svMisc)
library(ncdf4)
library(viridis)
#library(cubeview)



############## B. SET WD AND GLOBAL PARAMETERS #################

# Set WD
WD <- "_____________"
setwd(WD)

# Set CRS of input data: WGS 1984
WGS1984 <- CRS("+init=epsg:4326")




################# C. BIRD AND VESSEL TRAJECTORIES ###############

# View study area with basemap in ArcGIS Pro

# Load in bird and vessel trajectories
Bdata <- read.csv("Rainier Clean Clipped/Rainier_Clip_Dissertation.csv")
Vdata <- read.csv("GFW Queried Data/GFW_5daytest_clean.csv")

# How many SubTraj?
TotalSubTrajNum <- max(Vdata$SubTraj)

# Subset bird dataframe
Bdf34 <- Bdata %>%
  mutate(ID = as.factor(ID), DateTime = ymd_hms(DateTime)) %>%
  filter(ID == 64334, month(DateTime) == 7, day(DateTime) >= 9, day(DateTime) <= 11)
Bdf36 <- Bdata %>%
  mutate(ID = as.factor(ID), DateTime = ymd_hms(DateTime)) %>%
  filter(ID == 64336, month(DateTime) == 7, day(DateTime) >= 8, day(DateTime) <= 12)

# Subset vessel dataframe
Vdf34 <- Vdata %>%
  mutate(birdID = as.factor(birdID), vesselID = as.factor(vesselID), DateTime = ymd_hms(DateTime)) %>%
  filter(birdID == 64334, month(DateTime) == 7, day(DateTime) >= 9, day(DateTime) <= 11)
Vdf36 <- Vdata %>%
  mutate(birdID = as.factor(birdID), vesselID = as.factor(vesselID), DateTime = ymd_hms(DateTime)) %>%
  filter(birdID == 64336, month(DateTime) == 7, day(DateTime) >= 8, day(DateTime) <= 12)

# Turn trajectories into move objects
Bmv34 <- move(x = Bdf34$Long,
              y = Bdf34$Lat,
              time = Bdf34$DateTime,
              proj = WGS1984,
              animal = Bdf34$ID,
              data = Bdf34,
              removeDuplicatedTimestamps = TRUE)
Bmv36 <- move(x = Bdf36$Long,
              y = Bdf36$Lat,
              time = Bdf36$DateTime,
              proj = WGS1984,
              animal = Bdf36$ID,
              data = Bdf36,
              removeDuplicatedTimestamps = TRUE)
Vmv34 <- move(x = Vdf34$Long,
              y = Vdf34$Lat,
              time = Vdf34$DateTime,
              proj = WGS1984,
              animal = Vdf34$vesselID,
              data = Vdf34,
              removeDuplicatedTimestamps = TRUE)
Vmv36 <- move(x = Vdf36$Long,
              y = Vdf36$Lat,
              time = Vdf36$DateTime,
              proj = WGS1984,
              animal = Vdf36$vesselID,
              data = Vdf36,
              removeDuplicatedTimestamps = TRUE)

# Plot 64334 bird and vessel trajectories
plot(Vmv34, xlab = "Longitude", ylab = "Latitude", type = "l", 
     col = palette(rainbow(11)), 
     pch = 16, lwd = 1.8)
lines(Bmv34, pch = 16, lwd = 3.5, col = "blue")
points(Bmv34, pch = 20, cex = 1.8)
legend("topright", 
       c("Toroa Blue-61b", sub(".", "Vessel ", levels(Vmv34@trackId))),
       pch = rep(NA, 12),
       col = c("blue", palette(rainbow(11))),
       lty = rep(1, 12), lwd = c(3.5, rep(1.8, 11)),
       cex = 0.8, text.col = "white", bg = "white")
legend("topright", 
       c("Toroa Blue-61b", sub(".", "Vessel ", levels(Vmv34@trackId))),
       pch = c(20, rep(NA, 11)),
       pt.cex = c(1.8, rep(NA, 11)),
       col = c("black", rep(NA, 11)),
       lty = c(NA, rep(1, 11)), lwd = c(NA, rep(1.8, 11)),
       cex = 0.8, text.col = "black", bty = "n")

# Plot 64336 bird and vessel trajectories
plot(Vmv36, xlab = "Longitude", ylab = "Latitude", type = "l", pch = 16, lwd = 1.8)
lines(Bmv36, pch = 16, lwd = 3.5, col = "blue")
points(Bmv36, pch = 20, cex = 1.8)
legend("topright", 
       c("Toroa Blue-07b", sub(".", "Vessel ", levels(Vmv36@trackId))),
       pch = rep(NA, 5),
       col = c("blue", palette(rainbow(4))),
       lty = rep(1, 5), lwd = c(3.5, rep(1.8, 4)),
       cex = 0.8, text.col = "white", bg = "white")
legend("topright", 
       c("Toroa Blue-07b", sub(".", "Vessel ", levels(Vmv36@trackId))),
       pch = c(20, rep(NA, 4)),
       pt.cex = c(1.8, rep(NA, 4)),
       col = c("black", rep(NA, 4)),
       lty = c(NA, rep(1, 4)), lwd = c(NA, rep(1.8, 4)),
       cex = 0.8, text.col = "black", bty = "n")

# Plot 64334 variance
ggplot(Bdf34, aes(x = DateTime, y = sigma1_2)) +
  geom_area(fill = "darkgrey") +
  theme_bw() +
  labs(x = "DateTime (UTC)", y = "dBBMM Motion Variance (m^2/s)", title = "Blue-61b") +
  theme(plot.title = element_text(vjust = - 10, hjust = 0.95, size = 18))

# Plot 64336 variance
ggplot(Bdf36, aes(x = DateTime, y = sigma1_2)) +
  geom_area(fill = "darkgrey") +
  theme_bw() +
  labs(x = "DateTime (UTC)", y = "dBBMM Motion Variance (m^2/s)", title = "Blue-07b") +
  theme(plot.title = element_text(vjust = - 10, hjust = 0.95, size = 18))




############### D. CREATION OF VOXEL-BASED PSTPS ##################

# Set WD
setwd(paste(WD, "PSTP Outputs FINAL2/", sep = ""))

# Read in raster bricks for 1 bird-vessel pairing (R64334, SubTraj 9a)
BSTP9a <- brick("BV PSTP Bricks/BSTP_vox1_SubTraj9a.nc")
VSTP9a <- brick("BV PSTP Bricks/VSTP_vox1_SubTraj9a.nc")


# Read brick and CSV for SubTraj 9a, vox5
BSTP9a <- brick("BV PSTP Bricks/BSTP_vox1_SubTraj9a.nc")
BSTP9a[BSTP9a == 0] <- NA
Bdf9a <- read.csv("B Subset Dataframes/B_RegTraj_vox1_SubTraj9a.csv")

# Show 3D space-time cube
#BSTP9a_rc <- reclassify(BSTP9a, cbind(0, NA))
#my_col <- rgb(0, 0, 0, maxColorValue = 255, alpha = 0, names = NULL)
#cubeview(BSTP9a_rc, na.color = my_col, legend = TRUE)
# Try raster to point for vox5, plot size of point by value

# Show 3D space-time cube (plotKML)
#library(plotKML)
#RBTS <- new("RasterBrickTimeSeries", rasters = BSTP9a_rc,
#            TimeSpan.begin = ymd_hms("2019/07/10 03:00:00") + minutes(seq(from = 5, to = 5*nlayers(BSTP9a_rc), by = 5)),
#            TimeSpan.end = ymd_hms("2019/07/10 03:00:00") + minutes(seq(from = 5, to = 5*nlayers(BSTP9a_rc), by = 5)))
#plotKML(RBTS)

# Show 3D space-time cube (manually)?


# Show facet plot of space-time disks of 1 bird-vessel pairing
DT_vector <- seq(from = 930, to = 970, by = 5)
plot(BSTP9a[[DT_vector]], xlab = "Easting", ylab = "Northing", col = viridis(256), 
     main = Bdf9a$DateTime[DT_vector])





################## E. VOX OVERLAP-BASED INTERACTION ###############################
# Output 1: For each bird, time series of att and enc probabilities, aggregated for all vessels
# Output 2: For each bird-vessel pair, map of att and enc probabilities

# Set WD
setwd(paste(WD, "PSTP Outputs FINAL2/", sep = ""))

# Figure out overall projection and extent?

# Load vox overlap interactions summary table
Int_vox_summary <- read.csv("Summary Tables/Interactions_vox_summary_table.csv")

# Get vector of SubTrajs with interaction
Interaction_SubTrajs <- Int_vox_summary$SubTraj[Int_vox_summary$Duration_Enc_lower + Int_vox_summary$Duration_Att_lower > 0]

# Get vector of birds with interaction
Interaction_Birds <- unique(Int_vox_summary$birdID[Int_vox_summary$Duration_Enc_lower + Int_vox_summary$Duration_Att_lower > 0])

### VOB Interaction Analysis: For each Bird with interaction...
for(i in Interaction_Birds){
  
  # Find the SubTrajs with interaction for that bird
  ST_1bird <- Int_vox_summary[(Int_vox_summary$birdID == i), 1]
  Int_ST_1bird <- ST_1bird[which(ST_1bird %in% Interaction_SubTrajs)]
  
  ### For each SubTraj with interaction...
  for(j in Int_ST_1bird){
    
    # Load vox5 and vox30 B Regular Trajectory CSVs
    BRT5 <- read.csv(paste("B Subset Dataframes/B_RegTraj_vox5_SubTraj", j, ".csv", sep = ""))
    BRT30 <- read.csv(paste("B Subset Dataframes/B_RegTraj_vox30_SubTraj", j, ".csv", sep = ""))
    
    # Subset columns, add empty Attendance column to vox30 df
    BRT30 <- BRT30[, c("DateTime", "Interpolated", "prob_enc")] %>%
      mutate(DateTime = round_date(ymd_hms(DateTime), unit = "5 minutes"), prob_att = NA) %>%
      relocate(prob_att, .before = prob_enc)
    
    # Subset columns, add empty Encounter column to vox5 df
    # Round the date to the nearest 5 minutes, to ensure that subsequent SubTrajs line up
    BRT5 <- BRT5[, c("DateTime", "Interpolated", "prob_att")] %>%
      mutate(DateTime = round_date(ymd_hms(DateTime), unit = "5 minutes"), prob_enc = NA) %>%
      relocate(prob_att, .before = prob_enc)    
    
    # Rbind the dfs, with attendance df on top
    BRT_combined <- rbind(BRT5, BRT30)
    
    # Rearrange by DateTime
    BRT_combined <- BRT_combined %>% arrange(DateTime)
    
    # Merge rows with same DateTimes from different dataframes
    for(k in 1:nrow(BRT_combined)){
      if(k > 1){
        if(is.na(BRT_combined[k, "prob_att"]) & (!is.na(BRT_combined[k-1, "prob_att"]))){
          BRT_combined[k, "prob_att"] <- BRT_combined[k-1, "prob_att"]
        }
      }
      if(k < nrow(BRT_combined)){
        if(is.na(BRT_combined[k, "prob_enc"]) & (!is.na(BRT_combined[k+1, "prob_enc"]))){
          BRT_combined[k, "prob_enc"] <- BRT_combined[k+1, "prob_enc"]
        }
      }
    }
    BRT_combined <- BRT_combined[!duplicated(BRT_combined$DateTime),]
    
    # Fill in encounter probabilities across finer-resolution time steps
    BRT_combined <- BRT_combined %>% fill(prob_att, .direction = "downup") %>% fill(prob_enc, .direction = "downup")
    
    # Plot time series of probability for this SubTraj
    time_series_SubTraj <- BRT_combined %>%
      ggplot(aes(x = ymd_hms(DateTime))) +
      theme_bw() +
      geom_line(aes(y = prob_att, colour = "prob_att"), size = 1.2) +
      geom_line(aes(y = prob_enc, colour = "prob_enc"), size = 1.2) +
      scale_color_manual(name = "Interaction Type", 
                         values = c("prob_att" = "red", "prob_enc" = "blue"),
                         labels = c("Attendance", "Encounter")) +
      geom_vline(aes(xintercept = ymd_hms(DateTime)), data = filter(BRT_combined, Interpolated == "Real")) + 
      labs(x = "DateTime (UTC)", y = "Probability of Interaction", 
           title = paste("VOB Interaction: Blue-61b, Vessel Sub-Trajectory", j, sep = " "))
    print(time_series_SubTraj)
    
    # Add this dataframe to a big single-bird interaction dataframe
    if(!exists("Bdf_int")){
      
      # Initialise dataframe
      Bdf_int <- BRT_combined
      
    } else if(nrow(Bdf_int) > 0){
      
      # Join to dataframe from previous iteration
      Bdf_int <- full_join(Bdf_int, BRT_combined, by = "DateTime")
      Bdf_int[is.na(Bdf_int)] <- 0
      
      # Union formula over time
      a <- Bdf_int$prob_att.x
      b <- Bdf_int$prob_att.y
      c <- Bdf_int$prob_enc.x
      d <- Bdf_int$prob_enc.y
      Bdf_int <- Bdf_int %>%
        mutate(Interpolated = Interpolated.x, prob_att = a + b - (a*b), prob_enc = c + d - (c*d))
      
      # Delete extra columns
      Bdf_int <- Bdf_int[, c("DateTime", "Interpolated", "prob_att", "prob_enc")]
      
    }
    
    # Load vox5 and vox30 raster bricks
    BVint5 <- brick(paste("Interaction Bricks/BVinteraction_vox5_SubTraj", j, ".nc", sep = ""))
    BVint30 <- brick(paste("Interaction Bricks/BVinteraction_vox30_SubTraj", j, ".nc", sep = ""))
    
    # Initialise raster for probability surfaces (jPPA map)
    attendance_map <- BVint5[[1]]
    encounter_map <- BVint30[[1]]
    
    # Union formula over space
    for(m in 2:nlayers(BVint5)){
      attendance_map <- attendance_map + BVint5[[m]] - (attendance_map * BVint5[[m]])
    }
    for(n in 2:nlayers(BVint30)){
      encounter_map <- encounter_map + BVint30[[n]] - (encounter_map * BVint30[[n]])
    }
    
    # Assuming projections and extents are the same, make cumulative jPPA raster layers
    if(!exists("att_map_cumulative")){
      
      # Initialise raster layers
      att_map_cumulative <- attendance_map
      enc_map_cumulative <- encounter_map
      
    } else if(exists("att_map_cumulative")){
      
      # Union formula over space
      a <- att_map_cumulative
      b <- attendance_map
      c <- enc_map_cumulative
      d <- encounter_map
      att_map_cumulative <- a + b - (a*b)
      enc_map_cumulative <- c + d - (c*d)
      
    } else {print("This didn't work.")}
    
    # Remove 0s
    attendance_map[attendance_map == 0] <- NA
    encounter_map[encounter_map == 0] <- NA
    
    # Plot maps of interaction probabilities for this SubTraj
    plot(attendance_map, xlab = "Easting", ylab = "Northing", 
         main = paste("Blue-61b, Vessel Sub-Trajectory", j, "VOB Attendance", sep = " "), 
         col = rocket(256))
    plot(encounter_map, xlab = "Easting", ylab = "Northing", 
         main = paste("Blue-61b, Vessel Sub-Trajectory", j, "VOB Encounter", sep = " "), 
         col = mako(256))
    
  }
  
  # Plot time series of Cumulative Probability for single bird
  time_series <- Bdf_int %>%
    ggplot(aes(x = ymd_hms(DateTime))) +
    theme_bw() +
    geom_line(aes(y = prob_att, colour = "prob_att"), size = 1.2) +
    geom_line(aes(y = prob_enc, colour = "prob_enc"), size = 1.2) +
    scale_color_manual(name = "Interaction Type", 
                       values = c("prob_att" = "red", "prob_enc" = "blue"),
                       labels = c("Attendance", "Encounter")) +
    geom_vline(aes(xintercept = ymd_hms(DateTime)), data = filter(Bdf_int, Interpolated == "Real")) + 
    labs(x = "DateTime (UTC)", y = "Probability of Interaction", title = "VOB Interaction: Cumulative")
  print(time_series)
  
  # Get trajectories of vessels that interact with this bird?
  #for(p in Int_ST_1bird){}
  
  # Assuming constant projection and extent, plot maps of Cumulative Probability for single bird
  att_map_cumulative[att_map_cumulative == 0] <- NA
  enc_map_cumulative[enc_map_cumulative == 0] <- NA
  plot(att_map_cumulative, xlab = "Easting", ylab = "Northing", main = "Cumulative", 
       col = rocket(256))
  plot(enc_map_cumulative, xlab = "Easting", ylab = "Northing", main = "Cumulative", 
       col = mako(256))
  
}






############# F. DIST THRESHOLD-BASED INTERACTION ################################

# Set WD to Summary Tables folder
setwd(paste(WD, "PSTP Outputs FINAL2/", sep = ""))

# Set probability threshold (above which an interaction event is considered to occur)
p0_interaction <- 0.5

### If necessary, remake distance interactions summary table
# Requires Int_vox_summary CSV and vector of Birds with interaction found from vox overlap (above)

# For each SubTraj with interaction
for(i in Interaction_Birds){
  
  # Find the SubTrajs with interaction for that bird
  ST_1bird <- Int_vox_summary[(Int_vox_summary$birdID == i), 1]
  Int_ST_1bird <- ST_1bird[which(ST_1bird %in% Interaction_SubTrajs)]
  
  ### For each SubTraj with interaction...
  for(j in Int_ST_1bird){
    
    # Load vox1 B Regular Trajectory CSV
    BRT1 <- read.csv(paste("B Subset Dataframes/B_RegTraj_vox1_SubTraj", j, "a.csv", sep = ""))
    
    # Calculate expected duration of interactions
    exp_dur_attendance <- sum(BRT1$prob_att >= p0_interaction)
    exp_dur_attendance_lower <- sum(BRT1$prob_att >= p0_lower)
    exp_dur_attendance_upper <- sum(BRT1$prob_att >= p0_upper)
    exp_dur_encounter <- sum(BRT1$prob_enc >= p0_interaction)
    exp_dur_encounter_lower <- sum(BRT1$prob_enc >= p0_lower)
    exp_dur_encounter_upper <- sum(BRT1$prob_enc >= p0_upper)
    
    # If a summary table already exists...
    if(file.exists("Summary Tables/Interactions_dist_summary_table.csv")){
      
      # Read CSV
      Interactions_dist_summary <- read.csv("Summary Tables/Interactions_dist_summary_table.csv")
      
      # Add row to table
      Interactions_dist_summary <- Interactions_dist_summary %>% 
        add_row(SubTraj = j, 
                vesselID = Int_vox_summary$vesselID[Int_vox_summary$SubTraj == j][1], 
                birdID = Int_vox_summary$birdID[Int_vox_summary$SubTraj == j][1],
                Duration_Encounter = exp_dur_encounter,
                Duration_Enc_lower = exp_dur_encounter_lower,
                Duration_Enc_upper = exp_dur_encounter_upper,
                Duration_Attendance = exp_dur_attendance,
                Duration_Att_lower = exp_dur_attendance_lower,
                Duration_Att_upper = exp_dur_attendance_upper,
                Proportion_Attendance = exp_dur_attendance / exp_dur_encounter)
      
      # Overwrite CSV
      write.csv(Interactions_dist_summary, "Summary Tables/Interactions_dist_summary_table.csv", row.names = FALSE)
      
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
                                              Proportion_Attendance = double())
      
      # Add row to table
      Interactions_dist_summary <- Interactions_dist_summary %>% 
        add_row(SubTraj = j, 
                vesselID = Int_vox_summary$vesselID[Int_vox_summary$SubTraj == j][1], 
                birdID = Int_vox_summary$birdID[Int_vox_summary$SubTraj == j][1],
                Duration_Encounter = exp_dur_encounter,
                Duration_Enc_lower = exp_dur_encounter_lower,
                Duration_Enc_upper = exp_dur_encounter_upper,
                Duration_Attendance = exp_dur_attendance,
                Duration_Att_lower = exp_dur_attendance_lower,
                Duration_Att_upper = exp_dur_attendance_upper,
                Proportion_Attendance = exp_dur_attendance / exp_dur_encounter)
      
      # Write to CSV
      write.csv(Interactions_dist_summary, "Summary Tables/Interactions_dist_summary_table.csv", row.names = FALSE)
      
    }
    
  }
  
}


### DTB Interaction Analysis

# Read in CSV of dist summary table
Int_dist_summary <- read.csv("Summary Tables/Interactions_dist_summary_table.csv")

#For each Bird with interaction...
for(i in Interaction_Birds){
  
  # Find the SubTrajs with interaction for that bird
  ST_1bird <- Int_dist_summary[(Int_dist_summary$birdID == i), 1]
  Int_ST_1bird <- ST_1bird[which(ST_1bird %in% Interaction_SubTrajs)]
  
  ### For each SubTraj with interaction...
  for(j in Int_ST_1bird){
    
    # Load vox1 B Regular Trajectory CSV
    BRT1 <- read.csv(paste("B Subset Dataframes/B_RegTraj_vox1_SubTraj", j, "a.csv", sep = ""))
    
    # Round the DateTime to the nearest minute
    BRT1 <- BRT1 %>% mutate(DateTime = round_date(ymd_hms(DateTime), unit = "minute"))
    
    # Plot time series of probability for this SubTraj
    time_series_SubTraj <- BRT1 %>%
      ggplot(aes(x = ymd_hms(DateTime))) +
      theme_bw() +
      geom_line(aes(y = prob_att, colour = "prob_att"), linewidth = 1.2) +
      geom_line(aes(y = prob_enc, colour = "prob_enc"), linewidth = 1.2) +
      scale_color_manual(name = "Interaction Type", 
                         values = c("prob_att" = "red", "prob_enc" = "blue"),
                         labels = c("Attendance", "Encounter")) +
      geom_vline(aes(xintercept = ymd_hms(DateTime)), data = filter(BRT1, Interpolated == "Real")) + 
      labs(x = "DateTime (UTC)", y = "Probability of Interaction", 
           title = paste("DTB Interaction: Blue-61b, Vessel Sub-Trajectory", j, sep = " "))
    print(time_series_SubTraj)
    
    # Add this dataframe to a big single-bird interaction dataframe
    if(!exists("Bdf_dist_int")){
      
      # Initialise dataframe
      Bdf_dist_int <- BRT1
      
    } else if(nrow(Bdf_dist_int) > 0){
      
      # Join to dataframe from previous iteration
      Bdf_dist_int <- full_join(Bdf_dist_int, BRT1, by = "DateTime")
      Bdf_dist_int[is.na(Bdf_dist_int)] <- 0
      
      # Union formula over time
      a <- Bdf_dist_int$prob_att.x
      b <- Bdf_dist_int$prob_att.y
      c <- Bdf_dist_int$prob_enc.x
      d <- Bdf_dist_int$prob_enc.y
      Bdf_dist_int <- Bdf_dist_int %>%
        mutate(Interpolated = Interpolated.x, prob_att = a + b - (a*b), prob_enc = c + d - (c*d))
      
      # Delete extra columns
      Bdf_dist_int <- Bdf_dist_int[, c("DateTime", "Interpolated", "prob_att", "prob_enc")]
      
    }
    
    # Load att and enc raster bricks
    BVdist5 <- brick(paste("Interaction Bricks/BVinteraction_dist5_SubTraj", j, "a.nc", sep = ""))
    BVdist30 <- brick(paste("Interaction Bricks/BVinteraction_dist30_SubTraj", j, "a.nc", sep = ""))
    
    # Initialise raster for probability surfaces (jPPA map)
    attendance_distmap <- BVdist5[[1]]
    encounter_distmap <- BVdist30[[1]]
    
    # Fill in NAs with zeros
    attendance_distmap[is.na(attendance_distmap)] <- 0
    encounter_distmap[is.na(encounter_distmap)] <- 0
    
    # Union formula over space
    for(m in 2:nlayers(BVdist5)){
      BVdist5[[m]][is.na(BVdist5[[m]])] <- 0
      attendance_distmap <- attendance_distmap + BVdist5[[m]] - (attendance_distmap * BVdist5[[m]])
    }
    for(n in 2:nlayers(BVdist30)){
      BVdist30[[n]][is.na(BVdist30[[n]])] <- 0
      encounter_distmap <- encounter_distmap + BVdist30[[n]] - (encounter_distmap * BVdist30[[n]])
    }
    
    # Assuming projections and extents are the same, make cumulative jPPA raster layers
    # In this case, not constant projection and extent. Can correct for this in the future.
    #if(!exists("att_distmap_cumulative")){
    #  
    #  # Initialise raster layers
    #  att_distmap_cumulative <- attendance_distmap
    #  enc_distmap_cumulative <- encounter_distmap
    #  
    #} else if(exists("att_map_cumulative")){
    #  
    #  # Union formula over space
    #  a <- att_distmap_cumulative
    #  b <- attendance_distmap
    #  c <- enc_distmap_cumulative
    #  d <- encounter_distmap
    #  att_distmap_cumulative <- a + b - (a*b)
    #  enc_distmap_cumulative <- c + d - (c*d)
    #  
    #} else {print("This didn't work.")}
    
    # Remove 0s
    attendance_distmap[attendance_distmap == 0] <- NA
    encounter_distmap[encounter_distmap == 0] <- NA
    
    # Plot maps of interaction probabilities for this SubTraj
    plot(attendance_distmap, xlab = "Easting", ylab = "Northing", 
         main = paste("Blue-61b, Vessel Sub-Trajectory", j, "DTB Attendance", sep = " "), 
         col = rocket(256))
    plot(encounter_distmap, xlab = "Easting", ylab = "Northing", 
         main = paste("Blue-61b, Vessel Sub-Trajectory", j, "DTB Encounter", sep = " "), 
         col = mako(256))
    
  }
  
  # Plot time series of Cumulative Probability for single bird
  time_series <- Bdf_dist_int %>%
    ggplot(aes(x = ymd_hms(DateTime))) +
    theme_bw() +
    geom_line(aes(y = prob_att, colour = "prob_att"), linewidth = 1.2) +
    geom_line(aes(y = prob_enc, colour = "prob_enc"), linewidth = 1.2) +
    scale_color_manual(name = "Interaction Type", 
                       values = c("prob_att" = "red", "prob_enc" = "blue"),
                       labels = c("Attendance", "Encounter")) +
    geom_vline(aes(xintercept = ymd_hms(DateTime)), data = filter(Bdf_dist_int, Interpolated == "Real")) + 
    labs(x = "DateTime (UTC)", y = "Probability of Interaction", title = "DTB Interaction: Cumulative")
  print(time_series)
  
  # Get trajectories of vessels that interact with this bird?
  #for(p in Int_ST_1bird){}
  
  # Assuming constant projection and extent, plot maps of Cumulative Probability for single bird
  # In this case, not constant projection and extent. Can correct for this in the future.
  #att_distmap_cumulative[att_distmap_cumulative == 0] <- NA
  #enc_distmap_cumulative[enc_distmap_cumulative == 0] <- NA
  #plot(att_distmap_cumulative, xlab = "Easting", ylab = "Northing", main = "DTB Attendance: Cumulative", 
  #     col = rocket(256))
  #plot(enc_distmap_cumulative, xlab = "Easting", ylab = "Northing", main = "DTB Encounter: Cumulative", 
  #     col = mako(256))
  
}







############ G. DURATION OF INTERACTION #################

# Set WD
setwd(paste(WD, "PSTP Outputs FINAL2/", sep = ""))

# Read CSVs
Interactions_vox_summary <- read.csv("Summary Tables/Interactions_vox_summary_table.csv")
Interactions_dist_summary <- read.csv("Summary Tables/Interactions_dist_summary_table.csv")

# TO REMOVE LATER: Add upper and lower probability thresholds to summary df
Interactions_vox_summary <- Interactions_vox_summary %>%
  mutate(Duration_Enc_lower = 0, Duration_Enc_upper = 0, Duration_Att_lower = 0, Duration_Att_upper = 0)
Interactions_vox_summary$Duration_Enc_lower[Interactions_vox_summary$SubTraj == 9] <- 600
Interactions_vox_summary$Duration_Enc_upper[Interactions_vox_summary$SubTraj == 9] <- 400
Interactions_vox_summary$Duration_Att_lower[Interactions_vox_summary$SubTraj == 9] <- 200
Interactions_vox_summary$Duration_Att_upper[Interactions_vox_summary$SubTraj == 9] <- 100
Interactions_vox_summary$Duration_Enc_lower[Interactions_vox_summary$SubTraj == 14] <- 900
Interactions_vox_summary$Duration_Enc_upper[Interactions_vox_summary$SubTraj == 14] <- 600
Interactions_vox_summary$Duration_Att_lower[Interactions_vox_summary$SubTraj == 14] <- 400
Interactions_vox_summary$Duration_Att_upper[Interactions_vox_summary$SubTraj == 14] <- 200


# Pivot the dataframes to get Int_type and Duration columns, add Method column
Int_vox_2 <- Interactions_vox_summary[!is.na(Interactions_vox_summary$Proportion_Attendance),] %>%
  rename(Encounter = Duration_Encounter, Attendance = Duration_Attendance) %>%
  pivot_longer(cols = c("Encounter", "Attendance"), 
               names_to = "Int_type", values_to = "Duration") %>%
  mutate(Method = "VOB", 
         Duration_lower = ifelse(Int_type == "Encounter", Duration_Enc_lower, Duration_Att_lower), 
         Duration_upper = ifelse(Int_type == "Encounter", Duration_Enc_upper, Duration_Att_upper)) %>%
  dplyr::select(-c(CompTime, Duration_Enc_lower, Duration_Enc_upper, Duration_Att_lower, Duration_Att_upper))
Int_dist_2 <- Interactions_dist_summary[!is.na(Interactions_dist_summary$Proportion_Attendance),] %>%
  rename(Encounter = Duration_Encounter, Attendance = Duration_Attendance) %>%
  pivot_longer(cols = c("Encounter", "Attendance"), 
               names_to = "Int_type", values_to = "Duration") %>%
  mutate(Method = "DTB", 
         Duration_lower = ifelse(Int_type == "Encounter", Duration_Enc_lower, Duration_Att_lower), 
         Duration_upper = ifelse(Int_type == "Encounter", Duration_Enc_upper, Duration_Att_upper)) %>%
  dplyr::select(-c(CompTime, Duration_Enc_lower, Duration_Enc_upper, Duration_Att_lower, Duration_Att_upper))

# Rbind the dataframes
Int_both_summary <- rbind(Int_vox_2, Int_dist_2)

# ORIGINAL: Plot a double bar chart
duration_chart <- Int_both_summary %>% ggplot(aes(x = Method, y = Duration, fill = Int_type)) +
  geom_bar(position = "identity", stat = "identity") +
  scale_fill_manual(values = c("#F8766D", "#619CFF")) +
  geom_text(aes(label = Duration), vjust = 1.5, colour = "black", size = 3.5) +
  facet_wrap( ~ SubTraj) +
  theme_bw() +
  labs(x = "Analytical Method", y = "Duration (minutes)") + 
  guides(fill = guide_legend(title = "Interaction Type"))
duration_chart

# NEW: Plot a bar chart with error brackets
# Make double bar chart, add error bars, vessel not SubTraj, hours not minutes
Int_both_summary %>% ggplot(aes(x = Method, y = Duration/60, fill = Int_type)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(ymin = (Duration_upper)/60, ymax = (Duration_lower)/60), width = .3, linewidth = .6, position = position_dodge(.9)) +
  scale_fill_manual(values = c("#F8766D", "#619CFF")) +
  #geom_text(aes(label = round(Duration/60, digits = 1)), position = position_dodge(.9), vjust = -0.6, hjust = -0.5, colour = "black", size = 3.5) +
  facet_wrap( ~ SubTraj) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(linewidth = .1, color = "grey95"),
        plot.title = element_text(size = 11, hjust = 0.5)) +
  labs(x = "Analytical Method", y = "Duration (hours)", title = "Vessel Sub-Trajectory") + 
  guides(fill = guide_legend(title = "Interaction Type"))



