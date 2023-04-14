####################################################################
#################### 2. CLEANING GFW DATA ##########################
####################################################################

# Takes in a CSV of raw AIS vessel tracks
# This case study obtained data from Global Fishing Watch through BirdLife International

####################### CONTENTS ###########################

# A. SETUP
# B. CLEAN DATA AND MAKE NEW COLUMNS
# C. INSPECT DATA


################ A. SETUP #######################

library(tidyverse)
library(lubridate)
library(data.table)


# Set WD to R Scripts and Data folder
setwd("__________")

# Function to calculate Distance
distance.lat.long <- function(x1, y1, x2, y2){
  # A function to return the distance in KMs separating two points on the Earth's surface
  # Uses the simple spherical model
  # Points MUST be decimalised degrees
  radians <- function(xin){
    tmp.rad <- (pi*xin)/180
    return(tmp.rad)
  } 
  r <- 6378.388
  tmp.trig = (sin(radians(y1))*sin(radians(y2)))+(cos(radians(y1))*cos(radians(y2))*cos(radians(x1)-radians(x2)))
  tmp.d = r*acos(tmp.trig)
  ifelse(tmp.d > 0, return(tmp.d), return(0)) # I altered this line of code
}




############## B. CLEAN DATA AND MAKE NEW COLUMNS ################
############### JUST 5 DAY TEST (8-12 JULY 2019) ###############

# Read CSV
gfw0 <- read.csv("GFW Queried Data/GFWquery_5daytest.csv")
#View(gfw0, "GFW Data")

# Add boolean fishing activity column
gfw1 <- gfw0 %>% mutate(fishing = ifelse(hours == fishing_hours, 1, 0))
sum(gfw1$fishing)
#View(gfw1, "GFW Data")

# Correct error with 1 vessel having 2 ssvids
gfw1$flag <- ifelse(gfw1$ssvid == 98346636, "CHN", gfw1$flag)
gfw1$ssvid <- ifelse(gfw1$ssvid == 98346636, 412420195, gfw1$ssvid)

# Subset data to relevant columns, arrange by birdID/flag/ssvid/DateTime
# dt_int unit is minutes
gfw2 <- gfw1[, c("flag", "ssvid", "timestamp", "lat", "lon", "fishing", "hours", "id", "index")]
gfw2 <- gfw2 %>%
  mutate(flag = ifelse(flag == "", "Unknown", flag)) %>%
  mutate(flag = as.factor(flag), ssvid = as.factor(ssvid), id = as.factor(id)) %>%
  mutate(timestamp = ymd_hms(timestamp, tz = "UTC"), hours = hours*60) %>%
  rename(DateTime = timestamp, Lat = lat, Long = lon, birdID = id, dt_int = hours) %>%
  arrange(birdID, flag, ssvid, DateTime)
#View(gfw2, "GFW Data")

# Remove duplicates
nrow(gfw2)
gfw3 <- gfw2[!duplicated(gfw2[,-9]),]
gfw3 <- gfw3[-which(gfw3$DateTime == lag(gfw3$DateTime)),]
nrow(gfw3) # This removes a lot
#View(gfw3)

# Add columns for distance interval (km) and speed (km/h)
gfw4 <- gfw3 %>%
  mutate(dist_int = distance.lat.long(Long, Lat, lag(Long), lag(Lat))) %>%
  mutate(speed = (dist_int / dt_int) * 60)

# Set distance interval and speed interval to 0 where dataframe swaps to new vessel
gfw4$dist_int[gfw4$ssvid != lag(gfw4$ssvid)] <- 0
gfw4$dist_int[1] <- 0

# Set first row to 0 distance and speed
gfw4$speed[gfw4$ssvid != lag(gfw4$ssvid)] <- 0
gfw4$speed[1] <- 0

# INSPECT SPEEDS (see below)

# Speed filter: remove positions with speeds over 40 km/h, then recalculate dt_int, dist_int, and speed
# 7 positions removed
gfw4 <- gfw4 %>%
  filter(speed <= 40) %>%
  mutate(dt_int = abs(as.numeric(difftime(DateTime, lag(DateTime), units = "mins")))) %>%
  mutate(dist_int = distance.lat.long(Long, Lat, lag(Long), lag(Lat))) %>%
  mutate(speed = (dist_int / dt_int) * 60)
gfw4$dist_int[gfw4$ssvid != lag(gfw4$ssvid)] <- 0
gfw4$dist_int[1] <- 0
gfw4$speed[gfw4$ssvid != lag(gfw4$ssvid)] <- 0
gfw4$speed[1] <- 0
gfw4$dt_int[gfw4$ssvid != lag(gfw4$ssvid)] <- 0
gfw4$dt_int[1] <- 0

# Make new column of vesselID (ssvid, but confidential)
vesselIDNum <- 1
gfw4 <- gfw4 %>% add_column(vesselID = NA)
gfw4[1, "vesselID"] <- 1
for(i in 2:nrow(gfw4)){
  if((gfw4[i, "birdID"] != gfw4[i-1, "birdID"]) | (gfw4[i, "ssvid"] != gfw4[i-1, "ssvid"])){
    vesselIDNum <- vesselIDNum + 1
    gfw4[i, "vesselID"] <- vesselIDNum
  } else {
    gfw4[i, "vesselID"] <- vesselIDNum
  }
}


# Make new column of sub-trajectories (subset by ssvid, no more than 4 hours between positions)
SubTrajNum <- 1
gfw4 <- gfw4 %>% add_column(SubTraj = NA)
gfw4[1, "SubTraj"] <- 1
for(i in 2:nrow(gfw4)){
  if((gfw4[i, "birdID"] != gfw4[i-1, "birdID"]) | (gfw4[i, "ssvid"] != gfw4[i-1, "ssvid"]) | (gfw4[i, "dt_int"] > 240)){
    SubTrajNum <- SubTrajNum + 1
    gfw4[i, "SubTraj"] <- SubTrajNum
  } else {
    gfw4[i, "SubTraj"] <- SubTrajNum
  }
}
View(gfw4)




# Write CSV
setwd("C:/Users/jonat/Desktop/University of Oxford BCM/Oxford BCM 2021-2022/BCM Dissertation/R Scripts and Data/GFW Queried Data/")
write.csv(gfw4, "GFW_5daytest_clean.csv")
setwd("C:/Users/jonat/Desktop/University of Oxford BCM/Oxford BCM 2021-2022/BCM Dissertation/R Scripts and Data")


# Write CSV for vessel near 64334 on 10 July 2019 AM?





################### C. INSPECT DATA ###################

nrow(gfw0)
length(which(duplicated(gfw0[,"timestamp"])))

levels(as.factor(gfw0$ssvid)) # 16 vessels in this 5 day period
levels(as.factor(gfw0$id)) # Overlap with birds 64334 and 64336 in this period

# Inspect sub trajectories
count(gfw4, SubTraj)

# How many vessels per flag state
gfw2 %>% 
  group_by(flag) %>%
  summarise(vessels = n_distinct(ssvid))

# Which vessels in each flag state, and how many positions each
gfw2 %>% count(flag, ssvid)

# Fishing activity: 84% of positions are fishing
sum(gfw3$fishing)/nrow(gfw3)



# Time interval distribution

gfw4 %>%
  summarise(min = min(dt_int), max = max(dt_int), mean = mean(dt_int), sd = sd(dt_int))

gfw3 %>%
  filter(dt_int < 20) %>%
  ggplot(aes(x = dt_int, fill = ssvid)) + 
  geom_histogram(binwidth = 1) +
  facet_grid(rows = vars(flag))

gfw3 %>%
  filter(dt_int > 0) %>%
  ggplot(aes(x = DateTime, y = dt_int, color = ssvid)) +
  geom_point(size = 0.5, position = "jitter") +
  facet_grid(rows = vars(flag))


# Speed distribution

gfw4 %>%
  summarise(min = min(speed), max = max(speed), mean = mean(speed), median = median(speed))

gfw4 %>%
  filter(speed < 50) %>%
  ggplot(aes(x = speed, fill = fishing)) + 
  geom_histogram(binwidth = 1) +
  facet_grid(rows = vars(flag))

gfw4 %>%
  filter(speed < 50) %>%
  ggplot(aes(x = DateTime, y = speed, color = fishing)) +
  geom_point(size = 0.5, position = "jitter") +
  facet_grid(rows = vars(flag))

gfw4[gfw4$speed > 40,] # 40 kph seems like a good speed filter


