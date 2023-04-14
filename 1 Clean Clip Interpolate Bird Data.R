#######################################################
# 1. CLEAN CLIP INTERPOLATE BIRD DATA

# Takes in CSVs of seabird tracks 
# This case study takes in Rainier GPS Toroa tracks, already with some minor cleaning
# Exception: "Birdhour daily 100km scale overlap analysis" takes in CSVs from ArcGIS Pro Join Features analysis

###################### CONTENTS ##############################

# A. Libraries and WD
# B. Clean data and make new columns
# C. Clip Rainier data, add birdhour column
# D. Clip to 2 birds over 5 days, find dBBMM variances
# E. Interpolate trajectory
# F. Clip interpolated trajectory
# G. Birdhour daily 100km scale overlap analysis
# H. Inspect and visualise data



################ A. LIBRARIES AND WD #######################

library(tidyverse)
library(lubridate)
library(data.table)
library(move)


# Set WD to R Scripts and Data folder
setwd("________")



############## B. CLEAN DATA AND MAKE NEW COLUMNS ################

# Bind Rainier CSVs together

files <- list.files(path = "Rainier Clean", pattern = ".csv", full.names = TRUE)
data <- do.call(rbind, lapply(files, fread))
#View(data, title = "All_Rainier")


# Remove duplicates
nrow(data)
df1 <- data[!duplicated(data[,-1]),]
nrow(df1)


# Filter to May-September 2019, add datetime interval column
df2 <- df1 %>%
  mutate(ID = as.factor(ID), LocQual = as.factor(LocQual)) %>%
  mutate(DateTime = ymd_hms(DateTime)) %>%
  filter(year(DateTime) == 2019, month(DateTime) >= 5, month(DateTime) <= 9) %>%
  mutate(dt_int = as.numeric(difftime(DateTime, lag(DateTime), units = "mins"))) %>%
  filter(dt_int != "NA")


# Set time interval to 0 where dataframe swaps to new bird
df2$dt_int[df2$ID != lag(df2$ID)] <- 0


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


# Add columns for distance interval (km) and speed (km/hr)
df3 <- df2 %>%
  mutate(dist_int = distance.lat.long(Long, Lat, lag(Long), lag(Lat))) %>%
  mutate(speed = (dist_int / dt_int) * 60) %>%
  filter(dist_int != "NaN")

# Set distance interval and speed interval to 0 where dataframe swaps to new bird
df3$dist_int[df3$ID != lag(df3$ID)] <- 0
df3$speed[df3$ID != lag(df3$ID)] <- 0

# View new data frame
#View(df3, title = "All_Rainier")





########### C. CLIP RAINIER DATA, ADD BIRDHOUR COLUMN ###################

# Remove points where:
#   Both dist_int AND lead(dist_int) exceed 150km
#   AND both dt_int AND lead(dt_int) exceed 180min

df4 <- df3 %>%
  filter((dist_int <= 150) |
           (lead(dist_int) <= 150) |
           (dt_int <= 180) |
           (lead(dt_int) <= 180))

# Remove points where dt_int is less than 2 mins and dist_int is 0
df5 <- df4[!((df4$dt_int < 2) & (df4$dist_int == 0)),]

# Recalculate dt_int, dist_int, and speed
df6 <- df5 %>%
  mutate(dt_int = as.numeric(difftime(DateTime, lag(DateTime), units = "mins"))) %>%
  mutate(dist_int = distance.lat.long(Long, Lat, lag(Long), lag(Lat))) %>%
  filter(dt_int != "NA", dist_int != "NaN") %>%
  mutate(speed = (dist_int / dt_int) * 60)

# Set dt_int, dist_int, speed to 0 where dataframe swaps to new bird
df6$dt_int[df6$ID != lag(df6$ID)] <- 0
df6$dist_int[df6$ID != lag(df6$ID)] <- 0
df6$speed[df6$ID != lag(df6$ID)] <- 0

# Calculate Birdhour
# Birdhour = (Prior time interval - Next time interval)/2
# I am not weighting by cohort for this analysis
df6 <- df6 %>%
  mutate(birdhour = (dt_int + lead(dt_int))/2/60)

# Birdhour = 0 at start of track
df6$birdhour[df6$ID != lag(df6$ID)] <- 0

# Inspect
nrow(df3)
nrow(df5)
nrow(df6)
nrow(df3) - nrow(df6)

# View clipped dataframe
#View(df6, title = "All_Rainier")




########## D. CLIP TO 2 BIRDS OVER 5 DAYS, FIND dBBMM VARIANCES ##################

# Set parameters for dynamic Brownian motion variance (Kranstauber et al. 2012)
dBBMM_window <- 7 # has to be odd, unit is locations, 40 * 7 = 4 hours 40 mins
dBBMM_margin <- 3 # Minimum 3 locations
sigma2_2 <- 12 # 12m telemetry error for Rainier GPS fixes
WGS1984 <- CRS("+init=epsg:4326") # Set CRS of input data: WGS 1984

# Clip df to 8-12 July 2019, split to separate dfs for 64334 and 64336
df7_64334 <- df6 %>%
  filter(month(DateTime) == 7, day(DateTime) >= 8, day(DateTime) <= 12) %>%
  filter((ID == 64334)) %>%
  mutate(sigma1_2 = NA)
df7_64336 <- df6 %>%
  filter(month(DateTime) == 7, day(DateTime) >= 8, day(DateTime) <= 12) %>%
  filter((ID == 64336)) %>%
  mutate(sigma1_2 = NA)

# Convert each to move object
df8_64334 <- move(x = df7_64334$Long,
                  y = df7_64334$Lat,
                  time = df7_64334$DateTime,
                  proj = WGS1984,
                  animal = df7_64334$ID,
                  data = df7_64334,
                  removeDuplicatedTimestamps = TRUE)
df8_64336 <- move(x = df7_64336$Long,
                  y = df7_64336$Lat,
                  time = df7_64336$DateTime,
                  proj = WGS1984,
                  animal = df7_64336$ID,
                  data = df7_64336,
                  removeDuplicatedTimestamps = TRUE)

# For each bird, transform coordinates to Azimuthal Equi-Distance projection centred at center of track (default)
df8_64334 <- spTransform(df8_64334, center = TRUE)
df8_64336 <- spTransform(df8_64336, center = TRUE)

# For each bird, create dBM variance S4 object
dBBMMvar_64334 <- brownian.motion.variance.dyn(object = df8_64334,
                                               location.error = sigma2_2,
                                               window.size = dBBMM_window,
                                               margin = dBBMM_margin)
dBBMMvar_64336 <- brownian.motion.variance.dyn(object = df8_64336,
                                               location.error = sigma2_2,
                                               window.size = dBBMM_window,
                                               margin = dBBMM_margin)


# For each bird, add variance column (default time units minutes, so convert to seconds)
vector_variances_64334 <- getMotionVariance(dBBMMvar_64334)
df7_64334$sigma1_2 <- vector_variances_64334/60
vector_variances_64336 <- getMotionVariance(dBBMMvar_64336)
df7_64336$sigma1_2 <- vector_variances_64336/60

# For each bird, fill in the variances at the margins
# (Introduces some error, but necessary for the whole trajectory to be usable)
df7_64334 <- df7_64334 %>%
  fill(sigma1_2, .direction = "down") %>% # Fill in variances for upcoming interpolated positions
  fill(sigma1_2, .direction = "up") # Then fill those that were missed at start of trajectory
df7_64336 <- df7_64336 %>%
  fill(sigma1_2, .direction = "down") %>% # Fill in variances for upcoming interpolated positions
  fill(sigma1_2, .direction = "up") # Then fill those that were missed at start of trajectory


# Rbind dfs
df9 <- rbind(df7_64334, df7_64336)

# Output to CSV
setwd("C:/Users/jonat/Desktop/University of Oxford BCM/Oxford BCM 2021-2022/BCM Dissertation/R Scripts and Data/Rainier Clean Clipped/")
write.csv(df9, "Rainier_Clip_Dissertation.csv")
setwd("C:/Users/jonat/Desktop/University of Oxford BCM/Oxford BCM 2021-2022/BCM Dissertation/R Scripts and Data")






################### E. INTERPOLATE TRAJECTORY ########################

# Inspect df6
df6

# CRS code: WGS1984
WGS1984 <- CRS("+init=epsg:4326")


# Function to interpolate a data.frame from above
# Arguments: df = pre-formatted non-interpolated data.frame; interp.time = interpolation time step in minutes
# Returns a data.frame with interpolated positions

interpolate.df <- function(df, interp.time){
  
  # Set input CRS
  WGS1984 <- CRS("+init=epsg:4326")
  
  # Create move object for bird
  mv <- move(x = df$Long,
             y = df$Lat,
             time = df$DateTime,
             proj = WGS1984,
             animal = df$ID,
             data = df,
             removeDuplicatedTimestamps = TRUE)
  
  
  
  # Interpolate trajectory: 1 position every 10 minutes
  mv1 <- interpolateTime(mv,
                         time = as.difftime(interp.time, units = "mins"),
                         spaceMethod = "greatcircle")
  
  # Transform into data.frame and return
  mvdf1 <- as.data.frame(mv1)
  
  # Subset to columns: coords.x1, coords.x2, timestamps, ID, Meta, locQual, sensor
  mvdf2 <- mvdf1[,c("coords.x1", "coords.x2", "timestamps", "ID", "Meta", "LocQual", "dt_int", "dist_int", "speed", "birdhour", "sensor")]
  
  # If sensor = "interpolateTime", change to "Interpolated", change locQual to "I"
  # Then change names of columns
  mvdf3 <- mvdf2
  mvdf3$LocQual <- as.character(mvdf3$LocQual)
  mvdf3$sensor <- as.character(mvdf3$sensor)
  for(i in 1:nrow(mvdf3)){
    if(mvdf3[i,"sensor"] == "interpolateTime"){
      mvdf3[i,"LocQual"] <- "I"
      mvdf3[i,"sensor"] <- "Interpolated"
    } else {
      mvdf3[i,"sensor"] <- "Real"
    }
  }
  mvdf3$LocQual <- as.factor(mvdf3$LocQual)
  mvdf3$sensor <- as.factor(mvdf3$sensor)
  mvdf3 <- mvdf3 %>% rename(Long = coords.x1, Lat = coords.x2, DateTime = timestamps, Interpolated = sensor)
  
  # Add in missing real positions: make data.frame of reals that resembles interpolated
  mvdf4 <- df[,c("Long", "Lat", "DateTime", "ID", "Meta", "LocQual", "dt_int", "dist_int", "speed", "birdhour")]
  mvdf4 <- mvdf4 %>% mutate(Interpolated = "Real")
  
  # Rbind
  mvdf5 <- rbind(mvdf3, mvdf4)
  
  # Reorder by DateTime, remove duplicate timestamps
  mvdf5 <- mvdf5[order(mvdf5[,"DateTime"]),]
  mvdf6 <- mvdf5[!duplicated(mvdf5$DateTime),]
  
  # Return interpolated data.frame
  return(mvdf6)
}


# Running interpolation function on each individual bird
# because for some reason "move" won't let me do them all at once
# and "moveStack" isn't working either

R1 <- df6 %>% filter(ID == 64331) %>% interpolate.df(interp.time = 10)
R2 <- df6 %>% filter(ID == 64334) %>% interpolate.df(interp.time = 10)
R3 <- df6 %>% filter(ID == 64335) %>% interpolate.df(interp.time = 10)
R4 <- df6 %>% filter(ID == 64336) %>% interpolate.df(interp.time = 10)
R5 <- df6 %>% filter(ID == 84933) %>% interpolate.df(interp.time = 10)
R6 <- df6 %>% filter(ID == 84934) %>% interpolate.df(interp.time = 10)

# Rbind all interpolated tracks into one data.frame
Int1 <- rbind(R1, R2, R3, R4, R5, R6)
View(Int1, title = "Rainier_Interpolated1")





################# F. CLIP INTERPOLATED TRAJECTORY ################

# To remove Interpolated positions that are between spatiotemporally far-apart Real positions
# Thresholds from Real positions: 150km, 3 hours (must exceed both thresholds)


# Mutate a new column "new_burst" that = TRUE if the previous Real was >150km away AND >3hr prior
Int2 <- Int1 %>% mutate(new_burst = FALSE)
Int2$new_burst <- ifelse((Int2$dist_int > 150) & (Int2$dt_int > 180), TRUE, FALSE)
Int2$new_burst <- replace_na(Int2$new_burst, FALSE)

# Checks that resulting dataframe makes sense
sum(Int2$new_burst == TRUE, na.rm = TRUE)
#which((Int2$dt_int > 180) & (Int2$dist_int > 150))
#which(Int2$new_burst == TRUE)
#Int2[275:285,]

# Mutate a new column "num_interp" that shows number of Interpolated before each Real
Int3 <- Int2 %>%
  group_by(ID) %>%
  group_by(temporary1 = cumsum(Int2$Interpolated == "Real") - (Int2$Interpolated == "Real")) %>%
  mutate(num_interp = ifelse(Interpolated == "Real", sum(Interpolated == "Interpolated"), 0)) %>%
  ungroup()

# If new_burst == TRUE, delete num_interp number of rows above it
Int4 <- Int3 %>% mutate(temporary2 = 0)
for(i in 1:nrow(Int4)){
  if(Int4$new_burst[i] == TRUE){
    Int4[((i-Int4$num_interp[i]):(i-1)),"temporary2"] <- 1
  }
}
Int5 <- Int4[!(Int4$temporary2 == 1),]

# Double check this worked
length(which(Int4$temporary2 == 1)) # Should be a big number, I got 25766
length(which(Int5$temporary2 == 1)) # Should be 0

# Subset to remove 1st, temporary, and "num_interp" columns, add an Index ID just in case, reorder columns
Int6 <- subset(Int5, select = -c(temporary1, temporary2, num_interp))
Int6 <- Int6 %>%
  mutate(Index = row_number()) %>%
  relocate(Index)
nrow(Int6)
View(Int6)



# Export as CSV, then reset working directory
setwd("C:/Users/jonat/Desktop/University of Oxford BCM/Oxford BCM 2021-2022/BCM Dissertation/R Scripts and Data/Interpolated Data/")
write.csv(Int6, "Rainier_Interpolated2.csv")
setwd("C:/Users/jonat/Desktop/University of Oxford BCM/Oxford BCM 2021-2022/BCM Dissertation/R Scripts and Data")


# Export as CSV: initial GFW Query (2019 July 8-12)
# Then take into GIS and clip it into study area
# (Eventually I should try to have code for this in R, if time)
Int7 <- Int6 %>% filter(month(DateTime) == 7, day(DateTime) >= 8, day(DateTime) <= 12)
nrow(Int7)
setwd("C:/Users/jonat/Desktop/University of Oxford BCM/Oxford BCM 2021-2022/BCM Dissertation/R Scripts and Data/Interpolated Data/")
write.csv(Int7, "Rainier_Interpolated2_8-12.7.csv")
setwd("C:/Users/jonat/Desktop/University of Oxford BCM/Oxford BCM 2021-2022/BCM Dissertation/R Scripts and Data")


# Export as CSV: for PSTP testing (just bird 64334, July 10)
Int8 <- Int6 %>% filter(ID == 64334, month(DateTime) == 7, day(DateTime) == 10)
setwd("C:/Users/jonat/Desktop/University of Oxford BCM/Oxford BCM 2021-2022/BCM Dissertation/R Scripts and Data/Interpolated Data/")
write.csv(Int8, "Rainier_Interp_64334_10.7.csv")
setwd("C:/Users/jonat/Desktop/University of Oxford BCM/Oxford BCM 2021-2022/BCM Dissertation/R Scripts and Data")






############# G. BIRDHOUR DAILY 100km SCALE OVERLAP ANALYSIS #################

# See Birdhour calculation above
# Used Join Features tool in ArcGIS Pro
# Was forced to do it month by month for the tool to run
# So I ended up with 5 CSVs, one for each month May-Sep 2019

library(tidyverse)
library(lubridate)
library(data.table)

# Bind Bird-Vessel Overlap CSVs together
setwd("C:/Users/jonat/Desktop/University of Oxford BCM/Oxford BCM 2021-2022/BCM Dissertation/R Scripts and Data")
files <- list.files(path = "Daily Overlap Scale Data", pattern = ".csv", full.names = TRUE)
BVOL1 <- do.call(rbind, lapply(files, fread))

# Subset to relevant columns, clean up
BVOL2 <- BVOL1[,c("DateTime", "Date", "Lat", "Long", "ID", "Meta", "LocQual", "dt_int", "dist_int", "speed", "birdhour", "SUM_Apparent_Fishing_hours", "COUNT")]
BVOL2 <- BVOL2 %>%
  rename(DailyEffort = SUM_Apparent_Fishing_hours, Num_GFWCells = COUNT) %>%
  mutate(ID = as.factor(ID), DateTime = ymd_hm(DateTime), Date = ymd(Date)) %>%  #6 DateTime cells failed to parse but doesn't matter for this
  arrange(Date, DateTime)

# Set NAs to 0
BVOL2$DailyEffort[is.na(BVOL2$DailyEffort)] <- 0

# Add column for hourly effort within 100km of each location (daily effort in hours / 24)
# Add column for overlap = hourly fishing effort * unweighted birdhour
BVOL3 <- BVOL2 %>%
  mutate(HourlyEffort = DailyEffort/24) %>%
  mutate(Overlap = HourlyEffort * birdhour)

# View Overlap dataframe. Inspection code below.
View(BVOL3, title = "Daily Overlap")






############## H. INSPECT AND VISUALISE DATA ####################

#Testing distance.lat.long function
x1 <- -171.8932
y1 <- -41.5735
x2 <- -171.8988
y2 <- -41.5751
distance.lat.long(x1, y1, x2, y2)
distance.lat.long(x2, y2, x1, y1)
distance.lat.long(y1, x1, y2, x2)
distance.lat.long(y2, x2, y1, x1)
distance.lat.long(y1, x2, y2, x1)


# Inspect time intervals

df9 %>%
  summarize(min = min(dt_int), max = max(dt_int), mean = mean(dt_int), sd = sd(dt_int))

df9 %>%
  filter(dt_int < 200, dt_int > 0) %>%
  ggplot(aes(x = dt_int, fill = ID)) + 
  geom_histogram(binwidth = 5) +
  facet_grid(rows = vars(ID))

df2 %>%
  filter(dt_int < 200, dt_int > 0) %>%
  ggplot(aes(x = DateTime1, y = dt_int, color = ID)) +
  geom_point(size = 0.5, position = "jitter") +
  facet_grid(rows = vars(ID))


# Inspect distance and speed
df3 %>%
  summarize(min = min(dist_int), max = max(dist_int), mean = mean(dist_int), median = median(dist_int))

df3 %>%
  summarize(min = min(speed), max = max(speed), mean = mean(speed), median = median(speed))  

df3 %>%
  filter(dist_int < 100, dt_int > 0) %>%
  ggplot(aes(x = dist_int, fill = ID)) + 
  geom_histogram(binwidth = 1) +
  facet_grid(rows = vars(ID))

df3 %>%
  filter(dist_int < 500, dt_int > 0) %>%
  ggplot(aes(x = DateTime, y = dist_int, color = ID)) +
  geom_point(size = 0.5, position = "jitter") +
  facet_grid(rows = vars(ID))

df3 %>%
  filter(speed < 50, dt_int > 0) %>%
  ggplot(aes(x = speed, fill = ID)) + 
  geom_histogram(binwidth = 1) +
  facet_grid(rows = vars(ID))

df3 %>%
  filter(speed < 100, dt_int > 0) %>%
  ggplot(aes(x = DateTime, y = speed, color = ID)) +
  geom_point(size = 0.5, position = "jitter") +
  facet_grid(rows = vars(ID))



# Inspect birdhour

df6 %>%
  group_by(ID) %>%
  summarize(min = min(birdhour), max = max(birdhour), mean = mean(birdhour), median = median(birdhour))



# Inspect new clipped dataframe

df6 %>%
  group_by(ID) %>%
  summarize(min = min(dt_int), max = max(dt_int), mean = mean(dt_int), median = median(dt_int))

df6 %>%
  group_by(ID) %>%
  summarize(min = min(dist_int), max = max(dist_int), mean = mean(dist_int), median = median(dist_int))

df6 %>%
  group_by(ID) %>%
  summarize(min = min(speed), max = max(speed), mean = mean(speed), median = median(speed))  



# Inspect interpolated trajectories

View(R1)
nrow(R1)
Rainier_Interpolated1
View(Rainier_Interpolated1)
nrow(Rainier_Interpolated1)
nrow(df6)
nrow(Rainier_Interpolated1) - nrow(df6)


# Figure out Argos vs GPS

which(Int6$LocQual != "G" & Int6$LocQual != "I")
Int6[4290:4305,]
Int6[7180:7195,]


# Inspect spatial and temporal conditions for clipping real and interpolated points

Int6[Int6$dt_int>180 & Int6$dist_int<150 & !is.na(Int6$dt_int),] #510 points over 3 hours away but within 150km
Int6[Int6$dt_int<180 & Int6$dist_int>150 & !is.na(Int6$dt_int),] #2 points within 3 hours but over 150km

nrow(df4[df4$dt_int>180 & df4$dist_int<150 & !is.na(df4$dt_int),]) #510 points before initial clipping as well
nrow(df4[df4$dt_int<180 & df4$dist_int>150 & !is.na(df4$dt_int),]) #2 points before initial clipping as well


# Inspect Bird-Vessel Overlap (BVOL) at daily 100km scale, based on unweighted birdhours

BVOL3 %>%
  group_by(ID) %>%
  summarize(min = min(Overlap), max = max(Overlap), mean = mean(Overlap), median = median(Overlap))

BVOL3 %>%
  group_by(month(Date)) %>%
  summarize(min = min(Overlap), max = max(Overlap), mean = mean(Overlap), median = median(Overlap))

BVOL3 %>%
  filter(Overlap > 0) %>%
  group_by(Date) %>%
  summarize(Total_Overlap = sum(Overlap)) %>%
  arrange(desc(Total_Overlap)) %>%
  print.data.frame()

BVOL3 %>%
  ggplot(aes(x = Date, y = Overlap, fill = as.factor(month(Date)))) + 
  geom_col()




