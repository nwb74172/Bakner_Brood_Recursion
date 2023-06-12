# Batch Recurse Code #
#-------------------------------------
#load libraries
#-------------------------------------
library(lubridate)
library(recurse)
library(move)
library(terra)
library(raster)
library(scales)
library(viridis)
library(lubridate)
library(reshape2)
library(maptools)
library(cluster)
library(sf)
library(tidyr)
library(data.table)
library(dplyr)
library(plyr)
library(readr)
library(data.table)

# path to folder that holds multiple .csv files (make sure to leave / at end)
# must be ran for each study site separately
folder <- "E:/Brood_Recursion/Raw_Data_Files/Raw_Minitrack_GPS_Data/WEBB/"  

# Create list of all .csv files in folder
file_list <- list.files(path=folder, pattern="*.csv")

# Begin loop
for (i in 1:length(file_list)){
  #read in .csv file
  df <- read.csv(paste(folder, file_list[i], sep=''))
  
  head(df)
  # make lat and lon numeric
  df$Lat <- as.numeric(df$Latitude)
  df$Lon <- as.numeric(df$Longitude)
  
  # Omit N/A from WEBB, SC data
  df <- na.omit(df)
  
  # Convert to POSIXct
  df$datetime <- as.POSIXct(strptime(paste(df$LMT_Date, df$LMT_Time),"%d.%m.%Y %H:%M"))
  xy <- df[,c(5,1)]
  
  # Creates spatial point for UTM conversion
  xy2 <- df[,c(12,11)]
  cord.dec <- SpatialPoints(cbind(xy2$Lon, xy2$Lat), proj4string = CRS("+proj=longlat"))
  
  # Transforming coordinate to UTM using EPSG=32615 for WGS=84, UTM Zone=15N)
  #Louisiana
  #cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=15N ellps=WGS84"))
  #SRS, CC, BFG, WEBB
  cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=17N ellps=WGS84"))
  #SLWMA
  #cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=16N ellps=WGS84"))
  
  # Put coordinates into date frame
  df2 <- as.data.frame(coordinates(cord.UTM))
  df2$datetime <- xy$datetime
  df2$ID <- xy$Collar_ID
  
  #-------------------------------------
  # Caluculate revisits meters radius
  #-------------------------------------
  
  # 100 m radius to find foraging sites
  visit <- getRecursions(df2, 250)
 
  # calculate mean time between visits for each location
  time <- visit$revisitStats %>% group_by(coordIdx) %>% 
    dplyr:: summarise(avg_between_days = mean(timeSinceLastVisit, na.rm=TRUE))

  # Pull values for each day
  map.df <- as(df2,'data.frame')
  map.df$revisits <- visit$revisits
  map.df$residenceTime <- visit$residenceTime
  map.df$avg_time_between <- time$avg_between_days 
  
  #create column with date
  map.df$date <- strptime(paste(map.df$datetime),"%Y-%m-%d")
  
  #add column of day of brooding
  map.df <- map.df %>% 
    mutate(day = dense_rank(date)) 
  
  #subset to 28 days brooding
  map.df <- subset(map.df, day <29)
  #Create name for output file
  
  fname <- as.character(paste(file_list[i]))
  
  #Create output file with above name
  
  write.csv(map.df, file = fname)
  #End for loop
  }

###############################
# Merge HMM to Recursion      #
###############################

mobile <- read.csv("E:/Brood_Recursion/Recursive_Analysis/Revisit_Stats_Full/All_Birds/All_Forage.csv")

hmm <- read.csv("E:/Brood_Recursion/HMM_Output/Final_HMM.csv")

mobile <- mobile[order(mobile$ID, mobile$datetime),]

hmm <- hmm[order(hmm$ID, hmm$time),]

# Deal with multiple broods
ID60346 <- subset(hmm, ID == 60346)
ID60343 <- subset(hmm, ID == 60343)
ID60273 <- subset(hmm, ID == 60273)
ID60327 <- subset(hmm, ID == 60327)
ID70578 <- subset(hmm, ID == 60578)

# Change IDs to second brood identifier
ID60346$ID <- 70346
ID60273$ID <- 70273
ID60343$ID <- 70343
ID60327$ID <- 70327
ID70578$ID <- 70578

# Put back into the HMM frame
hmm <- rbind(hmm,ID60327,ID70578,ID60273,ID60343,ID60346)


# Create date classes
mobile$date <- as.POSIXct(mobile$date, format = "%m/%d/%Y")
hmm$date <- as.POSIXct(hmm$time, format = "%Y-%m-%d")

# Get just dates considered brooding
mobile$concatenated <- paste(mobile$date, mobile$ID)
hmm$concatenated <- paste(hmm$date, hmm$ID)

# Match columns by ID and Date time to merge columns to main database
mobile$state <- hmm$state[sapply(mobile$concatenated, function(x) match(x, hmm$concatenated))]

# subset to just mobile
mobile_df <- subset(mobile, state==3)

###############################################################################
######################### Run for restricted state ############################
###############################################################################

# Begin loop
for (i in 1:length(file_list)){
  #read in .csv file
  df <- read.csv(paste(folder, file_list[i], sep=''))
  
  head(df)
  # make lat and lon numeric
  df$Lat <- as.numeric(df$Latitude)
  df$Lon <- as.numeric(df$Longitude)
 
  # Omit N/A from WEBB, SC data
  df <- na.omit(df)
  
  # Convert to POSIXct
  df$datetime <- as.POSIXct(strptime(paste(df$LMT_Date, df$LMT_Time),"%d.%m.%Y %H:%M"))
  xy <- df[,c(5,1)]
  
  # Creates spatial point for UTM conversion
  xy2 <- df[,c(12,11)]
  head(xy2)
  cord.dec <- SpatialPoints(cbind(xy2$Lon, xy2$Lat), proj4string = CRS("+proj=longlat"))
  
  # Transforming coordinate to UTM using EPSG=32615 for WGS=84, UTM Zone=15N)
  #Louisiana
  #cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=15N ellps=WGS84"))
  #SRS, CC, BFG, WEBB
  cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=17N ellps=WGS84"))
  #SLWMA
  #cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=16N ellps=WGS84"))
  
  # Put coordinates into date frame
  df2 <- as.data.frame(coordinates(cord.UTM))
  df2$datetime <- xy$datetime
  df2$ID <- xy$Collar_ID
  
  #-------------------------------------
  # Caluculate revisits meter radius
  #-------------------------------------
  
  # 100 m radius to find foraging sites
  
  visit <- getRecursions(df2, 90)
  
  # calculate mean time between visits for each location
  time <- visit$revisitStats %>% group_by(coordIdx) %>% 
    dplyr:: summarise(avg_between_days = mean(timeSinceLastVisit, na.rm=TRUE))
  
  # Pull values for each day
  
  map.df <- as(df2,'data.frame')
  map.df$revisits <- visit$revisits
  map.df$residenceTime <- visit$residenceTime
  map.df$avg_time_between <- time$avg_between_days 
  
  #create column with date
  map.df$date <- strptime(paste(map.df$datetime),"%Y-%m-%d")
  
  #add column of day of brooding
  map.df <- map.df %>% 
    mutate(day = dense_rank(date)) 
  
  #subset to 28 days brooding
  map.df <- subset(map.df, day <29)
  #Create name for output file
  
  fname <- as.character(paste(file_list[i]))
  
  #Create output file with above name
  
  write.csv(map.df, file = fname)
  #End for loop
  #End for loop
}

###############################
# Merge HMM to recursion      #
###############################

restricted <- read.csv("E:/Brood_Recursion/Recursive_Analysis/Revisit_Stats_Full/All_Birds/All_Forage.csv")

# Create date classes
restricted$date <- as.POSIXct(restricted$date, format = "%m/%d/%Y")
hmm$date <- as.POSIXct(hmm$time, format = "%Y-%m-%d")

# Get just dates considered brooding
restricted$concatenated <- paste(restricted$date, restricted$ID)
hmm$concatenated <- paste(hmm$date, hmm$ID)

# Match columns by ID and Date time to merge columns to main database
restricted$state <- hmm$state[sapply(restricted$concatenated, function(x) match(x, hmm$concatenated))]

# subset to just restricted movements
restricted_df <- subset(restricted, state==1 | state==2)

## Create behavioral state recursion master sheet
master <- rbind(forage_df, walk_df)

write.csv(master, "E:/Brood_Recursion/Recurse_State_Master_Sheet.csv")
