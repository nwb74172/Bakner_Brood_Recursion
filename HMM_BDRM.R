# Plotting package
library(ggplot2)
# Movement modelling packages
library(momentuHMM)
library(foieGras)
library(adehabitatLT)
library(lubridate)
library(dplyr)
library(splitstackshape)
# GIS packages
library(sf)
library(sp)
# Find means and sd for states
library(mixtools)
# Load functions for later
source("E:/Brood_Recursion/HMM_Files/function.R")

# Load data from Movebank URL
raw <- read.csv("E:/Brood_Recursion/Clean_CSV_Full_Season_Brooding_Females/All_Broods_Cleaned_Edited/All_Birds.csv")


raw$time <- as.POSIXct(strptime(paste(raw$LMT_Date, raw$LMT_Time),"%Y-%m-%d %H"))
raw <- subset(raw, raw$Latitude > 0)
raw <- subset(raw, DOP < 7)
unique(raw$Site)

srs_raw <- subset(raw, Site == "SRS")

webb_raw <- subset(raw, Site == "WEBB")

cc_raw <- subset(raw, Site == "CC")

bfg_raw <- subset(raw, Site == "BFG")

slwma_raw <- subset(raw, Site == "SLWMA")

knf_raw <- subset(raw, Site == "KNF")

# Keep relevant columns: ID, time, longitude, latitude, temperature
data <- srs_raw[, c(11, 12, 3, 4)]
colnames(data) <- c("ID", "time", "lat", "lon")
data$time <- as.POSIXct(data$time)

data1 <- webb_raw[, c(11, 12, 3, 4)]
colnames(data1) <- c("ID", "time", "lat", "lon")
data1$time <- as.POSIXct(data1$time)

data2 <- cc_raw[, c(11, 12, 3, 4)]
colnames(data2) <- c("ID", "time", "lat", "lon")
data2$time <- as.POSIXct(data2$time)

data3 <- bfg_raw[, c(11, 12, 3, 4)]
colnames(data3) <- c("ID", "time", "lat", "lon")
data3$time <- as.POSIXct(data3$time)

data4 <- slwma_raw[, c(11, 12, 3, 4)]
colnames(data4) <- c("ID", "time", "lat", "lon")
data4$time <- as.POSIXct(data4$time)

data5 <- knf_raw[, c(11, 12, 3, 4)]
colnames(data5) <- c("ID", "time", "lat", "lon")
data5$time <- as.POSIXct(data5$time)

# Just keep 2000 observations to save time with model fitting
#data <- data_all[c(which(data_all$ID == unique(data_all$ID)[5])[1:1000],
#                   which(data_all$ID == unique(data_all$ID)[6])[1:1000]),]
# Subset to March 15
#create column with date

data$date <- format(as.Date(data$time), "%m-%d")
data$date <- as.Date( as.character(data$date), "%m-%d")
data <- subset(data, date > as.Date("2022-03-15"))
data$date <- NULL

data1$date <- format(as.Date(data1$time), "%m-%d")
data1$date <- as.Date( as.character(data1$date), "%m-%d")
data1 <- subset(data1, date > as.Date("2022-03-15"))
data1$date <- NULL

data2$date <- format(as.Date(data2$time), "%m-%d")
data2$date <- as.Date( as.character(data2$date), "%m-%d")
data2 <- subset(data2, date > as.Date("2022-03-15"))
data2$date <- NULL

data3$date <- format(as.Date(data3$time), "%m-%d")
data3$date <- as.Date( as.character(data3$date), "%m-%d")
data3 <- subset(data3, date > as.Date("2022-03-15"))
data3$date <- NULL

data4$date <- format(as.Date(data4$time), "%m-%d")
data4$date <- as.Date( as.character(data4$date), "%m-%d")
data4 <- subset(data4, date > as.Date("2022-03-15"))
data4$date <- NULL

data5$date <- format(as.Date(data5$time), "%m-%d")
data5$date <- as.Date( as.character(data5$date), "%m-%d")
data5 <- subset(data5, date > as.Date("2022-03-15"))
data5$date <- NULL

# Project to UTM
llcoord <- st_as_sf(data[, c("lon", "lat")], coords = c("lon", "lat"), 
                    crs = CRS("+proj=longlat +datum=WGS84"))
utmcoord <- st_transform(llcoord, crs = CRS("+proj=utm +zone=17 +datum=WGS84"))

llcoord1 <- st_as_sf(data1[, c("lon", "lat")], coords = c("lon", "lat"), 
                    crs = CRS("+proj=longlat +datum=WGS84"))
utmcoord1 <- st_transform(llcoord1, crs = CRS("+proj=utm +zone=17 +datum=WGS84"))

llcoord2 <- st_as_sf(data2[, c("lon", "lat")], coords = c("lon", "lat"), 
                     crs = CRS("+proj=longlat +datum=WGS84"))
utmcoord2 <- st_transform(llcoord2, crs = CRS("+proj=utm +zone=17 +datum=WGS84"))

llcoord3 <- st_as_sf(data3[, c("lon", "lat")], coords = c("lon", "lat"), 
                     crs = CRS("+proj=longlat +datum=WGS84"))
utmcoord3 <- st_transform(llcoord3, crs = CRS("+proj=utm +zone=17 +datum=WGS84"))

llcoord4 <- st_as_sf(data4[, c("lon", "lat")], coords = c("lon", "lat"), 
                     crs = CRS("+proj=longlat +datum=WGS84"))
utmcoord4 <- st_transform(llcoord4, crs = CRS("+proj=utm +zone=16 +datum=WGS84"))

llcoord5 <- st_as_sf(data5[, c("lon", "lat")], coords = c("lon", "lat"), 
                     crs = CRS("+proj=longlat +datum=WGS84"))
utmcoord5 <- st_transform(llcoord5, crs = CRS("+proj=utm +zone=15 +datum=WGS84"))

# Add Easting-Northing to data (in km)
data[, c("x", "y")] <- st_coordinates(utmcoord)/1000

data1[, c("x", "y")] <- st_coordinates(utmcoord1)/1000

data2[, c("x", "y")] <- st_coordinates(utmcoord2)/1000

data3[, c("x", "y")] <- st_coordinates(utmcoord3)/1000

data4[, c("x", "y")] <- st_coordinates(utmcoord4)/1000

data5[, c("x", "y")] <- st_coordinates(utmcoord5)/1000

data_split <- split_at_gap(data = data, max_gap = 2*60, shortest_track = 2*60)

data_split1 <- split_at_gap(data = data1, max_gap = 2*60, shortest_track = 2*60)

data_split2 <- split_at_gap(data = data2, max_gap = 2*60, shortest_track = 2*60)

data_split3 <- split_at_gap(data = data3, max_gap = 2*60, shortest_track = 2*60)

data_split4 <- split_at_gap(data = data4, max_gap = 2*60, shortest_track = 2*60)

data_split5 <- split_at_gap(data = data5, max_gap = 2*60, shortest_track = 2*60)

#datatest <- data[!duplicated(data[c(1,2)]),]
# remove duplicate times
data_split <- data_split %>%
  distinct(ID, time, .keep_all = TRUE)

data_split1 <- data_split1 %>%
  distinct(ID, time, .keep_all = TRUE)

data_split2 <- data_split2 %>%
  distinct(ID, time, .keep_all = TRUE)

data_split3 <- data_split3 %>%
  distinct(ID, time, .keep_all = TRUE)

data_split4 <- data_split4 %>%
  distinct(ID, time, .keep_all = TRUE)

data_split5 <- data_split5 %>%
  distinct(ID, time, .keep_all = TRUE)

# Create adehabitat trajectory padded with NAs
data_ade <- setNA(ltraj = as.ltraj(xy = data_split[, c("x", "y")], 
                                   date = data_split$time, 
                                   id = data_split$ID),
                  date.ref = data_split$time[1], 
                  dt = 60, tol = 5, units = "min")

data_ade1 <- setNA(ltraj = as.ltraj(xy = data_split1[, c("x", "y")], 
                                   date = data_split1$time, 
                                   id = data_split1$ID),
                  date.ref = data_split1$time[1], 
                  dt = 60, tol = 5, units = "min")

data_ade2 <- setNA(ltraj = as.ltraj(xy = data_split2[, c("x", "y")], 
                                   date = data_split2$time, 
                                   id = data_split2$ID),
                  date.ref = data_split2$time[1], 
                  dt = 60, tol = 5, units = "min")

data_ade3 <- setNA(ltraj = as.ltraj(xy = data_split3[, c("x", "y")], 
                                   date = data_split3$time, 
                                   id = data_split3$ID),
                  date.ref = data_split3$time[1], 
                  dt = 60, tol = 5, units = "min")

data_ade4 <- setNA(ltraj = as.ltraj(xy = data_split4[, c("x", "y")], 
                                   date = data_split4$time, 
                                   id = data_split4$ID),
                  date.ref = data_split4$time[1], 
                  dt = 60, tol = 5, units = "min")

data_ade5 <- setNA(ltraj = as.ltraj(xy = data_split5[, c("x", "y")], 
                                   date = data_split5$time, 
                                   id = data_split5$ID),
                  date.ref = data_split5$time[1], 
                  dt = 60, tol = 5, units = "min")

# Transform back to dataframe
data_na <- ld(data_ade)[, c("id", "x", "y", "date")]
colnames(data_na) <- c("ID", "x", "y", "time")

data_na1 <- ld(data_ade1)[, c("id", "x", "y", "date")]
colnames(data_na1) <- c("ID", "x", "y", "time")

data_na2 <- ld(data_ade2)[, c("id", "x", "y", "date")]
colnames(data_na2) <- c("ID", "x", "y", "time")

data_na3 <- ld(data_ade3)[, c("id", "x", "y", "date")]
colnames(data_na3) <- c("ID", "x", "y", "time")

data_na4 <- ld(data_ade4)[, c("id", "x", "y", "date")]
colnames(data_na4) <- c("ID", "x", "y", "time")

data_na5 <- ld(data_ade5)[, c("id", "x", "y", "date")]
colnames(data_na5) <- c("ID", "x", "y", "time")

# Combine all the data

final_data <- rbind(data_na, data_na1, data_na2, data_na3, data_na4, data_na5)
head(final_data)

# Prepare data for HMM (compute step lengths and turning angles)

movement_df <- prepData(final_data, type = "UTM")
summary(movement_df$step)

# remove points greater than 8.5

movement_df$angle = ifelse(movement_df$step > 4.5, NA, movement_df$angle) 
movement_df$step = ifelse(movement_df$step > 4.5, NA, movement_df$step)

# Explore data
hist(movement_df$step)
hist(movement_df$angle, breaks = seq(-pi, pi, length = 15))

# Observation distributions (step lengths and turning angles)
dist <- list(step = "gamma", angle = "vm")

# Fit a 3-state HMM

hmm <- fitHMM(data = movement_df,
              nbStates = 3,
              dist = list(step = "gamma", angle = "vm"),
              estAngleMean = list(angle = TRUE),
              Par0 = list(step = c(mean_1 = 0.027, 
                                   mean_2 = 0.150,
                                   mean_3 = 0.400,
                                   sd_1 = 0.027, 
                                   sd_2 = 0.150,
                                   sd_3 = 1.00,
                                   zeromass_1 = 0, 
                                   zeromass_2 = 0,
                                   zeromass_3 = 0),
                          angle = c(mean_1 = 3.14, 
                                    mean_2 = 2.5,
                                    mean_2 = 0.001,
                                    concentration_1 = 0.1,
                                    concentration_2 = 0.5,
                                    concentration_3 = 0.99)))

# Print parameter estimates
hmm

# Save most likely state sequences from 2-state and 3-state models
movement_df$state <- factor(viterbi(hmm))

# Split ID column to just get ID
splittest <- cSplit(srs,c("ID"), c("-"))
splittest$ID <- splittest$ID_1
splittest$ID_2 <- NULL
splittest$ID_1 <- NULL

# write csv
write.csv(splittest, "D:/Brood_Recursion/HMM_Output/Final_HMM.csv")


