############################
#####LOAD PACKAGES##########
############################
library(dplyr)
library(tidyr)
library(plyr)
library(amt)
library(adehabitatHR)
library(adehabitatLT)
library(circular)
library(raster)
library(rgeos)
library(geosphere)
library(rgdal)
library(sf)
memory.limit(100000)

# All birds
dat <- read.csv("E:/Brood_Recursion/Recurse_State_Master_Sheet.csv")
head(dat)
#miss <- subset(dat, ID == 61138)
# One bird
indiv <- subset(dat, site=="WEBB")

# Remove birds with no enough data
#indiv <-indiv[!(indiv$ID== 60328),]

# Set the date and time
indiv$datetime <- as.POSIXct(indiv$datetime, format = "%m/%d/%Y %H:%M")
indiv$yr <- format(indiv$datetime, format = "%Y")
head(indiv)

# Subset by year
#indiv <- subset(indiv, yr=="2015")

# Make a spatial data frame
df.sp <- indiv
coordinates(df.sp) <- c("coords.x1","coords.x2")

# Define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB SRS CC 
proj4string(df.sp) <- CRS("+proj=utm +zone=17N ellps=WGS84")

# Transfer the projection
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))

# Add to the data frame
turkey.sp <- as.data.frame(df.sp)
head(turkey.sp)
indiv$easting <- turkey.sp$coords.x1
indiv$northing <- turkey.sp$coords.x2
head(indiv)

# Get data set up for running the loop
#ids <- unique(indiv$ID)
#memory = random_memory=vector("list",length(ids))
#names(memory) <- ids

#Convert to a tibble file
df.t <- as_tibble(indiv)

# Now we can create a track from our cleaned tibble
# Louisiana
#trk <- df.t %>% make_track(easting, northing, datetime, burst_ = ID, crs = CRS("+proj=utm +zone=15N ellps=WGS84"), all_cols = TRUE)
# WEBB CC SRS
trk <- df.t %>% make_track(easting, northing, datetime, burst_ = ID, crs = CRS("+proj=utm +zone=17N ellps=WGS84"), all_cols = TRUE)

head(trk)
trk$burst_ <- trk$ID
summarize_sampling_rate(trk)

# Nest track by ID
trk <- trk %>% nest(data = -"ID")
ssf1 <- trk %>% mutate(steps = map(data, function(x) 
  x %>% steps_by_burst(keep_cols = 'end')))

# Unnest to get complete data frame and look at plot of step length
ssf1 <- ssf1 %>% dplyr::select(ID, steps) %>% tidyr::unnest(cols = steps)

# Get random steps
ssf1 <-ssf1 %>% random_steps(n = 100)
head(ssf1)

# make into a data frame
ran <- as.data.frame(ssf1)

# Get data set up for running the loop
ids <- unique(ran$ID)
memory = random_memory=vector("list",length(ids))
names(memory) <- ids

# Make a spatial data frame
df.sp <- ran
head(df.sp)
coordinates(df.sp) <- c("x2_","y2_")

# define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB, SRS, CC
proj4string(df.sp) <- CRS("+proj=utm +zone=17N ellps=WGS84")
#df.sp$dt_ <- NULL

# transfer the projection
# Louisiana
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))
# WEBB, SC, SRS
df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=17N ellps=WGS84"))


for(j in 1:length(ids)){ 
  temp <- indiv[indiv$ID==ids[j],]
  temp2 <- df.sp[df.sp$ID==ids[j],]
  unique(temp$ID)
  unique(temp2$ID)
  # Make a trajectory
  traj <- as.ltraj(temp[,c("easting","northing")],id=temp$ID,date=temp$datetime)
  
  # BRB.D() found average D = 1, Tmax is based on mean residence time, radius 
  # is based on behavioral state, hmin is the location error, max time until left
  # consider it an hour based on GPS unit.
  
  # Memory estimation type ID = residence time, type RD = revisits
  # Behavioral State Walking
  id_rmw <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="ID", radius =250,maxt =3600,grid=200)
  id_rmw1 <- getvolumeUD(id_rmw, standardize = TRUE)
  # Louisiana
  #proj4string(id_rmw1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB SC CC
  proj4string(id_rmw1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  plot(id_rmw1)
  rast1 <- raster(as(id_rmw1,"SpatialPixelsDataFrame"))
  
  rd_rmw <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="RD", radius =250,maxt =3600,grid=200)
  rd_rmw1 <- getvolumeUD(rd_rmw, standardize = TRUE)
  plot(rd_rmw1)
  # Louisiana
  # proj4string(rd_rmw1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(rd_rmw1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  rast2 <- raster(as(rd_rmw1,"SpatialPixelsDataFrame"))
  
  # Write raster files for each individual
  writeRaster(rast1, filename=paste0("E:/Brood_Covariates/Cognitive_Map/ID/WEBB", "_", ids[j], "_residence_time_walking.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(rast2, filename=paste0("E:/Brood_Covariates/Cognitive_Map/RD/WEBB", "_", ids[j], "_revist_walking.tif"), format="GTiff", overwrite=TRUE)
  
  # Behavioral State Foraging
  id_rmf <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="ID", radius =90,maxt =3600,grid=200)
  id_rmf1 <- getvolumeUD(id_rmf, standardize = TRUE)
  # Louisiana
  #proj4string(id_rmf1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(id_rmf1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  plot(id_rmf1)
  rast3 <- raster(as(id_rmf1,"SpatialPixelsDataFrame"))
  
  rd_rmf <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="RD", radius =90,maxt =3600,grid=200)
  rd_rmf1 <- getvolumeUD(rd_rmf, standardize = TRUE)
  plot(rd_rmf1)
  # Louisiana
  #proj4string(rd_rmf1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(rd_rmf1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  rast4 <- raster(as(rd_rmf1,"SpatialPixelsDataFrame"))
  
  # Write raster files for each individual
  writeRaster(rast3, filename=paste0("E:/Brood_Covariates/Cognitive_Map/ID/WEBB", "_", ids[j], "_residence_time_foraging.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(rast4, filename=paste0("E:/Brood_Covariates/Cognitive_Map/RD/WEBB", "_", ids[j], "_revist_foraging.tif"), format="GTiff", overwrite=TRUE)
  
  # extract values from ID and RD
  # Behavioral state walking
  temp2$residence_raster_walking <- raster::extract(rast1, temp2)
  temp2$revisit_raster_walking <- raster::extract(rast2, temp2)
  
  temp2$residence_raster_foraging <- raster::extract(rast3, temp2)
  temp2$revisit_raster_foraging <- raster::extract(rast4, temp2)
  
  temp2 <- as.data.frame(temp2)
  
  assign(paste0(temp2$site,temp2$ID), temp2)
}

df <- rbind( WEBB60187, WEBB60191, WEBB60196,
             WEBB60205, WEBB60216, WEBB60217, WEBB60221, WEBB60644,WEBB60647, WEBB60652, WEBB60658,
             WEBB60667, WEBB60672, WEBB60676, WEBB60682, WEBB60688, WEBB60690, WEBB60693,
             WEBB60696, WEBB60918, WEBB60919, WEBB61138)
head(df)
summary(df)

df <- as.data.frame(df)
head(df)
df$residence_raster_walking[is.na(df$residence_raster_walking)] <- 100
df$residence_raster_foraging[is.na(df$residence_raster_foraging)] <- 100

df$revisit_raster_walking[is.na(df$revisit_raster_walking)] <- 100
df$revisit_raster_foraging[is.na(df$revisit_raster_foraging)] <- 100

summary(df)
#df.test <- subset(df, case_=="TRUE")
#df.test <- subset(df.test, state==4)
#summary(df.test)
#library(ggplot2)





######## Step Selection Function #########
# Make a spatial data frame
df.sp <- df
head(df.sp)
coordinates(df.sp) <- c("x2_","y2_")

# Define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB SRS CC 
proj4string(df.sp) <- CRS("+proj=utm +zone=17N ellps=WGS84")

# Transfer the projection
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))

# subset out years
ssf15 <- subset(df.sp, yr==2015)
ssf16 <- subset(df.sp, yr==2016)
ssf18 <- subset(df.sp, yr==2018)

# extract covariates
proads <- raster("E:/Brood_Covariates/Roads/Webb_roads/Webb_prime_rd_dist.tif")
sroads <- raster("E:/Brood_Covariates/Roads/Webb_roads/Webb_2nd_rd_dist.tif")
dist_hw <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_HW/dist_hw_15.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_infra <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_infra/dist_infra_15.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_mix <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_Mixed/dist_mix_15.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_open <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_Open/dist_open_15.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_pine <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_Pine/dist_pine_15.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_shrub <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_Shrub/dist_shrub_15.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_water <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_Water/dist_water_15.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
ndvi <- raster("E:/Brood_Covariates/NDVI/Webb/Webb_2015/ndvi_WEBB_2015.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)


ssf15$proads <- raster::extract(proads, ssf15) 
ssf15$sroads <- raster::extract(sroads, ssf15)
ssf15$dist_hw <- raster::extract(dist_hw, ssf15)
ssf15$dist_infra <- raster::extract(dist_infra, ssf15)
ssf15$dist_mixed <- raster::extract(dist_mix, ssf15)
ssf15$dist_open <- raster::extract(dist_open, ssf15)
ssf15$dist_pine <- raster::extract(dist_pine, ssf15)
ssf15$dist_shrub <- raster::extract(dist_shrub, ssf15)
ssf15$dist_water <- raster::extract(dist_water, ssf15)
ssf15$ndvi <- raster::extract(ndvi, ssf15)
ssf15$cos_ta_ <- cos(ssf15$ta_)
ssf15$log_sl_ <- log(ssf15$sl_)


walking15 <- as.data.frame(subset(ssf15, state==3))
restricted15 <- as.data.frame(subset(ssf15, state==1 | state==2))


# extract covariates
proads <- raster("E:/Brood_Covariates/Roads/Webb_roads/Webb_prime_rd_dist.tif")
sroads <- raster("E:/Brood_Covariates/Roads/Webb_roads/Webb_2nd_rd_dist.tif")
dist_hw <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_HW/dist_hw_16.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_infra <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_infra/dist_infra_16.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_mix <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_Mixed/dist_mix_16.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_open <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_Open/dist_open_16.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_pine <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_Pine/dist_pine_16.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_shrub <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_Shrub/dist_shrub_16.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_water <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_Water/dist_water_16.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
ndvi <- raster("E:/Brood_Covariates/NDVI/Webb/Webb_2016/ndvi_WEBB_2016.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)


ssf16$proads <- raster::extract(proads, ssf16) 
ssf16$sroads <- raster::extract(sroads, ssf16)
ssf16$dist_hw <- raster::extract(dist_hw, ssf16)
ssf16$dist_infra <- raster::extract(dist_infra, ssf16)
ssf16$dist_mixed <- raster::extract(dist_mix, ssf16)
ssf16$dist_open <- raster::extract(dist_open, ssf16)
ssf16$dist_pine <- raster::extract(dist_pine, ssf16)
ssf16$dist_shrub <- raster::extract(dist_shrub, ssf16)
ssf16$dist_water <- raster::extract(dist_water, ssf16)
ssf16$ndvi <- raster::extract(ndvi, ssf16)
ssf16$cos_ta_ <- cos(ssf16$ta_)
ssf16$log_sl_ <- log(ssf16$sl_)

#walking16 <- as.data.frame(subset(ssf16, state==1))
#restricted16 <- as.data.frame(subset(ssf16, state==2 | state==3))
#roost16 <- as.data.frame(subset(ssf16, state==4))

# extract covariates
proads <- raster("E:/Brood_Covariates/Roads/Webb_roads/Webb_prime_rd_dist.tif")
sroads <- raster("E:/Brood_Covariates/Roads/Webb_roads/Webb_2nd_rd_dist.tif")
dist_hw <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_HW/dist_hw_18.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_infra <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_infra/dist_infra_17.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_mix <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_Mixed/dist_mix_18.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_open <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_Open/dist_open_18.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_pine <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_Pine/dist_pine_17.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_shrub <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_Shrub/dist_shrub_18.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_water <- raster("E:/Brood_Covariates/Landcover/Landcover_Webb/Final_Dist_Water/dist_water_18.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
ndvi <- raster("E:/Brood_Covariates/NDVI/Webb/Webb_2018/ndvi_WEBB_2018.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)


ssf18$proads <- raster::extract(proads, ssf18) 
ssf18$sroads <- raster::extract(sroads, ssf18)
ssf18$dist_hw <- raster::extract(dist_hw, ssf18)
ssf18$dist_infra <- raster::extract(dist_infra, ssf18)
ssf18$dist_mixed <- raster::extract(dist_mix, ssf18)
ssf18$dist_open <- raster::extract(dist_open, ssf18)
ssf18$dist_pine <- raster::extract(dist_pine, ssf18)
ssf18$dist_shrub <- raster::extract(dist_shrub, ssf18)
ssf18$dist_water <- raster::extract(dist_water, ssf18)
ssf18$ndvi <- raster::extract(ndvi, ssf18)
ssf18$cos_ta_ <- cos(ssf18$ta_)
ssf18$log_sl_ <- log(ssf18$sl_)

#walking16 <- as.data.frame(subset(ssf16, state==1))
#restricted16 <- as.data.frame(subset(ssf16, state==2 | state==3))
#roost16 <- as.data.frame(subset(ssf16, state==4))

############## Cedar Creek ############################
# All birds
dat <- read.csv("E:/Brood_Recursion/Recurse_State_Master_Sheet.csv")
head(dat)

# One bird
indiv <- subset(dat, site=="CC")
unique(indiv$ID)
# Remove birds with no enough data
#indiv <-indiv[!(indiv$ID== 60328),]

# Set the date and time
indiv$datetime <- as.POSIXct(indiv$datetime, format = "%m/%d/%Y %H:%M")
indiv$yr <- format(indiv$datetime, format = "%Y")
head(indiv)

# Subset by year
#indiv <- subset(indiv, yr=="2015")

# Make a spatial data frame
df.sp <- indiv
coordinates(df.sp) <- c("coords.x1","coords.x2")

# Define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB SRS CC 
proj4string(df.sp) <- CRS("+proj=utm +zone=17N ellps=WGS84")

# Transfer the projection
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))

# Add to the data frame
turkey.sp <- as.data.frame(df.sp)
head(turkey.sp)
indiv$easting <- turkey.sp$coords.x1
indiv$northing <- turkey.sp$coords.x2
head(indiv)

# Get data set up for running the loop
ids <- unique(indiv$ID)
memory = random_memory=vector("list",length(ids))
names(memory) <- ids

#Convert to a tibble file
df.t <- as_tibble(indiv)

# Now we can create a track from our cleaned tibble
# Louisiana
#trk <- df.t %>% make_track(easting, northing, datetime, burst_ = ID, crs = CRS("+proj=utm +zone=15N ellps=WGS84"), all_cols = TRUE)
# WEBB CC SRS
trk <- df.t %>% make_track(easting, northing, datetime, burst_ = ID, crs = CRS("+proj=utm +zone=17N ellps=WGS84"), all_cols = TRUE)

head(trk)
trk$burst_ <- trk$ID
summarize_sampling_rate(trk)

# Nest track by ID
trk <- trk %>% nest(data = -"ID")
ssf1 <- trk %>% mutate(steps = map(data, function(x) 
  x %>% steps_by_burst(keep_cols = 'end')))

# Unnest to get complete data frame and look at plot of step length
ssf1 <- ssf1 %>% dplyr::select(ID, steps) %>% tidyr::unnest(cols = steps)

# Get random steps
ssf1 <-ssf1 %>% random_steps(n = 100)
head(ssf1)
library(amt)
?random_steps
# make into a data frame
ran <- as.data.frame(ssf1)

# Make a spatial data frame
df.sp <- ran
head(ran)
coordinates(df.sp) <- c("x2_","y2_")

# define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB, SRS, CC
proj4string(df.sp) <- CRS("+proj=utm +zone=17N ellps=WGS84")
#df.sp$dt_ <- NULL

# transfer the projection
# Louisiana
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))
# WEBB, SC, SRS
df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=17N ellps=WGS84"))

# Set up code to create ID and RD for each individual
for(j in 1:length(ids)){ 
  temp <- indiv[indiv$ID==ids[j],]
  temp2 <- df.sp[df.sp$ID==ids[j],]
  unique(temp$ID)
  unique(temp2$ID)
  # Make a trajectory
  traj <- as.ltraj(temp[,c("easting","northing")],id=temp$ID,date=temp$datetime)
  
  # BRB.D() found average D = 1, Tmax is based on mean residence time, radius 
  # is based on behavioral state, hmin is the location error, max time until left
  # consider it an hour based on GPS unit.
  
  # Memory estimation type ID = residence time, type RD = revisits
  # Behavioral State Walking
  id_rmw <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="ID", radius =250,maxt =3600,grid=200)
  id_rmw1 <- getvolumeUD(id_rmw, standardize = TRUE)
  # Louisiana
  #proj4string(id_rmw1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB SC CC
  proj4string(id_rmw1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  plot(id_rmw1)
  rast1 <- raster(as(id_rmw1,"SpatialPixelsDataFrame"))
  
  rd_rmw <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="RD", radius =250,maxt =3600,grid=200)
  rd_rmw1 <- getvolumeUD(rd_rmw, standardize = TRUE)
  plot(rd_rmw1)
  # Louisiana
  # proj4string(rd_rmw1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(rd_rmw1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  rast2 <- raster(as(rd_rmw1,"SpatialPixelsDataFrame"))
  
  # Write raster files for each individual
  writeRaster(rast1, filename=paste0("E:/Brood_Covariates/Cognitive_Map/ID/CC", "_", ids[j], "_residence_time_walking.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(rast2, filename=paste0("E:/Brood_Covariates/Cognitive_Map/RD/CC", "_", ids[j], "_revist_walking.tif"), format="GTiff", overwrite=TRUE)
  
  # Behavioral State Foraging
  id_rmf <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="ID", radius =90,maxt =3600,grid=200)
  id_rmf1 <- getvolumeUD(id_rmf, standardize = TRUE)
  # Louisiana
  #proj4string(id_rmf1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(id_rmf1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  plot(id_rmf1)
  rast3 <- raster(as(id_rmf1,"SpatialPixelsDataFrame"))
  
  rd_rmf <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="RD", radius =90,maxt =3600,grid=200)
  rd_rmf1 <- getvolumeUD(rd_rmf, standardize = TRUE)
  plot(rd_rmf1)
  # Louisiana
  #proj4string(rd_rmf1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(rd_rmf1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  rast4 <- raster(as(rd_rmf1,"SpatialPixelsDataFrame"))
  
  # Write raster files for each individual
  writeRaster(rast3, filename=paste0("E:/Brood_Covariates/Cognitive_Map/ID/CC", "_", ids[j], "_residence_time_foraging.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(rast4, filename=paste0("E:/Brood_Covariates/Cognitive_Map/RD/CC", "_", ids[j], "_revist_foraging.tif"), format="GTiff", overwrite=TRUE)
  
  # extract values from ID and RD
  # Behavioral state walking
  temp2$residence_raster_walking <- raster::extract(rast1, temp2)
  temp2$revisit_raster_walking <- raster::extract(rast2, temp2)
  
  temp2$residence_raster_foraging <- raster::extract(rast3, temp2)
  temp2$revisit_raster_foraging <- raster::extract(rast4, temp2)
  
  temp2 <- as.data.frame(temp2)
  
  assign(paste0(temp2$site,temp2$ID), temp2)
}

df <- rbind( CC46213, CC46218, CC46222, CC46233, CC46237, CC46260, CC46563, CC46579,
             CC46580, CC46582, CC46583, CC46598, CC46611, CC47444, CC47447, CC47462,
             CC47469, CC47474, CC47475, CC61061, CC61068, CC61069, CC61077, CC61078,
             CC61079, CC61084, CC61089, CC61079, CC61078, CC61079, CC61084, CC61089,
             CC61097, CC61099, CC61102, CC61108, CC61111, CC61113, CC61254, CC61269,
             CC61301)
head(df)
summary(df)

df <- as.data.frame(df)

df$residence_raster_walking[is.na(df$residence_raster_walking)] <- 100
df$residence_raster_foraging[is.na(df$residence_raster_foraging)] <- 100

df$revisit_raster_walking[is.na(df$revisit_raster_walking)] <- 100
df$revisit_raster_foraging[is.na(df$revisit_raster_foraging)] <- 100

summary(df)

######## Step Selection Function #########
# Make a spatial data frame
df.sp <- df
head(df.sp)
coordinates(df.sp) <- c("x2_","y2_")

# Define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB SRS CC 
proj4string(df.sp) <- CRS("+proj=utm +zone=17N ellps=WGS84")

# Transfer the projection
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))

# subset out years
ccssf17 <- subset(df.sp, yr==2017)
ccssf18 <- subset(df.sp, yr==2018)
ccssf19 <- subset(df.sp, yr==2019)
ccssf20 <- subset(df.sp, yr==2020)
ccssf21 <- subset(df.sp, yr==2021)

# extract covariates
proads <- raster("E:/Brood_Covariates/Roads/CC_roads/PrimeRdRast(1).tif")
sroads <- raster("E:/Brood_Covariates/Roads/CC_roads/SecondRdRast(1).tif")
dist_hw <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_HW/dist_hw_17.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_infra <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_infra/dist_infra_17.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_mix <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Mixed/dist_mixed_17.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_open <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Open/dis_open_17.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_pine <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Pine/dist_pine_17.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_shrub <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Shrub/dist_shrub_17.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_water <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Water/dist_water_17.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
ndvi <- raster("E:/Brood_Covariates/NDVI/CC/CC_2017/ndvi_cc_2017.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)


ccssf17$proads <- raster::extract(proads, ccssf17) 
ccssf17$sroads <- raster::extract(sroads, ccssf17)
ccssf17$dist_hw <- raster::extract(dist_hw, ccssf17)
ccssf17$dist_infra <- raster::extract(dist_infra, ccssf17)
ccssf17$dist_mixed <- raster::extract(dist_mix, ccssf17)
ccssf17$dist_open <- raster::extract(dist_open, ccssf17)
ccssf17$dist_pine <- raster::extract(dist_pine, ccssf17)
ccssf17$dist_shrub <- raster::extract(dist_shrub, ccssf17)
ccssf17$dist_water <- raster::extract(dist_water, ccssf17)
ccssf17$ndvi <- raster::extract(ndvi, ccssf17)
ccssf17$cos_ta_ <- cos(ccssf17$ta_)
ccssf17$log_sl_ <- log(ccssf17$sl_)


#walking15 <- as.data.frame(subset(ssf15, state==3))
#restricted15 <- as.data.frame(subset(ssf15, state==1 | state==2))
#roost15 <- as.data.frame(subset(ssf15, state==4))

# extract covariates
proads <- raster("E:/Brood_Covariates/Roads/CC_roads/PrimeRdRast(1).tif")
sroads <- raster("E:/Brood_Covariates/Roads/CC_roads/SecondRdRast(1).tif")
dist_hw <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_HW/dist_hw_18.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_infra <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_infra/dist_infra_18.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_mix <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Mixed/dist_mixed_18.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_open <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Open/dist_open_18.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_pine <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Pine/dist_pine_18.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_shrub <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Shrub/dist_shrub_18.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_water <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Water/dist_water_18.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
ndvi <- raster("E:/Brood_Covariates/NDVI/CC/CC_2018/ndvi_cc_2018.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)


ccssf18$proads <- raster::extract(proads, ccssf18) 
ccssf18$sroads <- raster::extract(sroads, ccssf18)
ccssf18$dist_hw <- raster::extract(dist_hw, ccssf18)
ccssf18$dist_infra <- raster::extract(dist_infra, ccssf18)
ccssf18$dist_mixed <- raster::extract(dist_mix, ccssf18)
ccssf18$dist_open <- raster::extract(dist_open, ccssf18)
ccssf18$dist_pine <- raster::extract(dist_pine, ccssf18)
ccssf18$dist_shrub <- raster::extract(dist_shrub, ccssf18)
ccssf18$dist_water <- raster::extract(dist_water, ccssf18)
ccssf18$ndvi <- raster::extract(ndvi, ccssf18)
ccssf18$cos_ta_ <- cos(ccssf18$ta_)
ccssf18$log_sl_ <- log(ccssf18$sl_)

#walking16 <- as.data.frame(subset(ssf16, state==1))
#restricted16 <- as.data.frame(subset(ssf16, state==2 | state==3))
#roost16 <- as.data.frame(subset(ssf16, state==4))

# extract covariates
proads <- raster("E:/Brood_Covariates/Roads/CC_roads/PrimeRdRast(1).tif")
sroads <- raster("E:/Brood_Covariates/Roads/CC_roads/SecondRdRast(1).tif")
dist_hw <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_HW/dist_hw_19.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_infra <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_infra/dist_infra_19.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_mix <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Mixed/dist_mixed_19.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_open <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Open/dist_open_19.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_pine <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Pine/dist_pine_19.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_shrub <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Shrub/dist_shrub_19.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_water <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Water/dist_water_19.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
ndvi <- raster("E:/Brood_Covariates/NDVI/CC/CC_2019/ndvi_cc_2019.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)


ccssf19$proads <- raster::extract(proads, ccssf19) 
ccssf19$sroads <- raster::extract(sroads, ccssf19)
ccssf19$dist_hw <- raster::extract(dist_hw, ccssf19)
ccssf19$dist_infra <- raster::extract(dist_infra, ccssf19)
ccssf19$dist_mixed <- raster::extract(dist_mix, ccssf19)
ccssf19$dist_open <- raster::extract(dist_open, ccssf19)
ccssf19$dist_pine <- raster::extract(dist_pine, ccssf19)
ccssf19$dist_shrub <- raster::extract(dist_shrub, ccssf19)
ccssf19$dist_water <- raster::extract(dist_water, ccssf19)
ccssf19$ndvi <- raster::extract(ndvi, ccssf19)
ccssf19$cos_ta_ <- cos(ccssf19$ta_)
ccssf19$log_sl_ <- log(ccssf19$sl_)

#walking16 <- as.data.frame(subset(ssf16, state==1))
#restricted16 <- as.data.frame(subset(ssf16, state==2 | state==3))
#roost16 <- as.data.frame(subset(ssf16, state==4))

# extract covariates
proads <- raster("E:/Brood_Covariates/Roads/CC_roads/PrimeRdRast(1).tif")
sroads <- raster("E:/Brood_Covariates/Roads/CC_roads/SecondRdRast(1).tif")
dist_hw <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_HW/dist_hw_20.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_infra <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_infra/dist_infra_20.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_mix <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Mixed/dist_mixed_20.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_open <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Open/dist_open_20.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_pine <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Pine/dist_pine_20.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_shrub <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Shrub/dist_shrub_20.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_water <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Water/dist_water_20.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
ndvi <- raster("E:/Brood_Covariates/NDVI/CC/CC_2020/ndvi_cc_2020.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)


ccssf20$proads <- raster::extract(proads, ccssf20) 
ccssf20$sroads <- raster::extract(sroads, ccssf20)
ccssf20$dist_hw <- raster::extract(dist_hw, ccssf20)
ccssf20$dist_infra <- raster::extract(dist_infra, ccssf20)
ccssf20$dist_mixed <- raster::extract(dist_mix, ccssf20)
ccssf20$dist_open <- raster::extract(dist_open, ccssf20)
ccssf20$dist_pine <- raster::extract(dist_pine, ccssf20)
ccssf20$dist_shrub <- raster::extract(dist_shrub, ccssf20)
ccssf20$dist_water <- raster::extract(dist_water, ccssf20)
ccssf20$ndvi <- raster::extract(ndvi, ccssf20)
ccssf20$cos_ta_ <- cos(ccssf20$ta_)
ccssf20$log_sl_ <- log(ccssf20$sl_)

#walking16 <- as.data.frame(subset(ssf16, state==1))
#restricted16 <- as.data.frame(subset(ssf16, state==2 | state==3))
#roost16 <- as.data.frame(subset(ssf16, state==4))

# extract covariates
proads <- raster("E:/Brood_Covariates/Roads/CC_roads/PrimeRdRast(1).tif")
sroads <- raster("E:/Brood_Covariates/Roads/CC_roads/SecondRdRast(1).tif")
dist_hw <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_HW/dist_hw_21.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_infra <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_infra/dist_infra_21.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_mix <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Mixed/dist_mixed_21.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_open <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Open/dist_open_21.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_pine <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Pine/dist_pine_21.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_shrub <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Shrub/dist_shrub_21.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_water <- raster("E:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Water/dist_water_21.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
ndvi <- raster("E:/Brood_Covariates/NDVI/CC/CC_2020/ndvi_cc_2020.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)


ccssf21$proads <- raster::extract(proads, ccssf21) 
ccssf21$sroads <- raster::extract(sroads, ccssf21)
ccssf21$dist_hw <- raster::extract(dist_hw, ccssf21)
ccssf21$dist_infra <- raster::extract(dist_infra, ccssf21)
ccssf21$dist_mixed <- raster::extract(dist_mix, ccssf21)
ccssf21$dist_open <- raster::extract(dist_open, ccssf21)
ccssf21$dist_pine <- raster::extract(dist_pine, ccssf21)
ccssf21$dist_shrub <- raster::extract(dist_shrub, ccssf21)
ccssf21$dist_water <- raster::extract(dist_water, ccssf21)
ccssf21$ndvi <- raster::extract(ndvi, ccssf21)
ccssf21$cos_ta_ <- cos(ccssf21$ta_)
ccssf21$log_sl_ <- log(ccssf21$sl_)
summary(ccssf21)
#walking16 <- as.data.frame(subset(ssf16, state==1))
#restricted16 <- as.data.frame(subset(ssf16, state==2 | state==3))
#roost16 <- as.data.frame(subset(ssf16, state==4))

############## BF Grant ############################
# All birds
dat <- read.csv("D:/Brood_Recursion/Recurse_State_Master_Sheet.csv")
head(dat)

# One bird
indiv <- subset(dat, site=="BFG")

# Remove birds with no enough data
#indiv <-indiv[!(indiv$ID== 60328),]

# Set the date and time
indiv$datetime <- as.POSIXct(indiv$datetime, format = "%m/%d/%Y %H:%M")
indiv$yr <- format(indiv$datetime, format = "%Y")
head(indiv)

# Subset by year
#indiv <- subset(indiv, yr=="2015")

# Make a spatial data frame
df.sp <- indiv
coordinates(df.sp) <- c("coords.x1","coords.x2")

# Define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB SRS CC 
proj4string(df.sp) <- CRS("+proj=utm +zone=17N ellps=WGS84")

# Transfer the projection
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))

# Add to the data frame
turkey.sp <- as.data.frame(df.sp)
head(turkey.sp)
indiv$easting <- turkey.sp$coords.x1
indiv$northing <- turkey.sp$coords.x2
head(indiv)

# Get data set up for running the loop
ids <- unique(indiv$ID)
memory = random_memory=vector("list",length(ids))
names(memory) <- ids

#Convert to a tibble file
df.t <- as_tibble(indiv)

# Now we can create a track from our cleaned tibble
# Louisiana
#trk <- df.t %>% make_track(easting, northing, datetime, burst_ = ID, crs = CRS("+proj=utm +zone=15N ellps=WGS84"), all_cols = TRUE)
# WEBB CC SRS
trk <- df.t %>% make_track(easting, northing, datetime, burst_ = ID, crs = CRS("+proj=utm +zone=17N ellps=WGS84"), all_cols = TRUE)

head(trk)
trk$burst_ <- trk$ID
summarize_sampling_rate(trk)

# Nest track by ID
trk <- trk %>% nest(data = -"ID")
ssf1 <- trk %>% mutate(steps = map(data, function(x) 
  x %>% steps_by_burst(keep_cols = 'end')))

# Unnest to get complete data frame and look at plot of step length
ssf1 <- ssf1 %>% dplyr::select(ID, steps) %>% tidyr::unnest(cols = steps)

# Get random steps
ssf1 <-ssf1 %>% random_steps(n = 100)
head(ssf1)

# make into a data frame
ran <- as.data.frame(ssf1)

# Make a spatial data frame
df.sp <- ran
head(ran)
coordinates(df.sp) <- c("x2_","y2_")

# define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB, SRS, CC
proj4string(df.sp) <- CRS("+proj=utm +zone=17N ellps=WGS84")
#df.sp$dt_ <- NULL

# transfer the projection
# Louisiana
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))
# WEBB, SC, SRS
df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=17N ellps=WGS84"))

# Set up code to create ID and RD for each individual
for(j in 1:length(ids)){ 
  temp <- indiv[indiv$ID==ids[j],]
  temp2 <- df.sp[df.sp$ID==ids[j],]
  unique(temp$ID)
  unique(temp2$ID)
  # Make a trajectory
  traj <- as.ltraj(temp[,c("easting","northing")],id=temp$ID,date=temp$datetime)
  
  # BRB.D() found average D = 1, Tmax is based on mean residence time, radius 
  # is based on behavioral state, hmin is the location error, max time until left
  # consider it an hour based on GPS unit.
  
  # Memory estimation type ID = residence time, type RD = revisits
  # Behavioral State Walking
  id_rmw <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="ID", radius =200,maxt =3600,grid=200)
  id_rmw1 <- getvolumeUD(id_rmw, standardize = TRUE)
  # Louisiana
  #proj4string(id_rmw1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB SC CC
  proj4string(id_rmw1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  plot(id_rmw1)
  rast1 <- raster(as(id_rmw1,"SpatialPixelsDataFrame"))
  
  rd_rmw <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="RD", radius =200,maxt =3600,grid=200)
  rd_rmw1 <- getvolumeUD(rd_rmw, standardize = TRUE)
  plot(rd_rmw1)
  # Louisiana
  # proj4string(rd_rmw1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(rd_rmw1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  rast2 <- raster(as(rd_rmw1,"SpatialPixelsDataFrame"))
  
  # Write raster files for each individual
  writeRaster(rast1, filename=paste0("D:/Brood_Covariates/Cognitive_Map/ID/BFG", "_", ids[j], "_residence_time_walking.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(rast2, filename=paste0("D:/Brood_Covariates/Cognitive_Map/RD/BFG", "_", ids[j], "_revist_walking.tif"), format="GTiff", overwrite=TRUE)
  
  # Behavioral State Foraging
  id_rmf <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="ID", radius =100,maxt =3600,grid=200)
  id_rmf1 <- getvolumeUD(id_rmf, standardize = TRUE)
  # Louisiana
  #proj4string(id_rmf1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(id_rmf1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  plot(id_rmf1)
  rast3 <- raster(as(id_rmf1,"SpatialPixelsDataFrame"))
  
  rd_rmf <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="RD", radius =100,maxt =3600,grid=200)
  rd_rmf1 <- getvolumeUD(rd_rmf, standardize = TRUE)
  plot(rd_rmf1)
  # Louisiana
  #proj4string(rd_rmf1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(rd_rmf1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  rast4 <- raster(as(rd_rmf1,"SpatialPixelsDataFrame"))
  
  # Write raster files for each individual
  writeRaster(rast3, filename=paste0("D:/Brood_Covariates/Cognitive_Map/ID/BFG", "_", ids[j], "_residence_time_foraging.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(rast4, filename=paste0("D:/Brood_Covariates/Cognitive_Map/RD/BFG", "_", ids[j], "_revist_foraging.tif"), format="GTiff", overwrite=TRUE)
  
  # Behavioral State Roosting
  id_rmr <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="ID", radius =27,maxt =3600,grid=200)
  id_rmr1 <- getvolumeUD(id_rmr, standardize = TRUE)
  # Louisiana
  # proj4string(id_rmr1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(id_rmr1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  plot(id_rmr1)
  rast5 <- raster(as(id_rmr1,"SpatialPixelsDataFrame"))
  
  rd_rmr <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="RD", radius =27,maxt =3600,grid=200)
  rd_rmr1 <- getvolumeUD(rd_rmr, standardize = TRUE)
  plot(rd_rmr1)
  # Louisiana
  #proj4string(rd_rmr1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(rd_rmr1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  rast6 <- raster(as(rd_rmr1,"SpatialPixelsDataFrame"))
  
  # Write raster files for each individual
  writeRaster(rast5, filename=paste0("D:/Brood_Covariates/Cognitive_Map/ID/BFG", "_", ids[j], "_residence_time_roosting.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(rast6, filename=paste0("D:/Brood_Covariates/Cognitive_Map/RD/BFG", "_", ids[j], "_revist_roosting.tif"), format="GTiff", overwrite=TRUE)
  
  # extract values from ID and RD
  # Behavioral state walking
  temp2$residence_raster_walking <- raster::extract(rast1, temp2)
  temp2$revisit_raster_walking <- raster::extract(rast2, temp2)
  
  temp2$residence_raster_foraging <- raster::extract(rast3, temp2)
  temp2$revisit_raster_foraging <- raster::extract(rast4, temp2)
  
  temp2$residence_raster_roosting <- raster::extract(rast5, temp2)
  temp2$revisit_raster_roosting <- raster::extract(rast6, temp2)
  
  temp2 <- as.data.frame(temp2)
  
  assign(paste0(temp2$site,temp2$ID), temp2)
}

df <- rbind( BFG46237, BFG46580, BFG46602, BFG61061, BFG61068, BFG61069,
             BFG61078, BFG61079, BFG61084, BFG61089, BFG61097, BFG61099,
             BFG61108, BFG61111, BFG61113)
head(df)
summary(df)

df <- as.data.frame(df)
head(df)
df$residence_raster_walking[is.na(df$residence_raster_walking)] <- 100
df$residence_raster_foraging[is.na(df$residence_raster_foraging)] <- 100
df$residence_raster_roosting[is.na(df$residence_raster_roosting)] <- 100
df$revisit_raster_walking[is.na(df$revisit_raster_walking)] <- 100
df$revisit_raster_foraging[is.na(df$revisit_raster_foraging)] <- 100
df$revisit_raster_roosting[is.na(df$revisit_raster_roosting)] <- 100
summary(df)

######## Step Selection Function #########
# Make a spatial data frame
df.sp <- df
head(df.sp)
coordinates(df.sp) <- c("x2_","y2_")

# Define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB SRS CC 
proj4string(df.sp) <- CRS("+proj=utm +zone=17N ellps=WGS84")

# Transfer the projection
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))

# subset out years
bfgssf17 <- subset(df.sp, yr==2017)
bfgssf18 <- subset(df.sp, yr==2018)
bfgssf19 <- subset(df.sp, yr==2019)
bfgssf20 <- subset(df.sp, yr==2020)
bfgssf21 <- subset(df.sp, yr==2021)

# extract covariates
proads <- raster("D:/Brood_Covariates/Roads/CC_roads/PrimeRdRast(1).tif")
sroads <- raster("D:/Brood_Covariates/Roads/CC_roads/SecondRdRast(1).tif")
dist_hw <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_HW/dist_hw_17.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_infra <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_infra/dist_infra_17.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_mix <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Mixed/dist_mixed_17.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_open <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Open/dis_open_17.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_pine <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Pine/dist_pine_17.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_shrub <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Shrub/dist_shrub_17.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_water <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Water/dist_water_17.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
ndvi <- raster("D:/Brood_Covariates/NDVI/CC/CC_2017/ndvi_cc_2017.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)


bfgssf17$proads <- raster::extract(proads, bfgssf17) 
bfgssf17$sroads <- raster::extract(sroads, bfgssf17)
bfgssf17$dist_hw <- raster::extract(dist_hw, bfgssf17)
bfgssf17$dist_infra <- raster::extract(dist_infra, bfgssf17)
bfgssf17$dist_mixed <- raster::extract(dist_mix, bfgssf17)
bfgssf17$dist_open <- raster::extract(dist_open, bfgssf17)
bfgssf17$dist_pine <- raster::extract(dist_pine, bfgssf17)
bfgssf17$dist_shrub <- raster::extract(dist_shrub, bfgssf17)
bfgssf17$dist_water <- raster::extract(dist_water, bfgssf17)
bfgssf17$ndvi <- raster::extract(ndvi, bfgssf17)
bfgssf17$cos_ta_ <- cos(bfgssf17$ta_)
bfgssf17$log_sl_ <- log(bfgssf17$sl_)

head(bfgssf17)
#walking15 <- as.data.frame(subset(ssf15, state==3))
#restricted15 <- as.data.frame(subset(ssf15, state==1 | state==2))
#roost15 <- as.data.frame(subset(ssf15, state==4))

# extract covariates
proads <- raster("D:/Brood_Covariates/Roads/CC_roads/PrimeRdRast(1).tif")
sroads <- raster("D:/Brood_Covariates/Roads/CC_roads/SecondRdRast(1).tif")
dist_hw <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_HW/dist_hw_18.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_infra <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_infra/dist_infra_18.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_mix <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Mixed/dist_mixed_18.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_open <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Open/dist_open_18.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_pine <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Pine/dist_pine_18.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_shrub <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Shrub/dist_shrub_18.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_water <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Water/dist_water_18.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
ndvi <- raster("D:/Brood_Covariates/NDVI/CC/CC_2018/ndvi_cc_2018.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)


bfgssf18$proads <- raster::extract(proads, bfgssf18) 
bfgssf18$sroads <- raster::extract(sroads, bfgssf18)
bfgssf18$dist_hw <- raster::extract(dist_hw, bfgssf18)
bfgssf18$dist_infra <- raster::extract(dist_infra, bfgssf18)
bfgssf18$dist_mixed <- raster::extract(dist_mix, bfgssf18)
bfgssf18$dist_open <- raster::extract(dist_open, bfgssf18)
bfgssf18$dist_pine <- raster::extract(dist_pine, bfgssf18)
bfgssf18$dist_shrub <- raster::extract(dist_shrub, bfgssf18)
bfgssf18$dist_water <- raster::extract(dist_water, bfgssf18)
bfgssf18$ndvi <- raster::extract(ndvi, bfgssf18)
bfgssf18$cos_ta_ <- cos(bfgssf18$ta_)
bfgssf18$log_sl_ <- log(bfgssf18$sl_)

#walking16 <- as.data.frame(subset(ssf16, state==1))
#restricted16 <- as.data.frame(subset(ssf16, state==2 | state==3))
#roost16 <- as.data.frame(subset(ssf16, state==4))

# extract covariates
proads <- raster("D:/Brood_Covariates/Roads/CC_roads/PrimeRdRast(1).tif")
sroads <- raster("D:/Brood_Covariates/Roads/CC_roads/SecondRdRast(1).tif")
dist_hw <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_HW/dist_hw_19.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_infra <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_infra/dist_infra_19.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_mix <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Mixed/dist_mixed_19.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_open <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Open/dist_open_19.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_pine <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Pine/dist_pine_19.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_shrub <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Shrub/dist_shrub_19.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_water <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Water/dist_water_19.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
ndvi <- raster("D:/Brood_Covariates/NDVI/CC/CC_2019/ndvi_cc_2019.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)


bfgssf19$proads <- raster::extract(proads, bfgssf19) 
bfgssf19$sroads <- raster::extract(sroads, bfgssf19)
bfgssf19$dist_hw <- raster::extract(dist_hw, bfgssf19)
bfgssf19$dist_infra <- raster::extract(dist_infra, bfgssf19)
bfgssf19$dist_mixed <- raster::extract(dist_mix, bfgssf19)
bfgssf19$dist_open <- raster::extract(dist_open, bfgssf19)
bfgssf19$dist_pine <- raster::extract(dist_pine, bfgssf19)
bfgssf19$dist_shrub <- raster::extract(dist_shrub, bfgssf19)
bfgssf19$dist_water <- raster::extract(dist_water, bfgssf19)
bfgssf19$ndvi <- raster::extract(ndvi, bfgssf19)
bfgssf19$cos_ta_ <- cos(bfgssf19$ta_)
bfgssf19$log_sl_ <- log(bfgssf19$sl_)



#walking16 <- as.data.frame(subset(ssf16, state==1))
#restricted16 <- as.data.frame(subset(ssf16, state==2 | state==3))
#roost16 <- as.data.frame(subset(ssf16, state==4))

# extract covariates
proads <- raster("D:/Brood_Covariates/Roads/CC_roads/PrimeRdRast(1).tif")
sroads <- raster("D:/Brood_Covariates/Roads/CC_roads/SecondRdRast(1).tif")
dist_hw <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_HW/dist_hw_20.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_infra <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_infra/dist_infra_20.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_mix <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Mixed/dist_mixed_20.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_open <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Open/dist_open_20.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_pine <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Pine/dist_pine_20.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_shrub <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Shrub/dist_shrub_20.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_water <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Water/dist_water_20.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
ndvi <- raster("D:/Brood_Covariates/NDVI/CC/CC_2020/ndvi_cc_2020.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)

bfgssf20$proads <- raster::extract(proads, bfgssf20) 
bfgssf20$sroads <- raster::extract(sroads, bfgssf20)
bfgssf20$dist_hw <- raster::extract(dist_hw, bfgssf20)
bfgssf20$dist_infra <- raster::extract(dist_infra, bfgssf20)
bfgssf20$dist_mixed <- raster::extract(dist_mix, bfgssf20)
bfgssf20$dist_open <- raster::extract(dist_open, bfgssf20)
bfgssf20$dist_pine <- raster::extract(dist_pine, bfgssf20)
bfgssf20$dist_shrub <- raster::extract(dist_shrub, bfgssf20)
bfgssf20$dist_water <- raster::extract(dist_water, bfgssf20)
bfgssf20$ndvi <- raster::extract(ndvi, bfgssf20)
bfgssf20$cos_ta_ <- cos(bfgssf20$ta_)
bfgssf20$log_sl_ <- log(bfgssf20$sl_)

dfbfgssf17 <- as.data.frame(bfgssf17)
dfbfgssf18 <- as.data.frame(bfgssf18)
dfbfgssf19 <- as.data.frame(bfgssf19)
dfbfgssf20 <- as.data.frame(bfgssf20)
dfbfgssf21 <- as.data.frame(bfgssf21)

masterbfg <- rbind(dfbfgssf17, dfbfgssf18, dfbfgssf19)
#walking16 <- as.data.frame(subset(ssf16, state==1))
#restricted16 <- as.data.frame(subset(ssf16, state==2 | state==3))
#roost16 <- as.data.frame(subset(ssf16, state==4))

# extract covariates
proads <- raster("D:/Brood_Covariates/Roads/CC_roads/PrimeRdRast(1).tif")
sroads <- raster("D:/Brood_Covariates/Roads/CC_roads/SecondRdRast(1).tif")
dist_hw <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_HW/dist_hw_21.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_infra <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_infra/dist_infra_21.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_mix <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Mixed/dist_mixed_21.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_open <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Open/dist_open_21.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_pine <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Pine/dist_pine_21.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_shrub <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Shrub/dist_shrub_21.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_water <- raster("D:/Brood_Covariates/Landcover/Landcover_CC/Final_Dist_Water/dist_water_21.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
ndvi <- raster("D:/Brood_Covariates/NDVI/CC/CC_2020/ndvi_cc_2020.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)


ccssf21$proads <- raster::extract(proads, ccssf21) 
ccssf21$sroads <- raster::extract(sroads, ccssf21)
ccssf21$dist_hw <- raster::extract(dist_hw, ccssf21)
ccssf21$dist_infra <- raster::extract(dist_infra, ccssf21)
ccssf21$dist_mixed <- raster::extract(dist_mix, ccssf21)
ccssf21$dist_open <- raster::extract(dist_open, ccssf21)
ccssf21$dist_pine <- raster::extract(dist_pine, ccssf21)
ccssf21$dist_shrub <- raster::extract(dist_shrub, ccssf21)
ccssf21$dist_water <- raster::extract(dist_water, ccssf21)
ccssf21$ndvi <- raster::extract(ndvi, ccssf21)
ccssf21$cos_ta_ <- cos(ccssf21$ta_)
ccssf21$log_sl_ <- log(ccssf21$sl_)

dfccssf17 <- as.data.frame(ccssf17)
dfccssf18 <- as.data.frame(ccssf18)
dfccssf19 <- as.data.frame(ccssf19)
dfccssf20 <- as.data.frame(ccssf20)
dfccssf21 <- as.data.frame(ccssf21)
####################################
head(dfccssf17)
head(dfccssf18)
head(dfccssf19)
head(dfccssf21)
mastercc <- rbind(dfccssf17, dfccssf18, dfccssf19, dfccssf20, dfccssf21)
#walking16 <- as.data.frame(subset(ssf16, state==1))
#restricted16 <- as.data.frame(subset(ssf16, state==2 | state==3))
#roost16 <- as.data.frame(subset(ssf16, state==4))

############## Silver Lake WMA ############################
# All birds
dat <- read.csv("D:/Brood_Recursion/Recurse_State_Master_Sheet.csv")
head(dat)
unique(dat$ID)
# One bird
indiv <- subset(dat, site=="SLWMA")
#indiv <- subset(dat, dat$ID==60578)

# Remove birds with no enough data
indiv <-indiv[!(indiv$ID== 60328),]
indiv <-indiv[!(indiv$ID== 60578),]
unique(indiv$ID)
# Set the date and time
indiv$datetime <- as.POSIXct(indiv$datetime, format = "%m/%d/%Y %H:%M")
indiv$yr <- format(indiv$datetime, format = "%Y")
head(indiv)

# Subset by year
#indiv <- subset(indiv, yr=="2015")

# Make a spatial data frame
df.sp <- indiv
coordinates(df.sp) <- c("coords.x1","coords.x2")

# Define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB SRS CC 
proj4string(df.sp) <- CRS("+proj=utm +zone=16N ellps=WGS84")

# Transfer the projection
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))

# Add to the data frame
turkey.sp <- as.data.frame(df.sp)
head(turkey.sp)
indiv$easting <- turkey.sp$coords.x1
indiv$northing <- turkey.sp$coords.x2
head(indiv)

# Get data set up for running the loop
ids <- unique(indiv$ID)
memory = random_memory=vector("list",length(ids))
names(memory) <- ids

#Convert to a tibble file
df.t <- as_tibble(indiv)

# Now we can create a track from our cleaned tibble
# Louisiana
#trk <- df.t %>% make_track(easting, northing, datetime, burst_ = ID, crs = CRS("+proj=utm +zone=15N ellps=WGS84"), all_cols = TRUE)
# WEBB CC SRS
trk <- df.t %>% make_track(easting, northing, datetime, burst_ = ID, crs = CRS("+proj=utm +zone=16N ellps=WGS84"), all_cols = TRUE)

head(trk)
trk$burst_ <- trk$ID
summarize_sampling_rate(trk)

# Nest track by ID
trk <- trk %>% nest(data = -"ID")
ssf1 <- trk %>% mutate(steps = map(data, function(x) 
  x %>% steps_by_burst(keep_cols = 'end')))

# Unnest to get complete data frame and look at plot of step length
ssf1 <- ssf1 %>% dplyr::select(ID, steps) %>% tidyr::unnest(cols = steps)

# Get random steps
ssf1 <-ssf1 %>% random_steps(n = 100)
head(ssf1)

# make into a data frame
ran <- as.data.frame(ssf1)

# Make a spatial data frame
df.sp <- ran
head(ran)
coordinates(df.sp) <- c("x2_","y2_")

# define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB, SRS, CC
proj4string(df.sp) <- CRS("+proj=utm +zone=16N ellps=WGS84")
#df.sp$dt_ <- NULL

# transfer the projection
# Louisiana
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))
# WEBB, SC, SRS
df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=16N ellps=WGS84"))

# Set up code to create ID and RD for each individual
for(j in 1:length(ids)){ 
  temp <- indiv[indiv$ID==ids[j],]
  temp2 <- df.sp[df.sp$ID==ids[j],]
  unique(temp$ID)
  unique(temp2$ID)
  # Make a trajectory
  traj <- as.ltraj(temp[,c("easting","northing")],id=temp$ID,date=temp$datetime)
  
  # BRB.D() found average D = 1, Tmax is based on mean residence time, radius 
  # is based on behavioral state, hmin is the location error, max time until left
  # consider it an hour based on GPS unit.
  
  # Memory estimation type ID = residence time, type RD = revisits
  # Behavioral State Walking
  id_rmw <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="ID", radius =200,maxt =3600,grid=200)
  id_rmw1 <- getvolumeUD(id_rmw, standardize = TRUE)
  # Louisiana
  #proj4string(id_rmw1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB SC CC
  proj4string(id_rmw1) <- CRS("+proj=utm +zone=16N ellps=WGS84")
  plot(id_rmw1)
  rast1 <- raster(as(id_rmw1,"SpatialPixelsDataFrame"))
  
  rd_rmw <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="RD", radius =200,maxt =3600,grid=200)
  rd_rmw1 <- getvolumeUD(rd_rmw, standardize = TRUE)
  plot(rd_rmw1)
  # Louisiana
  # proj4string(rd_rmw1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(rd_rmw1) <- CRS("+proj=utm +zone=16N ellps=WGS84")
  rast2 <- raster(as(rd_rmw1,"SpatialPixelsDataFrame"))
  
  # Write raster files for each individual
  writeRaster(rast1, filename=paste0("D:/Brood_Covariates/Cognitive_Map/ID/SLWMA", "_", ids[j], "_residence_time_walking.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(rast2, filename=paste0("D:/Brood_Covariates/Cognitive_Map/RD/SLWMA", "_", ids[j], "_revist_walking.tif"), format="GTiff", overwrite=TRUE)
  
  # Behavioral State Foraging
  id_rmf <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="ID", radius =100,maxt =3600,grid=200)
  id_rmf1 <- getvolumeUD(id_rmf, standardize = TRUE)
  # Louisiana
  #proj4string(id_rmf1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(id_rmf1) <- CRS("+proj=utm +zone=16N ellps=WGS84")
  plot(id_rmf1)
  rast3 <- raster(as(id_rmf1,"SpatialPixelsDataFrame"))
  
  rd_rmf <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="RD", radius =100,maxt =3600,grid=200)
  rd_rmf1 <- getvolumeUD(rd_rmf, standardize = TRUE)
  plot(rd_rmf1)
  # Louisiana
  #proj4string(rd_rmf1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(rd_rmf1) <- CRS("+proj=utm +zone=16N ellps=WGS84")
  rast4 <- raster(as(rd_rmf1,"SpatialPixelsDataFrame"))
  
  # Write raster files for each individual
  writeRaster(rast3, filename=paste0("D:/Brood_Covariates/Cognitive_Map/ID/SLWMA", "_", ids[j], "_residence_time_foraging.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(rast4, filename=paste0("D:/Brood_Covariates/Cognitive_Map/RD/SLWMA", "_", ids[j], "_revist_foraging.tif"), format="GTiff", overwrite=TRUE)
  
  # Behavioral State Roosting
  id_rmr <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="ID", radius =27,maxt =3600,grid=200)
  id_rmr1 <- getvolumeUD(id_rmr, standardize = TRUE)
  # Louisiana
  # proj4string(id_rmr1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(id_rmr1) <- CRS("+proj=utm +zone=16N ellps=WGS84")
  plot(id_rmr1)
  rast5 <- raster(as(id_rmr1,"SpatialPixelsDataFrame"))
  
  rd_rmr <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="RD", radius =27,maxt =3600,grid=200)
  rd_rmr1 <- getvolumeUD(rd_rmr, standardize = TRUE)
  plot(rd_rmr1)
  # Louisiana
  #proj4string(rd_rmr1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(rd_rmr1) <- CRS("+proj=utm +zone=16N ellps=WGS84")
  rast6 <- raster(as(rd_rmr1,"SpatialPixelsDataFrame"))
  
  # Write raster files for each individual
  writeRaster(rast5, filename=paste0("D:/Brood_Covariates/Cognitive_Map/ID/SLWMA", "_", ids[j], "_residence_time_roosting.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(rast6, filename=paste0("D:/Brood_Covariates/Cognitive_Map/RD/SLWMA", "_", ids[j], "_revist_roosting.tif"), format="GTiff", overwrite=TRUE)
  
  # extract values from ID and RD
  # Behavioral state walking
  temp2$residence_raster_walking <- raster::extract(rast1, temp2)
  temp2$revisit_raster_walking <- raster::extract(rast2, temp2)
  
  temp2$residence_raster_foraging <- raster::extract(rast3, temp2)
  temp2$revisit_raster_foraging <- raster::extract(rast4, temp2)
  
  temp2$residence_raster_roosting <- raster::extract(rast5, temp2)
  temp2$revisit_raster_roosting <- raster::extract(rast6, temp2)
  
  temp2 <- as.data.frame(temp2)
  
  assign(paste0(temp2$site,temp2$ID), temp2)
}

df <- rbind( SLWMA60164, SLWMA60168, SLWMA60170, SLWMA60180, SLWMA60182,
             SLWMA60267, SLWMA60286, SLWMA60288, SLWMA60392, SLWMA60398,
             SLWMA60570, SLWMA60572, SLWMA60574, SLWMA60580, SLWMA60839,
             SLWMA60842, SLWMA60850, SLWMA60852, SLWMA60856, SLWMA70578)
head(df)
summary(df)

df <- as.data.frame(df)
head(df)
df$residence_raster_walking[is.na(df$residence_raster_walking)] <- 100
df$residence_raster_foraging[is.na(df$residence_raster_foraging)] <- 100
df$residence_raster_roosting[is.na(df$residence_raster_roosting)] <- 100
df$revisit_raster_walking[is.na(df$revisit_raster_walking)] <- 100
df$revisit_raster_foraging[is.na(df$revisit_raster_foraging)] <- 100
df$revisit_raster_roosting[is.na(df$revisit_raster_roosting)] <- 100
summary(df)

######## Step Selection Function #########
# Make a spatial data frame
df.sp <- df
head(df.sp)
coordinates(df.sp) <- c("x2_","y2_")

# Define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB SRS CC 
proj4string(df.sp) <- CRS("+proj=utm +zone=16N ellps=WGS84")

# Transfer the projection
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))

# subset out years
slwmassf15 <- subset(df.sp, yr==2015)
slwmassf16 <- subset(df.sp, yr==2016)


# extract covariates
proads <- raster("D:/Brood_Covariates/Roads/SLWMA_roads/prime_rd.tif")
sroads <- raster("D:/Brood_Covariates/Roads/SLWMA_roads/.tif")
dist_hw <- raster("D:/Brood_Covariates/Landcover/Landcover_SLWMA/Final_Dist_HW/dist_hw_15.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=16N ellps=WGS84", res = 30)
dist_infra <- raster("D:/Brood_Covariates/Landcover/Landcover_SLWMA/Final_Dist_infra/dist_infra_15.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=16N ellps=WGS84", res = 30)
dist_mix <- raster("D:/Brood_Covariates/Landcover/Landcover_SLWMA/Final_Dist_Mixed/dist_mix_15.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=16N ellps=WGS84", res = 30)
dist_open <- raster("D:/Brood_Covariates/Landcover/Landcover_SLWMA/Final_Dist_Open/dist_open_15.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=16N ellps=WGS84", res = 30)
dist_pine <- raster("D:/Brood_Covariates/Landcover/Landcover_SLWMA/Final_Dist_Pine/dist_pine_15.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=16N ellps=WGS84", res = 30)
dist_shrub <- raster("D:/Brood_Covariates/Landcover/Landcover_SLWMA/Final_Dist_Shrub/dist_shrub_15.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=16N ellps=WGS84", res = 30)
dist_water <- raster("D:/Brood_Covariates/Landcover/Landcover_SLWMA/Final_Dist_Water/dist_water_15.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=16N ellps=WGS84", res = 30)
ndvi <- raster("D:/Brood_Covariates/NDVI/SLWMA/SLWMA_2015/ndvi_sl_2015.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=16N ellps=WGS84", res = 30)


slwmassf15$proads <- raster::extract(proads, slwmassf15) 
slwmassf15$sroads <- raster::extract(sroads, slwmassf15)
slwmassf15$dist_hw <- raster::extract(dist_hw, slwmassf15)
slwmassf15$dist_infra <- raster::extract(dist_infra, slwmassf15)
slwmassf15$dist_mixed <- raster::extract(dist_mix, slwmassf15)
slwmassf15$dist_open <- raster::extract(dist_open, slwmassf15)
slwmassf15$dist_pine <- raster::extract(dist_pine, slwmassf15)
slwmassf15$dist_shrub <- raster::extract(dist_shrub, slwmassf15)
slwmassf15$dist_water <- raster::extract(dist_water, slwmassf15)
slwmassf15$ndvi <- raster::extract(ndvi, slwmassf15)
slwmassf15$cos_ta_ <- cos(slwmassf15$ta_)
slwmassf15$log_sl_ <- log(slwmassf15$sl_)


#walking15 <- as.data.frame(subset(ssf15, state==3))
#restricted15 <- as.data.frame(subset(ssf15, state==1 | state==2))
#roost15 <- as.data.frame(subset(ssf15, state==4))

# extract covariates
proads <- raster("D:/Brood_Covariates/Roads/SLWMA_roads/prime_rd.tif")
sroads <- raster("D:/Brood_Covariates/Roads/SLWMA_roads/.tif")
dist_hw <- raster("D:/Brood_Covariates/Landcover/Landcover_SLWMA/Final_Dist_HW/dist_hw_16.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=16N ellps=WGS84", res = 30)
dist_infra <- raster("D:/Brood_Covariates/Landcover/Landcover_SLWMA/Final_Dist_infra/dist_infra_16.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=16N ellps=WGS84", res = 30)
dist_mix <- raster("D:/Brood_Covariates/Landcover/Landcover_SLWMA/Final_Dist_Mixed/dist_mix_16.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=16N ellps=WGS84", res = 30)
dist_open <- raster("D:/Brood_Covariates/Landcover/Landcover_SLWMA/Final_Dist_Open/dist_open_16.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=16N ellps=WGS84", res = 30)
dist_pine <- raster("D:/Brood_Covariates/Landcover/Landcover_SLWMA/Final_Dist_Pine/dist_pine_16.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=16N ellps=WGS84", res = 30)
dist_shrub <- raster("D:/Brood_Covariates/Landcover/Landcover_SLWMA/Final_Dist_Shrub/dist_shrub_16.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=16N ellps=WGS84", res = 30)
dist_water <- raster("D:/Brood_Covariates/Landcover/Landcover_SLWMA/Final_Dist_Water/dist_water_16.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=16N ellps=WGS84", res = 30)
ndvi <- raster("D:/Brood_Covariates/NDVI/SLWMA/SLWMA_2016/ndvi_sl_2016.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=16N ellps=WGS84", res = 30)


slwmassf16$proads <- raster::extract(proads, slwmassf16) 
slwmassf16$sroads <- raster::extract(sroads, slwmassf16)
slwmassf16$dist_hw <- raster::extract(dist_hw, slwmassf16)
slwmassf16$dist_infra <- raster::extract(dist_infra, slwmassf16)
slwmassf16$dist_mixed <- raster::extract(dist_mix, slwmassf16)
slwmassf16$dist_open <- raster::extract(dist_open, slwmassf16)
slwmassf16$dist_pine <- raster::extract(dist_pine, slwmassf16)
slwmassf16$dist_shrub <- raster::extract(dist_shrub, slwmassf16)
slwmassf16$dist_water <- raster::extract(dist_water, slwmassf16)
slwmassf16$ndvi <- raster::extract(ndvi, slwmassf16)
slwmassf16$cos_ta_ <- cos(slwmassf16$ta_)
slwmassf16$log_sl_ <- log(slwmassf16$sl_)

#walking16 <- as.data.frame(subset(ssf16, state==1))
#restricted16 <- as.data.frame(subset(ssf16, state==2 | state==3))
#roost16 <- as.data.frame(subset(ssf16, state==4))
dfslwma15 <- as.data.frame(slwmassf15)
dfslwma16 <- as.data.frame(slwmassf16)
masterslwma <- rbind(dfslwma15,dfslwma16)

############## Kistachie National Forest ############################
# All birds
dat <- read.csv("E:/Brood_Recursion/Recurse_State_Master_Sheet.csv")
unique(dat$ID)

# One bird
indiv <- subset(dat, site=="KNF")
#indiv <- indiv[!(indiv$ID==60578),]

# Remove birds with no enough data
#indiv <-indiv[!(indiv$ID== 60328),]

# Set the date and time
indiv$datetime <- as.POSIXct(indiv$datetime, format = "%m/%d/%Y %H:%M")
indiv$yr <- format(indiv$datetime, format = "%Y")
head(indiv)

# Subset by year
#indiv <- subset(indiv, yr=="2015")

# Make a spatial data frame
df.sp <- indiv
coordinates(df.sp) <- c("coords.x1","coords.x2")

# Define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB SRS CC 
proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")

# Transfer the projection
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))

# Add to the data frame
turkey.sp <- as.data.frame(df.sp)
head(turkey.sp)
indiv$easting <- turkey.sp$coords.x1
indiv$northing <- turkey.sp$coords.x2
head(indiv)

# Get data set up for running the loop
ids <- unique(indiv$ID)
memory = random_memory=vector("list",length(ids))
names(memory) <- ids

#Convert to a tibble file
df.t <- as_tibble(indiv)

# Now we can create a track from our cleaned tibble
# Louisiana
#trk <- df.t %>% make_track(easting, northing, datetime, burst_ = ID, crs = CRS("+proj=utm +zone=15N ellps=WGS84"), all_cols = TRUE)
# WEBB CC SRS
trk <- df.t %>% make_track(easting, northing, datetime, burst_ = ID, crs = CRS("+proj=utm +zone=15N ellps=WGS84"), all_cols = TRUE)

head(trk)
trk$burst_ <- trk$ID
summarize_sampling_rate(trk)

# Nest track by ID
trk <- trk %>% nest(data = -"ID")
ssf1 <- trk %>% mutate(steps = map(data, function(x) 
  x %>% steps_by_burst(keep_cols = 'end')))

# Unnest to get complete data frame and look at plot of step length
ssf1 <- ssf1 %>% dplyr::select(ID, steps) %>% tidyr::unnest(cols = steps)

# Get random steps
ssf1 <-ssf1 %>% random_steps(n = 50)
head(ssf1)

# make into a data frame
ran <- as.data.frame(ssf1)

# Make a spatial data frame
df.sp <- ran
head(ran)
coordinates(df.sp) <- c("x2_","y2_")

# define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB, SRS, CC
proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
#df.sp$dt_ <- NULL

# transfer the projection
# Louisiana
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))
# WEBB, SC, SRS
df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))

# Set up code to create ID and RD for each individual
for(j in 1:length(ids)){ 
  temp <- indiv[indiv$ID==ids[j],]
  temp2 <- df.sp[df.sp$ID==ids[j],]
  unique(temp$ID)
  unique(temp2$ID)
  # Make a trajectory
  traj <- as.ltraj(temp[,c("easting","northing")],id=temp$ID,date=temp$datetime)
  
  # BRB.D() found average D = 1, Tmax is based on mean residence time, radius 
  # is based on behavioral state, hmin is the location error, max time until left
  # consider it an hour based on GPS unit.
  
  # Memory estimation type ID = residence time, type RD = revisits
  # Behavioral State Walking
  id_rmw <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="ID", radius =200,maxt =3600,grid=200)
  id_rmw1 <- getvolumeUD(id_rmw, standardize = TRUE)
  # Louisiana
  #proj4string(id_rmw1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB SC CC
  proj4string(id_rmw1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  plot(id_rmw1)
  rast1 <- raster(as(id_rmw1,"SpatialPixelsDataFrame"))
  
  rd_rmw <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="RD", radius =200,maxt =3600,grid=200)
  rd_rmw1 <- getvolumeUD(rd_rmw, standardize = TRUE)
  plot(rd_rmw1)
  # Louisiana
  # proj4string(rd_rmw1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(rd_rmw1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  rast2 <- raster(as(rd_rmw1,"SpatialPixelsDataFrame"))
  
  # Write raster files for each individual
  writeRaster(rast1, filename=paste0("D:/Brood_Covariates/Cognitive_Map/ID/KNF", "_", ids[j], "_residence_time_walking.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(rast2, filename=paste0("D:/Brood_Covariates/Cognitive_Map/RD/KNF", "_", ids[j], "_revist_walking.tif"), format="GTiff", overwrite=TRUE)
  
  # Behavioral State Foraging
  id_rmf <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="ID", radius =100,maxt =3600,grid=200)
  id_rmf1 <- getvolumeUD(id_rmf, standardize = TRUE)
  # Louisiana
  #proj4string(id_rmf1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(id_rmf1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  plot(id_rmf1)
  rast3 <- raster(as(id_rmf1,"SpatialPixelsDataFrame"))
  
  rd_rmf <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="RD", radius =100,maxt =3600,grid=200)
  rd_rmf1 <- getvolumeUD(rd_rmf, standardize = TRUE)
  plot(rd_rmf1)
  # Louisiana
  #proj4string(rd_rmf1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(rd_rmf1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  rast4 <- raster(as(rd_rmf1,"SpatialPixelsDataFrame"))
  
  # Write raster files for each individual
  writeRaster(rast3, filename=paste0("D:/Brood_Covariates/Cognitive_Map/ID/KNF", "_", ids[j], "_residence_time_foraging.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(rast4, filename=paste0("D:/Brood_Covariates/Cognitive_Map/RD/KNF", "_", ids[j], "_revist_foraging.tif"), format="GTiff", overwrite=TRUE)
  
  # Behavioral State Roosting
  id_rmr <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="ID", radius =27,maxt =3600,grid=200)
  id_rmr1 <- getvolumeUD(id_rmr, standardize = TRUE)
  # Louisiana
  # proj4string(id_rmr1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(id_rmr1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  plot(id_rmr1)
  rast5 <- raster(as(id_rmr1,"SpatialPixelsDataFrame"))
  
  rd_rmr <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="RD", radius =27,maxt =3600,grid=200)
  rd_rmr1 <- getvolumeUD(rd_rmr, standardize = TRUE)
  plot(rd_rmr1)
  # Louisiana
  #proj4string(rd_rmr1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(rd_rmr1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  rast6 <- raster(as(rd_rmr1,"SpatialPixelsDataFrame"))
  
  # Write raster files for each individual
  writeRaster(rast5, filename=paste0("D:/Brood_Covariates/Cognitive_Map/ID/KNF", "_", ids[j], "_residence_time_roosting.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(rast6, filename=paste0("D:/Brood_Covariates/Cognitive_Map/RD/KNF", "_", ids[j], "_revist_roosting.tif"), format="GTiff", overwrite=TRUE)
  
  # extract values from ID and RD
  # Behavioral state walking
  temp2$residence_raster_walking <- raster::extract(rast1, temp2)
  temp2$revisit_raster_walking <- raster::extract(rast2, temp2)
  
  temp2$residence_raster_foraging <- raster::extract(rast3, temp2)
  temp2$revisit_raster_foraging <- raster::extract(rast4, temp2)
  
  temp2$residence_raster_roosting <- raster::extract(rast5, temp2)
  temp2$revisit_raster_roosting <- raster::extract(rast6, temp2)
  
  temp2 <- as.data.frame(temp2)
  
  assign(paste0(temp2$site,temp2$ID), temp2)
}
#### Where I stopped!####################################
df <- rbind(KNF45, KNF46193, KNF46196, KNF46473, KNF46474, KNF46475, KNF46476,
            KNF46739, KNF46890, KNF46905, KNF46927, KNF46938, KNF46939, KNF46941,
            KNF46946, KNF46962, KNF47366, KNF59, KNF60273, KNF60327, KNF60328,KNF60338,
            KNF60343, KNF60345, KNF60346, KNF61121, KNF61154, KNF61174, KNF70327)
head(df)
summary(df)

df <- as.data.frame(df)
head(df)
df$residence_raster_walking[is.na(df$residence_raster_walking)] <- 100
df$residence_raster_foraging[is.na(df$residence_raster_foraging)] <- 100
df$residence_raster_roosting[is.na(df$residence_raster_roosting)] <- 100
df$revisit_raster_walking[is.na(df$revisit_raster_walking)] <- 100
df$revisit_raster_foraging[is.na(df$revisit_raster_foraging)] <- 100
df$revisit_raster_roosting[is.na(df$revisit_raster_roosting)] <- 100
summary(df)

######## Step Selection Function #########
df <- read.csv("C:/Users/nwb74172/Desktop/knf_recurse_covariates.csv")
df$X.2 <- NULL
# Make a spatial data frame
df.sp <- df
head(df.sp)
coordinates(df.sp) <- c("x2_","y2_")

# Define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB SRS CC 
proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")

# Transfer the projection
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))

# subset out years
knfssf14 <- subset(df.sp, yr==2014)
knfssf15 <- subset(df.sp, yr==2015)
knfssf18 <- subset(df.sp, yr==2018)
knfssf19 <- subset(df.sp, yr==2019)
knfssf20 <- subset(df.sp, yr==2020)
knfssf21 <- subset(df.sp, yr==2021)

# extract covariates
proads <- raster("E:/Brood_Covariates/Roads/KNF_roads/Use_Prime.tif")
proads <- projectRaster(proads, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
sroads <- raster("E:/Brood_Covariates/Roads/KNF_roads/Use_Secondary.tif")
sroads <- projectRaster(sroads, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_hw <- raster("E:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_HW/dist_hw_14.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_infra <- raster("E:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_infra/dist_infra_14.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_mix <- raster("E:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Mixed/dist_mix_14.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_open <- raster("E:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Open/dist_open_14.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_pine <- raster("E:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Pine/dist_pine_14.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_shrub <- raster("E:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Shrub/dist_shrub_14.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_water <- raster("E:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Water/dist_water_14.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
ndvi <- raster("E:/Brood_Covariates/NDVI/KNF/KNF_2014/ndvi_knf_2014.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)


knfssf14$proads <- raster::extract(proads, knfssf14) 
knfssf14$sroads <- raster::extract(sroads, knfssf14)
knfssf14$dist_hw <- raster::extract(dist_hw, knfssf14)
knfssf14$dist_infra <- raster::extract(dist_infra, knfssf14)
knfssf14$dist_mixed <- raster::extract(dist_mix, knfssf14)
knfssf14$dist_open <- raster::extract(dist_open, knfssf14)
knfssf14$dist_pine <- raster::extract(dist_pine, knfssf14)
knfssf14$dist_shrub <- raster::extract(dist_shrub, knfssf14)
knfssf14$dist_water <- raster::extract(dist_water, knfssf14)
knfssf14$ndvi <- raster::extract(ndvi, knfssf14)
knfssf14$cos_ta_ <- cos(knfssf14$ta_)
knfssf14$log_sl_ <- log(knfssf14$sl_)
summary(knfssf14)

#walking15 <- as.data.frame(subset(ssf15, state==3))
#restricted15 <- as.data.frame(subset(ssf15, state==1 | state==2))
#roost15 <- as.data.frame(subset(ssf15, state==4))

# extract covariates
proads <- raster("D:/Brood_Covariates/Roads/KNF_roads/Use_Prime.tif")
proads <- projectRaster(proads, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
sroads <- raster("D:/Brood_Covariates/Roads/KNF_roads/Use_Secondary.tif")
sroads <- projectRaster(sroads, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_hw <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_HW/dist_hw_15.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_infra <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_infra/dist_infra_15.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_mix <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Mixed/dist_mix_15.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_open <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Open/dist_open_15.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_pine <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Pine/dist_pine_15.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_shrub <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Shrub/dist_shrub_15.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_water <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Water/dist_water_15.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
ndvi <- raster("D:/Brood_Covariates/NDVI/KNF/KNF_2015/ndvi_knf_2015.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)


knfssf15$proads <- raster::extract(proads, knfssf15) 
knfssf15$sroads <- raster::extract(sroads, knfssf15)
knfssf15$dist_hw <- raster::extract(dist_hw, knfssf15)
knfssf15$dist_infra <- raster::extract(dist_infra, knfssf15)
knfssf15$dist_mixed <- raster::extract(dist_mix, knfssf15)
knfssf15$dist_open <- raster::extract(dist_open, knfssf15)
knfssf15$dist_pine <- raster::extract(dist_pine, knfssf15)
knfssf15$dist_shrub <- raster::extract(dist_shrub, knfssf15)
knfssf15$dist_water <- raster::extract(dist_water, knfssf15)
knfssf15$ndvi <- raster::extract(ndvi, knfssf15)
knfssf15$cos_ta_ <- cos(knfssf15$ta_)
knfssf15$log_sl_ <- log(knfssf15$sl_)
####### Fix distance landcover ###################
summary(knfssf15)
#walking16 <- as.data.frame(subset(ssf16, state==1))
#restricted16 <- as.data.frame(subset(ssf16, state==2 | state==3))
#roost16 <- as.data.frame(subset(ssf16, state==4))

# extract covariates
proads <- raster("D:/Brood_Covariates/Roads/KNF_roads/Use_Prime.tif")
proads <- projectRaster(proads, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
sroads <- raster("D:/Brood_Covariates/Roads/KNF_roads/Use_Secondary.tif")
sroads <- projectRaster(sroads, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_hw <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_HW/dist_hw_18.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_infra <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_infra/dist_infra_18.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_mix <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Mixed/dist_mix_18.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_open <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Open/dist_open_18.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_pine <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Pine/dist_pine_18.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_shrub <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Shrub/dist_shrub_18.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_water <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Water/dist_water_18.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
ndvi <- raster("D:/Brood_Covariates/NDVI/KNF/KNF_2018/ndvi_knf_2018.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)



knfssf18$proads <- raster::extract(proads, knfssf18) 
knfssf18$sroads <- raster::extract(sroads, knfssf18)
knfssf18$dist_hw <- raster::extract(dist_hw, knfssf18)
knfssf18$dist_infra <- raster::extract(dist_infra, knfssf18)
knfssf18$dist_mixed <- raster::extract(dist_mix, knfssf18)
knfssf18$dist_open <- raster::extract(dist_open, knfssf18)
knfssf18$dist_pine <- raster::extract(dist_pine, knfssf18)
knfssf18$dist_shrub <- raster::extract(dist_shrub, knfssf18)
knfssf18$dist_water <- raster::extract(dist_water, knfssf18)
knfssf18$ndvi <- raster::extract(ndvi, knfssf18)
knfssf18$cos_ta_ <- cos(knfssf18$ta_)
knfssf18$log_sl_ <- log(knfssf18$sl_)
summary(knfssf18)

#### Fix landcover maps#######

#walking16 <- as.data.frame(subset(ssf16, state==1))
#restricted16 <- as.data.frame(subset(ssf16, state==2 | state==3))
#roost16 <- as.data.frame(subset(ssf16, state==4))

# extract covariates
proads <- raster("D:/Brood_Covariates/Roads/KNF_roads/Use_Prime.tif")
proads <- projectRaster(proads, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
sroads <- raster("D:/Brood_Covariates/Roads/KNF_roads/Use_Secondary.tif")
sroads <- projectRaster(sroads, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_hw <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_HW/dist_hw_19.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_infra <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_infra/dist_infra_19.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_mix <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Mixed/dist_mix_19.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_open <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Open/dist_open_19.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_pine <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Pine/dist_pine_19.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_shrub <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Shrub/dist_shrub_19.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_water <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Water/dist_water_19.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
ndvi <- raster("D:/Brood_Covariates/NDVI/KNF/KNF_2019/ndvi_knf_2019.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)


knfssf19$proads <- raster::extract(proads, knfssf19) 
knfssf19$sroads <- raster::extract(sroads, knfssf19)
knfssf19$dist_hw <- raster::extract(dist_hw, knfssf19)
knfssf19$dist_infra <- raster::extract(dist_infra, knfssf19)
knfssf19$dist_mixed <- raster::extract(dist_mix, knfssf19)
knfssf19$dist_open <- raster::extract(dist_open, knfssf19)
knfssf19$dist_pine <- raster::extract(dist_pine, knfssf19)
knfssf19$dist_shrub <- raster::extract(dist_shrub, knfssf19)
knfssf19$dist_water <- raster::extract(dist_water, knfssf19)
knfssf19$ndvi <- raster::extract(ndvi, knfssf19)
knfssf19$cos_ta_ <- cos(knfssf19$ta_)
knfssf19$log_sl_ <- log(knfssf19$sl_)


#walking16 <- as.data.frame(subset(ssf16, state==1))
#restricted16 <- as.data.frame(subset(ssf16, state==2 | state==3))
#roost16 <- as.data.frame(subset(ssf16, state==4))

# extract covariates
proads <- raster("D:/Brood_Covariates/Roads/KNF_roads/Use_Prime.tif")
proads <- projectRaster(proads, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
sroads <- raster("D:/Brood_Covariates/Roads/KNF_roads/Use_Secondary.tif")
sroads <- projectRaster(sroads, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_hw <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_HW/dist_hw_20.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_infra <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_infra/dist_infra_20.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_mix <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Mixed/dist_mix_20.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_open <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Open/dist_open_20.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_pine <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Pine/dist_pine_20.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_shrub <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Shrub/dist_shrub_20.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_water <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Water/dist_water_20.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
ndvi <- raster("D:/Brood_Covariates/NDVI/KNF/KNF_2020/ndvi_knf_2020.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)


knfssf20$proads <- raster::extract(proads, knfssf20) 
knfssf20$sroads <- raster::extract(sroads, knfssf20)
knfssf20$dist_hw <- raster::extract(dist_hw, knfssf20)
knfssf20$dist_infra <- raster::extract(dist_infra, knfssf20)
knfssf20$dist_mixed <- raster::extract(dist_mix, knfssf20)
knfssf20$dist_open <- raster::extract(dist_open, knfssf20)
knfssf20$dist_pine <- raster::extract(dist_pine, knfssf20)
knfssf20$dist_shrub <- raster::extract(dist_shrub, knfssf20)
knfssf20$dist_water <- raster::extract(dist_water, knfssf20)
knfssf20$ndvi <- raster::extract(ndvi, knfssf20)
knfssf20$cos_ta_ <- cos(knfssf20$ta_)
knfssf20$log_sl_ <- log(knfssf20$sl_)
summary(knfssf20)
#walking16 <- as.data.frame(subset(ssf16, state==1))
#restricted16 <- as.data.frame(subset(ssf16, state==2 | state==3))
#roost16 <- as.data.frame(subset(ssf16, state==4))

# extract covariates
proads <- raster("D:/Brood_Covariates/Roads/KNF_roads/Use_Prime.tif")
proads <- projectRaster(proads, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
sroads <- raster("D:/Brood_Covariates/Roads/KNF_roads/Use_Secondary.tif")
sroads <- projectRaster(sroads, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_hw <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_HW/dist_hw_21.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_infra <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_infra/dist_infra_21.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_mix <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Mixed/dist_mix_21.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_open <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Open/dist_open_21.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_pine <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Pine/dist_pine_21.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_shrub <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Shrub/dist_shrub_21.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
dist_water <- raster("D:/Brood_Covariates/Landcover/Landcover_KNF/Final_Dist_Water/dist_water_21.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)
ndvi <- raster("D:/Brood_Covariates/NDVI/KNF/KNF_2021/ndvi_knf_2021.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=15N ellps=WGS84", res = 30)


knfssf21$proads <- raster::extract(proads, knfssf21) 
knfssf21$sroads <- raster::extract(sroads, knfssf21)
knfssf21$dist_hw <- raster::extract(dist_hw, knfssf21)
knfssf21$dist_infra <- raster::extract(dist_infra, knfssf21)
knfssf21$dist_mixed <- raster::extract(dist_mix, knfssf21)
knfssf21$dist_open <- raster::extract(dist_open, knfssf21)
knfssf21$dist_pine <- raster::extract(dist_pine, knfssf21)
knfssf21$dist_shrub <- raster::extract(dist_shrub, knfssf21)
knfssf21$dist_water <- raster::extract(dist_water, knfssf21)
knfssf21$ndvi <- raster::extract(ndvi, knfssf21)
knfssf21$cos_ta_ <- cos(knfssf21$ta_)
knfssf21$log_sl_ <- log(knfssf21$sl_)
summary(knfssf21)
masterknf <- masterknf[masterknf$yr != "2015", ]
unique(masterknf$yr)
masterknf <- rbind(knfssf14, knfssf15, knfssf18, knfssf19, knfssf20, knfssf21)

knfssf15 <- as.data.frame(knfssf15)
masterknf <- rbind(masterknf, knfssf15)
masterknf <- as.data.frame(masterknf)

###### knf 2018 and 2015 need landcover done ###

masterwebb <- rbind(ssf15, ssf16, ssf18)
masterwebb <- as.data.frame(masterwebb)


############## SRS ############################
# All birds
dat <- read.csv("E:/Brood_Recursion/Recurse_State_Master_Sheet.csv")
head(dat)

# One bird
indiv <- subset(dat, site=="SRS")
#indiv <- indiv[!(indiv$ID==60578),]

# Remove birds with no enough data
#indiv <-indiv[!(indiv$ID== 60328),]

# Set the date and time
indiv$datetime <- as.POSIXct(indiv$datetime, format = "%m/%d/%Y %H:%M")
indiv$yr <- format(indiv$datetime, format = "%Y")
head(indiv)

# Subset by year
#indiv <- subset(indiv, yr=="2015")

# Make a spatial data frame
df.sp <- indiv
coordinates(df.sp) <- c("coords.x1","coords.x2")

# Define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB SRS CC 
proj4string(df.sp) <- CRS("+proj=utm +zone=17N ellps=WGS84")

# Transfer the projection
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))

# Add to the data frame
turkey.sp <- as.data.frame(df.sp)
head(turkey.sp)
indiv$easting <- turkey.sp$coords.x1
indiv$northing <- turkey.sp$coords.x2
head(indiv)

# Get data set up for running the loop
ids <- unique(indiv$ID)
memory = random_memory=vector("list",length(ids))
names(memory) <- ids

#Convert to a tibble file
df.t <- as_tibble(indiv)

# Now we can create a track from our cleaned tibble
# Louisiana
#trk <- df.t %>% make_track(easting, northing, datetime, burst_ = ID, crs = CRS("+proj=utm +zone=15N ellps=WGS84"), all_cols = TRUE)
# WEBB CC SRS
trk <- df.t %>% make_track(easting, northing, datetime, burst_ = ID, crs = CRS("+proj=utm +zone=17N ellps=WGS84"), all_cols = TRUE)

head(trk)
trk$burst_ <- trk$ID
summarize_sampling_rate(trk)

# Nest track by ID
trk <- trk %>% nest(data = -"ID")
ssf1 <- trk %>% mutate(steps = map(data, function(x) 
  x %>% steps_by_burst(keep_cols = 'end')))

# Unnest to get complete data frame and look at plot of step length
ssf1 <- ssf1 %>% dplyr::select(ID, steps) %>% tidyr::unnest(cols = steps)

# Get random steps
ssf1 <-ssf1 %>% random_steps(n = 100)
head(ssf1)

# make into a data frame
ran <- as.data.frame(ssf1)

# Make a spatial data frame
df.sp <- ran
head(ran)
coordinates(df.sp) <- c("x2_","y2_")

# define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB, SRS, CC
proj4string(df.sp) <- CRS("+proj=utm +zone=17N ellps=WGS84")
#df.sp$dt_ <- NULL

# transfer the projection
# Louisiana
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))
# WEBB, SC, SRS
df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=17N ellps=WGS84"))

# Set up code to create ID and RD for each individual
for(j in 1:length(ids)){ 
  temp <- indiv[indiv$ID==ids[j],]
  temp2 <- df.sp[df.sp$ID==ids[j],]
  unique(temp$ID)
  unique(temp2$ID)
  # Make a trajectory
  traj <- as.ltraj(temp[,c("easting","northing")],id=temp$ID,date=temp$datetime)
  
  # BRB.D() found average D = 1, Tmax is based on mean residence time, radius 
  # is based on behavioral state, hmin is the location error, max time until left
  # consider it an hour based on GPS unit.
  
  # Memory estimation type ID = residence time, type RD = revisits
  # Behavioral State Walking
  id_rmw <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="ID", radius =250,maxt =3600,grid=200)
  id_rmw1 <- getvolumeUD(id_rmw, standardize = TRUE)
  # Louisiana
  #proj4string(id_rmw1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB SC CC
  proj4string(id_rmw1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  plot(id_rmw1)
  rast1 <- raster(as(id_rmw1,"SpatialPixelsDataFrame"))
  
  rd_rmw <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="RD", radius =250,maxt =3600,grid=200)
  rd_rmw1 <- getvolumeUD(rd_rmw, standardize = TRUE)
  plot(rd_rmw1)
  # Louisiana
  # proj4string(rd_rmw1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(rd_rmw1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  rast2 <- raster(as(rd_rmw1,"SpatialPixelsDataFrame"))
  
  # Write raster files for each individual
  writeRaster(rast1, filename=paste0("E:/Brood_Covariates/Cognitive_Map/ID/SRS", "_", ids[j], "_residence_time_walking.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(rast2, filename=paste0("E:/Brood_Covariates/Cognitive_Map/RD/SRS", "_", ids[j], "_revist_walking.tif"), format="GTiff", overwrite=TRUE)
  
  # Behavioral State Foraging
  id_rmf <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="ID", radius =90,maxt =3600,grid=200)
  id_rmf1 <- getvolumeUD(id_rmf, standardize = TRUE)
  # Louisiana
  #proj4string(id_rmf1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(id_rmf1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  plot(id_rmf1)
  rast3 <- raster(as(id_rmf1,"SpatialPixelsDataFrame"))
  
  rd_rmf <- BRB(traj, D=1, Tmax=139968, Lmin=15, hmin=27, type="RD", radius =90,maxt =3600,grid=200)
  rd_rmf1 <- getvolumeUD(rd_rmf, standardize = TRUE)
  plot(rd_rmf1)
  # Louisiana
  #proj4string(rd_rmf1) <- CRS("+proj=utm +zone=15N ellps=WGS84")
  # WEBB CC SRS
  proj4string(rd_rmf1) <- CRS("+proj=utm +zone=17N ellps=WGS84")
  rast4 <- raster(as(rd_rmf1,"SpatialPixelsDataFrame"))
  
  # Write raster files for each individual
  writeRaster(rast3, filename=paste0("E:/Brood_Covariates/Cognitive_Map/ID/SRS", "_", ids[j], "_residence_time_foraging.tif"), format="GTiff", overwrite=TRUE)
  writeRaster(rast4, filename=paste0("E:/Brood_Covariates/Cognitive_Map/RD/SRS", "_", ids[j], "_revist_foraging.tif"), format="GTiff", overwrite=TRUE)
  
  # extract values from ID and RD
  # Behavioral state walking
  temp2$residence_raster_walking <- raster::extract(rast1, temp2)
  temp2$revisit_raster_walking <- raster::extract(rast2, temp2)
  
  temp2$residence_raster_foraging <- raster::extract(rast3, temp2)
  temp2$revisit_raster_foraging <- raster::extract(rast4, temp2)
  
  temp2 <- as.data.frame(temp2)
  
  assign(paste0(temp2$site,temp2$ID), temp2)
}

df <- rbind(SRS47387, SRS47424, SRS47425)
head(df)
summary(df)

df <- as.data.frame(df)
head(df)
df$residence_raster_walking[is.na(df$residence_raster_walking)] <- 100
df$residence_raster_foraging[is.na(df$residence_raster_foraging)] <- 100

df$revisit_raster_walking[is.na(df$revisit_raster_walking)] <- 100
df$revisit_raster_foraging[is.na(df$revisit_raster_foraging)] <- 100

summary(df)

######## Step Selection Function #########
# Make a spatial data frame
df.sp <- df
head(df.sp)
coordinates(df.sp) <- c("x2_","y2_")

# Define crs for df.sp
# Louisiana
#proj4string(df.sp) <- CRS("+proj=utm +zone=15N ellps=WGS84")
# WEBB SRS CC 
proj4string(df.sp) <- CRS("+proj=utm +zone=17N ellps=WGS84")

# Transfer the projection
#df.sp <- spTransform(df.sp, CRS("+proj=utm +zone=15N ellps=WGS84"))

# subset out years
srsssf21 <- subset(df.sp, yr==2021)

# extract covariates
proads <- raster("E:/Brood_Covariates/Roads/SRS_roads/primary_rd1.tif")
proads <- projectRaster(proads, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
sroads <- raster("E:/Brood_Covariates/Roads/SRS_roads/secondary_rd1.tif")
sroads <- projectRaster(sroads, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_hw <- raster("E:/Brood_Covariates/Landcover/Landcover_SRS/Final_Dist_HW/dist_hw_21.tif")
dist_hw <- projectRaster(dist_hw, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_infra <- raster("E:/Brood_Covariates/Landcover/Landcover_SRS/Final_Dist_infra/dist_infra_21.tif")
dist_infra <- projectRaster(dist_infra, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_mix <- raster("E:/Brood_Covariates/Landcover/Landcover_SRS/Final_Dist_Mixed/dist_mix_21.tif")
dist_mix <- projectRaster(dist_mix, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_open <- raster("E:/Brood_Covariates/Landcover/Landcover_SRS/Final_Dist_Open/dist_open_21.tif")
dist_open <- projectRaster(dist_open, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_pine <- raster("E:/Brood_Covariates/Landcover/Landcover_SRS/Final_Dist_Pine/dist_pine_21.tif")
dist_pine <- projectRaster(dist_pine, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_shrub <- raster("E:/Brood_Covariates/Landcover/Landcover_SRS/Final_Dist_Shrub/dist_shrub_21.tif")
dist_shrub <- projectRaster(dist_shrub, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
dist_water <- raster("E:/Brood_Covariates/Landcover/Landcover_SRS/Final_Dist_Water/dist_water_21.tif")
dist_water <- projectRaster(dist_water, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)
ndvi <- raster("E:/Brood_Covariates/NDVI/SRS/SRS_2021/ndvi_srs_2021.tif")
ndvi <- projectRaster(ndvi, crs="+proj=utm +zone=17N ellps=WGS84", res = 30)


srsssf21$proads <- raster::extract(proads, srsssf21) 
srsssf21$sroads <- raster::extract(sroads, srsssf21)
srsssf21$dist_hw <- raster::extract(dist_hw, srsssf21)
srsssf21$dist_infra <- raster::extract(dist_infra, srsssf21)
srsssf21$dist_mixed <- raster::extract(dist_mix, srsssf21)
srsssf21$dist_open <- raster::extract(dist_open, srsssf21)
srsssf21$dist_pine <- raster::extract(dist_pine, srsssf21)
srsssf21$dist_shrub <- raster::extract(dist_shrub, srsssf21)
srsssf21$dist_water <- raster::extract(dist_water, srsssf21)
srsssf21$ndvi <- raster::extract(ndvi, srsssf21)
srsssf21$cos_ta_ <- cos(srsssf21$ta_)
srsssf21$log_sl_ <- log(srsssf21$sl_)
summary(srsssf21)
plot(dist_water)
plot(srsssf21, add=TRUE)

mastersrs <- as.data.frame(srsssf21)

##### Create master sheet of covariates #####
masterknf <- read.csv("C:/Users/nwb74172/Desktop/knf_covariates_final_use.csv")
masterknf$X.2 <- NULL
masterknf$direction_p <- NULL
str(masterknf)
str(masterslwma)
str(masterwebb)

masterknf$t1_<- as.POSIXct(masterknf$t1_, format = "%m/%d/%Y %H:%M")
masterslwma$t1_<- as.POSIXct(masterslwma$t1_, format = "%m/%d/%Y %H:%M")
masterwebb$t1_<- as.POSIXct(masterwebb$t1_, format = "%m/%d/%Y %H:%M")

masterknf$t2_<- as.POSIXct(masterknf$t2_, format = "%m/%d/%Y %H:%M")
masterslwma$t2_<- as.POSIXct(masterslwma$t2_, format = "%m/%d/%Y %H:%M")
masterwebb$t2_<- as.POSIXct(masterwebb$t2_, format = "%m/%d/%Y %H:%M")

all_birds_covariates <- rbind(mastercc, masterslwma, masterknf, masterwebb, mastersrs)

dat <- dat[order(dat$ID),]
unique(dat$ID)

all_birds_covariates <- all_birds_covariates[order(all_birds_covariates$ID),]
unique(all_birds_covariates$ID)

write.csv(all_birds_covariates, "E:/Brood_Recursion/covariates_master_sheet_100R.csv")
summary(all_birds_covariates)

#masterslwma <- read.csv("C:/Users/nwb74172/Desktop/Covariate_SLWMA.csv")
#masterslwma$X.2 <- NULL
#head(masterslwma)

#masterwebb<- read.csv("C:/Users/nwb74172/Desktop/covariates_master_sheet_100R.csv")
#masterwebb$X.2 <- NULL
#head(masterwebb)

#all_birds_covariates <- rbind(masterslwma, masterwebb)
#write.csv(all_birds_covariates, "C:/Users/nwb74172/Desktop/No_KNF.csv")
############### Begin the model ##############################

# Subset to different states
all_birds_covariates100 <- read.csv("D:/Brood_Recursion/covariates_master_sheet_100R.csv")
#unique(all_birds_covariates100$ID)
dfw <- as.data.frame(subset(all_birds_covariates100, state==3))
restricted <- as.data.frame(subset(all_birds_covariates100, state==1 | state==2))
#roost <- as.data.frame(subset(all_birds_covariates100, state==4))

###### Run the mixed effect model ###########
library(glmmTMB)
#head(walking)
#spearmanR <- cor(tw[,c(27, 28, 35:46)])
#spearmanR

#### ground roosting and tree roosting
dfw$phase[dfw$day < 15] <- "ground"
dfw$phase[dfw$day > 14] <- "tree"
gw <- subset(dfw, phase == "ground")
tw <- subset(dfw, phase == "tree")

# Make animal id a factor
gw$ID <- as.factor(gw$ID)

# order stratum-ID
d.map <- data.frame(NA_ID=unique(gw$step_id_),str_ID=1:length(unique(gw$step_id_)))
gw$str_ID <- d.map[match(gw$step_id_,d.map$NA_ID),"str_ID"]
gw <- gw[order(gw$str_ID),]
head(dfw)

# Scale covariates
gw$proads.s <- scale(gw$proads) 
gw$sroads.s <- scale(gw$sroads)
gw$dist_hw.s <- scale(gw$dist_hw)
gw$dist_infra.s <- scale(gw$dist_infra)
gw$dist_mixed.s <- scale(gw$dist_mix)
gw$dist_open.s <- scale(gw$dist_open)
gw$dist_pine.s <- scale(gw$dist_pine)
gw$dist_shrub.s <- scale(gw$dist_shrub)
gw$dist_water.s <- scale(gw$dist_water)
gw$ndvi.s <- scale(gw$ndvi)
#gw$cos_ta_.s <- scale(gw$ta_)
#gw$log_sl_.s <- scale(gw$sl_)
gw$residence_raster_walking.s <- scale(gw$residence_raster_walking)
gw$revisit_raster_walking.s <- scale(gw$revisit_raster_walking)
str(gw)
summary(gw)

tw$proads.s <- scale(tw$proads) 
tw$sroads.s <- scale(tw$sroads)
tw$dist_hw.s <- scale(tw$dist_hw)
tw$dist_infra.s <- scale(tw$dist_infra)
tw$dist_mixed.s <- scale(tw$dist_mix)
tw$dist_open.s <- scale(tw$dist_open)
tw$dist_pine.s <- scale(tw$dist_pine)
tw$dist_shrub.s <- scale(tw$dist_shrub)
tw$dist_water.s <- scale(tw$dist_water)
tw$ndvi.s <- scale(tw$ndvi)
#tw$cos_ta_.s <- scale(tw$ta_)
#tw$log_sl_.s <- scale(tw$sl_)
tw$residence_raster_walking.s <- scale(tw$residence_raster_walking)
tw$revisit_raster_walking.s <- scale(tw$revisit_raster_walking)
str(tw)
summary(tw)
head(gw)

gw <- gw[!is.na(gw$proads.s),]
gw <- gw[!is.na(gw$sroads.s),]
cor(gw[, c(43,44,48:58)], method = "pearson")

tw <- tw[!is.na(tw$proads.s),]
tw <- tw[!is.na(tw$sroads.s),]
cor(tw[, c(43,44,47:57)])
# Walking ground roost Null model
Null.gw.100 <- glmmTMB(case_ ~ + 1 + (1|str_ID) + (1|ID), family=poisson, data=gw ,map=list(theta=factor(c(NA,NA))),
                                                  start=list(theta=c(log(1e3),log(1e3))))       

summary(Null.gw.100)

# Walking ground roost revisit: recursion habitat model
glmm.TMB.random.walking_ground_revisit100 <- glmmTMB(case_ ~ -1 + sroads.s + sroads.s:revisit_raster_walking.s +
                                           dist_hw.s + dist_hw.s:revisit_raster_walking.s + 
                                           dist_open.s + dist_open.s:revisit_raster_walking.s +
                                           dist_mixed.s + dist_mixed.s:revisit_raster_walking.s +
                                           #dist_water.s + dist_water.s:revisit_raster_walking.s +
                                           dist_pine.s + dist_pine.s:revisit_raster_walking.s +
                                           dist_shrub.s + dist_shrub.s:revisit_raster_walking.s +
                                           ndvi.s + ndvi.s:revisit_raster_walking.s +
                                           cos_ta_ + log_sl_ + revisit_raster_walking.s +
                                            (1|str_ID) + (1|ID) + (0 + sroads.s| ID) + 
                                            (0+ dist_hw.s | ID) + (0 + dist_open.s | ID) + 
                                            (0 + dist_pine.s | ID) + (0 + dist_shrub.s | ID)+ 
                                            (0 + dist_mixed.s | ID) + (0 + ndvi.s | ID)+ 
                                            #(0 + dist_water.s | ID) + 
                                            (0 + revisit_raster_walking.s | ID),
                                             family=poisson, data=gw,
                                             map=list(theta=factor(c(NA,NA,1:8))),
                                             start=list(theta=c(log(1e3),log(1e3),0,0,0,0,0,0,0,0)))

summary(glmm.TMB.random.walking_ground_revisit100)

# Walking ground roost revisit: resource only model
glmm.TMB.random.walking_ground_revisit_resource_only100 <- glmmTMB(case_ ~ -1 + #dist_water.s + 
                                                    sroads.s + dist_hw.s + dist_open.s + dist_pine.s + 
                                                    dist_shrub.s + ndvi.s + dist_mixed.s +
                                                    cos_ta_ + log_sl_ +  
                                                    (1|str_ID) + (1|ID) + (0 + dist_mixed.s| ID) +
                                                    (0 + sroads.s| ID) + (0+ dist_hw.s | ID) + 
                                                    (0 + dist_open.s | ID) + (0 + dist_pine.s | ID) + 
                                                    (0 + dist_shrub.s | ID) + 
                                                    #(0 + dist_water.s | ID) +
                                                    (0 + ndvi.s | ID),
                                                  family=poisson, data=gw,
                                                  map=list(theta=factor(c(NA,NA,1:7))),
                                                  start=list(theta=c(log(1e3),log(1e3),0,0,0,0,0,0,0)))

summary(glmm.TMB.random.walking_ground_revisit_resource_only100)

#Walking ground roost revisit: resource only model
glmm.TMB.random.walking_ground_revisit_recursion_only100 <- glmmTMB(case_ ~ -1 + cos_ta_ + log_sl_ + 
                                                                  revisit_raster_walking.s + (1|str_ID) +
                                                                    (1|ID) + (0 + revisit_raster_walking.s | ID),
                                                                family=poisson, data=gw,
                                                                map=list(theta=factor(c(NA,NA,1))),
                                                                start=list(theta=c(log(1e3),log(1e3),0)))

summary(glmm.TMB.random.walking_ground_revisit_recursion_only100)

# Walking ground roost residence time
glmm.TMB.random.walking_ground_residence100 <- glmmTMB(case_ ~ -1 + #dist_water.s + dist_water.s:residence_raster_walking.s +
                                            sroads.s + sroads.s:residence_raster_walking.s +
                                            dist_hw.s + dist_hw.s:residence_raster_walking.s + 
                                            dist_open.s + dist_open.s:residence_raster_walking.s +
                                            dist_pine.s + dist_pine.s:residence_raster_walking.s +
                                            dist_shrub.s + dist_shrub.s:residence_raster_walking.s +
                                            dist_mixed.s + dist_mixed.s:residence_raster_walking.s +
                                            ndvi.s + ndvi.s:residence_raster_walking.s +
                                            cos_ta_ + log_sl_ + residence_raster_walking.s +
                                            (1|str_ID) + (1|ID) + 
                                            (0 + sroads.s| ID) + 
                                            (0 + dist_hw.s | ID) + 
                                            #(0 + dist_water.s | ID) +
                                            (0 + dist_mixed.s | ID) +
                                            (0 + dist_open.s | ID) + 
                                            (0 + dist_pine.s | ID) + 
                                            (0 + dist_shrub.s | ID) + 
                                            (0 + ndvi.s | ID)+ 
                                            (0 + residence_raster_walking.s | ID),
                                          family=poisson, data=gw,
                                          map=list(theta=factor(c(NA,NA,1:8))),
                                          start=list(theta=c(log(1e3),log(1e3),0,0,0,0,0,0,0,0)))

summary(glmm.TMB.random.walking_ground_residence100)

# Walking ground roost residence time: resource only
glmm.TMB.random.walking_ground_residence_resource_only100 <- glmmTMB(case_ ~ -1 + #dist_water.s + 
                                                      sroads.s + dist_hw.s +  
                                                      dist_open.s + dist_pine.s +
                                                      dist_mixed.s +
                                                      dist_shrub.s + ndvi.s + 
                                                      cos_ta_ + log_sl_ +
                                                      (1|str_ID) + (1|ID) +
                                                      #(0 + dist_water.s| ID) + 
                                                      (0 + sroads.s| ID) +
                                                      (0 + dist_mixed.s| ID) +
                                                      (0 + dist_hw.s | ID) +
                                                      (0 + dist_open.s | ID) + 
                                                      (0 + dist_pine.s | ID) + 
                                                      (0 + dist_shrub.s | ID) + 
                                                      (0 + ndvi.s | ID),
                                                    family=poisson, data=gw,
                                                    map=list(theta=factor(c(NA,NA,1:7))),
                                                    start=list(theta=c(log(1e3),log(1e3),0,0,0,0,0,0,0)))


summary(glmm.TMB.random.walking_ground_residence_resource_only100)

# Walking ground roost residence time: recursion only
glmm.TMB.random.walking_ground_residence_recursion_only100 <- glmmTMB(case_ ~ -1 + cos_ta_ + log_sl_ + 
                                                    residence_raster_walking.s +
                                                    (1|str_ID) + (1|ID) +
                                                    (0 + residence_raster_walking.s | ID),
                                                    family=poisson, data=gw,
                                                    map=list(theta=factor(c(NA,NA, 1))),
                                                    start=list(theta=c(log(1e3),log(1e3),0)))

summary(glmm.TMB.random.walking_ground_residence_recursion_only100)

# Walking tree roost revisit
# Make animal id a factor
tw$ID <- as.factor(tw$ID)

# order stratum-ID
d.map <- data.frame(NA_ID=unique(tw$step_id_),str_ID=1:length(unique(tw$step_id_)))
tw$str_ID <- d.map[match(tw$step_id_,d.map$NA_ID),"str_ID"]
tw <- tw[order(tw$str_ID),]
head(tw)
summary(tw)

#Tree roosting and walking model: global
glmm.TMB.random.walking_tree_revisit100 <- glmmTMB(case_ ~ -1 + #dist_water.s + dist_water.s:revisit_raster_walking.s +
                                                    sroads.s + sroads.s:revisit_raster_walking.s +
                                                    dist_hw.s + dist_hw.s:revisit_raster_walking.s +
                                                    dist_mixed.s + dist_mixed.s:revisit_raster_walking.s +
                                                    dist_open.s + dist_open.s:revisit_raster_walking.s +
                                                    dist_pine.s + dist_pine.s:revisit_raster_walking.s +
                                                    dist_shrub.s + dist_shrub.s:revisit_raster_walking.s +
                                                    ndvi.s + ndvi.s:revisit_raster_walking.s +
                                                    cos_ta_ + log_sl_ + revisit_raster_walking.s +
                                                    (1|str_ID) + (1|ID) +
                                                    #(0 + dist_water.s| ID) 
                                                    (0 + sroads.s| ID) + 
                                                    (0+ dist_hw.s | ID) + (0 + dist_open.s | ID) + 
                                                    (0 + dist_pine.s | ID)+ (0 + dist_shrub.s | ID) +
                                                    (0 + ndvi.s | ID)+ (0 + dist_mixed.s | ID)+
                                                    (0 + revisit_raster_walking.s | ID),
                                                  family=poisson,
                                                  data=tw,
                                                  map=list(theta=factor(c(NA,NA,1:8))),
                                                  start=list(theta=c(log(1e3),log(1e3),0,0,0,0,0,0,0,0)))

summary(glmm.TMB.random.walking_tree_revisit100)

#Tree roosting and walking model: resource only
glmm.TMB.random.walking_tree_revisit_resource_only100 <- glmmTMB(case_ ~ -1 + #dist_water.s + 
                                                  sroads.s + 
                                                  dist_hw.s +  dist_open.s +  dist_mixed.s +
                                                  dist_pine.s + dist_shrub.s +  
                                                  ndvi.s +
                                                  cos_ta_ + log_sl_ + 
                                                  #(0 + dist_water.s| ID) + 
                                                  (0 + sroads.s| ID) + 
                                                  (0+ dist_hw.s | ID) + (0 + dist_open.s | ID) + 
                                                  (0 + dist_pine.s | ID)+ (0 + dist_shrub.s | ID) + 
                                                  (0 + ndvi.s | ID)+(0 + dist_mixed.s | ID)+
                                                  (1|str_ID) + (1|ID),
                                                family=poisson,
                                                data=tw,
                                                map=list(theta=factor(c(1:7, NA, NA))),
                                                start=list(theta=c(0,0,0,0,0,0,0,log(1e3),log(1e3))))

summary(glmm.TMB.random.walking_tree_revisit_resource_only100)

#Tree roosting and walking model: recursion only
glmm.TMB.random.walking_tree_revisit_recursion_only100 <- glmmTMB(case_ ~ -1 + cos_ta_ + log_sl_ + 
                                                                revisit_raster_walking.s +
                                                                (1|str_ID) + (1|ID) +
                                                                (0 + revisit_raster_walking.s | ID),
                                                              family=poisson,
                                                              data=tw,
                                                              map=list(theta=factor(c(NA, NA, 1))),
                                                              start=list(theta=c(log(1e3),log(1e3), 0)))

summary(glmm.TMB.random.walking_tree_revisit_recursion_only100)

# Walking tree roost residence time: global
glmm.TMB.random.walking_tree_residence100 <- glmmTMB(case_ ~ -1 + #dist_water.s + dist_water.s:residence_raster_walking.s +
                                                      sroads.s + sroads.s:residence_raster_walking.s +
                                                      dist_hw.s + dist_hw.s:residence_raster_walking.s +
                                                      dist_mixed.s + dist_mixed.s:residence_raster_walking.s +
                                                      dist_open.s + dist_open.s:residence_raster_walking.s +
                                                      dist_pine.s + dist_pine.s:residence_raster_walking.s +
                                                      dist_shrub.s + dist_shrub.s:residence_raster_walking.s +
                                                      ndvi.s + ndvi.s:residence_raster_walking.s +
                                                      cos_ta_ + #cos_ta_.s:residence_raster_walking.s +
                                                      log_sl_ + #log_sl_.s:residence_raster_walking.s +
                                                      residence_raster_walking.s +
                                                      (1|str_ID) + (1|ID) +
                                                      #(0 + dist_water.s| ID) + #(0 + proads.s:residence_raster_walking.s | ID) +
                                                      (0 + sroads.s| ID) + #(0 + sroads.s:residence_raster_walking.s | ID) +
                                                      (0 + dist_hw.s | ID) + #(0 + dist_hw.s:residence_raster_walking.s | ID) +
                                                      (0 + dist_mixed.s | ID) +
                                                      (0 + dist_open.s | ID) + #(0 + dist_open.s:residence_raster_walking.s | ID) +
                                                      (0 + dist_pine.s | ID) + #(0 + dist_pine.s:residence_raster_walking.s | ID) +
                                                      (0 + dist_shrub.s | ID) + #(0 + dist_shrub.s:residence_raster_walking.s | ID)+
                                                      (0 + ndvi.s | ID)+ #(0 + ndvi.s:residence_raster_walking.s | ID)+
                                                      (0 + residence_raster_walking.s | ID),
                                                    data = tw,
                                                    family = poisson,
                                                    map=list(theta=factor(c(NA,NA,1:8))),
                                                    start=list(theta=c(log(1e3),log(1e3),0,0,0,0,0,0,0,0)))

summary(glmm.TMB.random.walking_tree_residence100)


# Walking tree roost residence time: resource only
glmm.TMB.random.walking_tree_residence_resource_only100 <- glmmTMB(case_ ~ -1 + dist_mixed.s + sroads.s + 
                                                    dist_hw.s + dist_open.s + 
                                                    dist_pine.s + dist_shrub.s + 
                                                    ndvi.s + 
                                                    cos_ta_ + log_sl_ + 
                                                    (1|str_ID) + (1|ID) +
                                                    (0 + sroads.s| ID) + 
                                                    (0 + dist_hw.s | ID) +
                                                    (0 + dist_mixed.s | ID) +
                                                    (0 + dist_open.s | ID) + 
                                                    (0 + dist_pine.s | ID) + 
                                                    (0 + dist_shrub.s | ID) +
                                                    (0 + ndvi.s | ID),
                                                  data = tw,
                                                  family = poisson,
                                                  map=list(theta=factor(c(NA,NA,1:7))),
                                                  start=list(theta=c(log(1e3),log(1e3),0,0,0,0,0,0,0)))

summary(glmm.TMB.random.walking_tree_residence_resource_only100)

# Walking tree roost residence time: recursion only
glmm.TMB.random.walking_tree_residence_recursion_only100 <- glmmTMB(case_ ~ -1 + residence_raster_walking.s +
                                                    cos_ta_ + log_sl_ +
                                                    (1|str_ID) + (1|ID) + (0 + residence_raster_walking.s | ID),
                                                  data = tw,
                                                  family = poisson,
                                                  map=list(theta=factor(c(NA,NA,1))),
                                                  start=list(theta=c(log(1e3),log(1e3),0)))

summary(glmm.TMB.random.walking_tree_residence_recursion_only100)

# Walking ground roost Null model
Null.tw.100 <- glmmTMB(case_ ~ + 1 + (1|str_ID) + (1|ID), family=poisson, data=tw ,map=list(theta=factor(c(NA,NA))),
                       start=list(theta=c(log(1e3),log(1e3))))

summary(Null.tw.100)
###### Restricted movements
#spearmanR <- cor(trestrict[,c(27, 28, 35:46)])
#spearmanR

## Restricted state

#### ground roosting and tree roosting
restricted$phase[restricted$day < 15] <- "ground"
restricted$phase[restricted$day > 14] <- "tree"
grestrict <- subset(restricted, phase == "ground")
trestrict <- subset(restricted, phase == "tree")

# Make animal id a factor
grestrict$ID <- as.factor(grestrict$ID)

# order stratum-ID
d.map <- data.frame(NA_ID=unique(grestrict$step_id_),str_ID=1:length(unique(grestrict$step_id_)))
grestrict$str_ID <- d.map[match(grestrict$step_id_,d.map$NA_ID),"str_ID"]
grestrict <- grestrict[order(grestrict$str_ID),]
head(grestrict)

# Scale covariates
grestrict$proads.s <- scale(grestrict$proads) 
grestrict$sroads.s <- scale(grestrict$sroads)
grestrict$dist_hw.s <- scale(grestrict$dist_hw)
grestrict$dist_infra.s <- scale(grestrict$dist_infra)
grestrict$dist_mixed.s <- scale(grestrict$dist_mix)
grestrict$dist_open.s <- scale(grestrict$dist_open)
grestrict$dist_pine.s <- scale(grestrict$dist_pine)
grestrict$dist_shrub.s <- scale(grestrict$dist_shrub)
grestrict$dist_water.s <- scale(grestrict$dist_water)
grestrict$ndvi.s <- scale(grestrict$ndvi)
#gw$cos_ta_.s <- scale(gw$ta_)
#gw$log_sl_.s <- scale(gw$sl_)
grestrict$residence_raster_foraging.s <- scale(grestrict$residence_raster_foraging)
grestrict$revisit_raster_foraging.s <- scale(grestrict$revisit_raster_foraging)
str(grestrict)
summary(grestrict)

trestrict$proads.s <- scale(trestrict$proads) 
trestrict$sroads.s <- scale(trestrict$sroads)
trestrict$dist_hw.s <- scale(trestrict$dist_hw)
trestrict$dist_infra.s <- scale(trestrict$dist_infra)
trestrict$dist_mixed.s <- scale(trestrict$dist_mix)
trestrict$dist_open.s <- scale(trestrict$dist_open)
trestrict$dist_pine.s <- scale(trestrict$dist_pine)
trestrict$dist_shrub.s <- scale(trestrict$dist_shrub)
trestrict$dist_water.s <- scale(trestrict$dist_water)
trestrict$ndvi.s <- scale(trestrict$ndvi)
#tw$cos_ta_.s <- scale(tw$ta_)
#tw$log_sl_.s <- scale(tw$sl_)
trestrict$residence_raster_foraging.s <- scale(trestrict$residence_raster_foraging)
trestrict$revisit_raster_foraging.s <- scale(trestrict$revisit_raster_foraging)
str(trestrict)
summary(trestrict)

grestrict <- grestrict[!is.na(grestrict$proads.s),]
grestrict <- grestrict[!is.na(grestrict$sroads.s),]
grestrict <- grestrict[!is.na(grestrict$revisit_raster_foraging.s),]
cor(grestrict[, c(43,44,47:58)], method = "pearson")

trestrict <- trestrict[!is.na(trestrict$proads.s),]
trestrict <- trestrict[!is.na(trestrict$sroads.s),]
trestrict <- trestrict[!is.na(trestrict$revisit_raster_foraging.s),]
cor(trestrict[, c(43,44,47:57)])
library(glmmTMB)
# Walking ground roost Null model
Null.grestrict.100 <- glmmTMB(case_ ~ + 1 + (1|str_ID) + (1|ID), family=poisson, data=grestrict ,map=list(theta=factor(c(NA,NA))),
                       start=list(theta=c(log(1e3),log(1e3))))
summary(Null.grestrict.100)

# Restricted ground roost revisit: Global
glmm.TMB.random.restrict_ground_revisit100 <- glmmTMB(case_ ~ -1 + dist_mixed.s + dist_mixed.s:revisit_raster_foraging.s +
                                                    sroads.s + sroads.s:revisit_raster_foraging.s +
                                                    dist_hw.s + dist_hw.s:revisit_raster_foraging.s +
                                                    dist_open.s + dist_open.s:revisit_raster_foraging.s +
                                                    dist_pine.s + dist_pine.s:revisit_raster_foraging.s +
                                                    dist_shrub.s + dist_shrub.s:revisit_raster_foraging.s +
                                                    ndvi.s + ndvi.s:revisit_raster_foraging.s +
                                                    cos_ta_ + #cos_ta_.s:revisit_raster_walking.s +
                                                    log_sl_ + #log_sl_.s:revisit_raster_walking.s +
                                                    revisit_raster_foraging.s +
                                                    (1|str_ID) + (1|ID) +
                                                    (0 + dist_mixed.s| ID) +
                                                    #(0 + proads.s:revisit_raster_walking.s | ID) +
                                                    (0 + sroads.s| ID) + #(0 + sroads.s:revisit_raster_walking.s | ID) +
                                                    (0 + dist_hw.s | ID) + #(0 + dist_hw.s:revisit_raster_walking.s | ID) +
                                                    (0 + dist_open.s | ID) + #(0 + dist_open.s:revisit_raster_walking.s | ID) +
                                                    (0 + dist_pine.s | ID) + #(0 + dist_pine.s:revisit_raster_walking.s | ID) +
                                                    (0 + dist_shrub.s | ID) + #(0 + dist_shrub.s:revisit_raster_walking.s | ID)+
                                                    (0 + ndvi.s | ID)+ #(0 + ndvi.s:revisit_raster_walking.s | ID)+
                                                    (0 + revisit_raster_foraging.s | ID),
                                                  #(0 + cos_ta_ | ID) +# (0 + cos_ta_.s:revisit_raster_walking.s | ID) +
                                                  #(0 + log_sl_ | ID), #+ (0 + log_sl_.s:revisit_raster_walking.s | ID),
                                                  family=poisson, data=grestrict,
                                                  map=list(theta=factor(c(NA,NA,1:8))),
                                                  start=list(theta=c(log(1e3),log(1e3),0,0,0,0,0,0,0,0)))

summary(glmm.TMB.random.restrict_ground_revisit100)

# Restricted ground roost revisit: resource only
glmm.TMB.random.restrict_ground_revisit_resource_only100 <- glmmTMB(case_ ~ -1 + dist_mixed.s + sroads.s + 
                                                     dist_hw.s + dist_open.s +
                                                     dist_pine.s + dist_shrub.s + 
                                                     ndvi.s + 
                                                     cos_ta_ + log_sl_ +
                                                     (1|str_ID) + (1|ID) +
                                                     (0 + dist_water.s| ID) +
                                                     (0 + sroads.s| ID) + 
                                                     (0 + dist_hw.s | ID) + 
                                                     (0 + dist_open.s | ID) + 
                                                     (0 + dist_pine.s | ID) + 
                                                     (0 + dist_shrub.s | ID) +
                                                     (0 + ndvi.s | ID),
                                                   family=poisson, data=grestrict,
                                                   map=list(theta=factor(c(NA,NA,1:7))),
                                                   start=list(theta=c(log(1e3),log(1e3),0,0,0,0,0,0,0)))

summary(glmm.TMB.random.restrict_ground_revisit_resource_only100)

# Restricted ground roost revisit: recursion only
glmm.TMB.random.restrict_ground_revisit_recursion_only100 <- glmmTMB(case_ ~ -1 + cos_ta_ + log_sl_ + 
                                                     revisit_raster_foraging.s +
                                                     (1|str_ID) + (1|ID) +
                                                     (0 + revisit_raster_foraging.s | ID),
                                                   family=poisson, data=grestrict,
                                                   map=list(theta=factor(c(NA,NA,1))),
                                                   start=list(theta=c(log(1e3),log(1e3),0)))

summary(glmm.TMB.random.restrict_ground_revisit_recursion_only100)

# Restricted ground roost residence time: Global

glmm.TMB.random.restrict_ground_residence100 <- glmmTMB(case_ ~ -1 + dist_mixed.s + dist_mixed.s:residence_raster_foraging.s +
                                                      sroads.s + sroads.s:residence_raster_foraging.s +
                                                      dist_hw.s + dist_hw.s:residence_raster_foraging.s +
                                                      dist_open.s + dist_open.s:residence_raster_foraging.s +
                                                      dist_pine.s + dist_pine.s:residence_raster_foraging.s +
                                                      dist_shrub.s + dist_shrub.s:residence_raster_foraging.s +
                                                      ndvi.s + ndvi.s:residence_raster_foraging.s +
                                                      cos_ta_ + #cos_ta_.s:residence_raster_walking.s +
                                                      log_sl_ + #log_sl_.s:residence_raster_walking.s +
                                                      residence_raster_foraging.s +
                                                      (1|str_ID) + (1|ID) +
                                                      (0 + dist_mixed.s| ID) + #(0 + proads.s:residence_raster_walking.s | ID) +
                                                      (0 + sroads.s| ID) + #(0 + sroads.s:residence_raster_walking.s | ID) +
                                                      (0+ dist_hw.s | ID) +  #(0 + dist_mixed.s:residence_raster_walking.s | ID)+
                                                      (0 + dist_open.s | ID) + #(0 + dist_open.s:residence_raster_walking.s | ID) +
                                                      (0 + dist_pine.s | ID) + #(0 + dist_pine.s:residence_raster_walking.s | ID) +
                                                      (0 + dist_shrub.s | ID) + #(0 + dist_water.s:residence_raster_walking.s | ID)+
                                                      (0 + ndvi.s | ID)+ #(0 + ndvi.s:residence_raster_walking.s | ID)+
                                                      (0 + residence_raster_foraging.s | ID),
                                                    #(0 + cos_ta_ | ID) + #(0 + cos_ta_.s:residence_raster_walking.s | ID),
                                                    #(0 + log_sl_ | ID),# (0 + log_sl_.s:residence_raster_walking.s | ID),
                                                    family=poisson, data=grestrict,
                                                    map=list(theta=factor(c(NA,NA,1:8))),
                                                    start=list(theta=c(log(1e3),log(1e3),0,0,0,0,0,0,0,0)))

summary(glmm.TMB.random.restrict_ground_residence100)

# Restricted ground roost residence time: resource only
glmm.TMB.random.restrict_ground_residence_resource_only100 <- glmmTMB(case_ ~ -1 + dist_mixed.s + 
                                                       sroads.s + dist_hw.s + 
                                                       dist_open.s +  dist_pine.s + 
                                                       dist_shrub.s +
                                                       ndvi.s + cos_ta_ + log_sl_ +
                                                       (1|str_ID) + (1|ID) +
                                                       (0 + sroads.s| ID) + 
                                                       (0 + dist_hw.s | ID) +
                                                       (0 + dist_mixed.s | ID) +
                                                       (0 + dist_open.s | ID) + 
                                                       (0 + dist_pine.s | ID) + 
                                                       (0 + dist_shrub.s | ID) + 
                                                       (0 + ndvi.s | ID),
                                                     family=poisson, data=grestrict,
                                                     map=list(theta=factor(c(NA,NA,1:7))),
                                                     start=list(theta=c(log(1e3),log(1e3),0,0,0,0,0,0,0)))

summary(glmm.TMB.random.restrict_ground_residence_resource_only100)

# Restricted ground roost residence time: recursion model
glmm.TMB.random.restrict_ground_residence_recursion_only100 <- glmmTMB(case_ ~ -1 + cos_ta_ + log_sl_ + 
                                                       residence_raster_foraging.s +
                                                       (1|str_ID) + (1|ID) +
                                                       (0 + residence_raster_foraging.s | ID),
                                                     family=poisson, data=grestrict,
                                                     map=list(theta=factor(c(NA,NA,1))),
                                                     start=list(theta=c(log(1e3),log(1e3),0)))

summary(glmm.TMB.random.restrict_ground_residence_recursion_only100)


# Walking tree roost revisit
# Make animal id a factor
trestrict$ID <- as.factor(trestrict$ID)

# order stratum-ID
d.map <- data.frame(NA_ID=unique(trestrict$step_id_),str_ID=1:length(unique(trestrict$step_id_)))
trestrict$str_ID <- d.map[match(trestrict$step_id_,d.map$NA_ID),"str_ID"]
trestrict<- trestrict[order(trestrict$str_ID),]
head(trestrict)
summary(trestrict)

# Walking ground roost Null model
Null.trestrict.100 <- glmmTMB(case_ ~ + 1 + (1|str_ID) + (1|ID), family=poisson, data=trestrict ,map=list(theta=factor(c(NA,NA))),
                              start=list(theta=c(log(1e3),log(1e3))))
summary(Null.trestrict.100)

# Restricted tree roost revisit: Global
glmm.TMB.random.restrict_tree_revisit100 <- glmmTMB(case_ ~ -1 + dist_mixed.s + dist_mixed.s:revisit_raster_foraging.s +
                                                  sroads.s + sroads.s:revisit_raster_foraging.s +
                                                  dist_hw.s + dist_hw.s:revisit_raster_foraging.s +
                                                  dist_open.s + dist_open.s:revisit_raster_foraging.s +
                                                  dist_pine.s + dist_pine.s:revisit_raster_foraging.s +
                                                  dist_shrub.s + dist_shrub.s:revisit_raster_foraging.s +
                                                  ndvi.s + ndvi.s:revisit_raster_foraging.s +
                                                  cos_ta_ + #cos_ta_.s:revisit_raster_walking.s +
                                                  log_sl_ + #log_sl_.s:revisit_raster_walking.s +
                                                  revisit_raster_foraging.s +
                                                  (0 + dist_mixed.s| ID) +
                                                  #(0 + proads.s:revisit_raster_walking.s | ID) +
                                                  (0 + sroads.s| ID) + #(0 + sroads.s:revisit_raster_walking.s | ID) +
                                                  (0+ dist_hw.s | ID) +# (0 + dist_hw.s:revisit_raster_walking.s | ID) +
                                                  # (0 + dist_mixed.s:revisit_raster_walking.s | ID) +
                                                  (0 + dist_open.s | ID) + 
                                                  #(0 + dist_open.s:revisit_raster_walking.s | ID) +
                                                  (0 + dist_pine.s | ID)+ #(0 + dist_pine.s:revisit_raster_walking.s | ID) +
                                                  (0 + dist_shrub.s | ID)+# (0 + dist_water.s:revisit_raster_walking.s | ID)+
                                                  (0 + ndvi.s | ID)+ #(0 + ndvi.s:revisit_raster_walking.s | ID)+
                                                  (0 + revisit_raster_foraging.s | ID)+
                                                  (1|str_ID) + (1|ID),
                                                family=poisson,
                                                data=trestrict,
                                                map=list(theta=factor(c(1:8, NA,NA))),
                                                start=list(theta=c(0,0,0,0,0,0,0,0,log(1e3),log(1e3))))

summary(glmm.TMB.random.restrict_tree_revisit100)

# Restricted tree roost revisit: Resource only
glmm.TMB.random.restrict_tree_revisit_resource_only100 <- glmmTMB(case_ ~ -1 + dist_mixed.s + sroads.s +
                                                   dist_hw.s + dist_open.s + 
                                                   dist_pine.s + dist_shrub.s + 
                                                   ndvi.s +
                                                   cos_ta_ + log_sl_ +
                                                   (0 + dist_mixed.s| ID) +
                                                   (0 + sroads.s| ID) +
                                                   (0 + dist_hw.s | ID) +
                                                   (0 + dist_open.s | ID) + 
                                                   (0 + dist_pine.s | ID)+ 
                                                   (0 + dist_shrub.s | ID)+
                                                   (0 + ndvi.s | ID)+ 
                                                   (1|str_ID) + (1|ID),
                                                 family=poisson,
                                                 data=trestrict,
                                                 map=list(theta=factor(c(1:7, NA,NA))),
                                                 start=list(theta=c(0,0,0,0,0,0,0,log(1e3),log(1e3))))

summary(glmm.TMB.random.restrict_tree_revisit_resource_only100)

# Restricted tree roost revisit: Recursion only
glmm.TMB.random.restrict_tree_revisit_recursion_only100 <- glmmTMB(case_ ~ -1 + cos_ta_ + 
                                                   log_sl_ + revisit_raster_foraging.s +
                                                   (0 + revisit_raster_foraging.s | ID)+
                                                   (1|str_ID) + (1|ID),
                                                 family=poisson,
                                                 data=trestrict,
                                                 map=list(theta=factor(c(1, NA, NA))),
                                                 start=list(theta=c(0,log(1e3),log(1e3))))

summary(glmm.TMB.random.restrict_tree_revisit_recursion_only100)

# Restricted tree roost residence time: Global
glmm.TMB.random.restrict_tree_residence100 <- glmmTMB(case_ ~ -1 + dist_mixed.s + dist_mixed.s:residence_raster_foraging.s +
                                                    sroads.s + sroads.s:residence_raster_foraging.s +
                                                    dist_hw.s + dist_hw.s:residence_raster_foraging.s +
                                                    dist_open.s + dist_open.s:residence_raster_foraging.s +
                                                    dist_pine.s + dist_pine.s:residence_raster_foraging.s +
                                                    dist_shrub.s + dist_shrub.s:residence_raster_foraging.s +
                                                    ndvi.s + ndvi.s:residence_raster_foraging.s +
                                                    cos_ta_ + #cos_ta_.s:residence_raster_walking.s +
                                                    log_sl_ + #log_sl_.s:residence_raster_walking.s +
                                                    residence_raster_foraging.s +
                                                    (1|str_ID) + (1|ID) +
                                                    (0 + dist_mixed.s| ID) + #(0 + proads.s:residence_raster_walking.s | ID) +
                                                    (0 + sroads.s| ID) + #(0 + sroads.s:residence_raster_walking.s | ID) +
                                                    (0+ dist_hw.s | ID) + #(0 + dist_mixed.s:residence_raster_walking.s | ID)+
                                                    (0 + dist_open.s | ID) + #(0 + dist_open.s:residence_raster_walking.s | ID) +
                                                    (0 + dist_pine.s | ID) + #(0 + dist_pine.s:residence_raster_walking.s | ID) +
                                                    (0 + dist_shrub.s | ID) +#(0 + dist_water.s:residence_raster_walking.s | ID)+
                                                    (0 + ndvi.s | ID)+ #(0 + ndvi.s:residence_raster_walking.s | ID)+
                                                    (0 + residence_raster_foraging.s | ID),
                                                  data = trestrict,
                                                  family = poisson,
                                                  map=list(theta=factor(c(NA,NA,1:8))),
                                                  start=list(theta=c(log(1e3),log(1e3),0,0,0,0,0,0,0,0)))

summary(glmm.TMB.random.restrict_tree_residence100)

# Restricted tree roost residence time: Resource only
glmm.TMB.random.restrict_tree_residence_resource_only100 <- glmmTMB(case_ ~ -1 + dist_mixed.s +
                                                     sroads.s + dist_hw.s + 
                                                     dist_open.s + dist_pine.s + 
                                                     dist_shrub.s + 
                                                     ndvi.s + cos_ta_ + 
                                                     log_sl_ + 
                                                     (1|str_ID) + (1|ID) +
                                                     (0 + dist_mixed.s| ID) + 
                                                     (0 + sroads.s| ID) + 
                                                     (0+ dist_hw.s | ID) + 
                                                     (0 + dist_open.s | ID) + 
                                                     (0 + dist_pine.s | ID) + 
                                                     (0 + dist_shrub.s | ID) +
                                                     (0 + ndvi.s | ID),
                                                   data = trestrict,
                                                   family = poisson,
                                                   map=list(theta=factor(c(NA,NA,1:7))),
                                                   start=list(theta=c(log(1e3),log(1e3),0,0,0,0,0,0,0)))

summary(glmm.TMB.random.restrict_tree_residence_resource_only100)

# Restricted tree roost residence time: Recursion only
glmm.TMB.random.restrict_tree_residence_recursion_only100 <- glmmTMB(case_ ~ -1 + cos_ta_ + log_sl_ + 
                                                     residence_raster_foraging.s +
                                                     (1|str_ID) + (1|ID) +
                                                     (0 + residence_raster_foraging.s | ID),
                                                   data = trestrict,
                                                   family = poisson,
                                                   map=list(theta=factor(c(NA,NA,1))),
                                                   start=list(theta=c(log(1e3),log(1e3),0)))

summary(glmm.TMB.random.restrict_tree_residence_recursion_only100)

# Model selection
library(AICcmodavg)
library(glmmTMB)
#### models for restricted ground revisit ####
Cand.modelsRGRev <- list(glmm.TMB.random.restrict_ground_revisit100,
                    glmm.TMB.random.restrict_ground_revisit_recursion_only100,
                    glmm.TMB.random.restrict_ground_revisit_resource_only100)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.modelsRGRev), sep = " ")

##generate AICc table
aictab(cand.set = Cand.modelsRGRev, modnames = Modnames, sort = TRUE)

##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.modelsRGRev, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

#### models for restricted ground residence ####
Cand.modelsRGRes <- list(glmm.TMB.random.restrict_ground_residence100, 
                         glmm.TMB.random.restrict_ground_residence_recursion_only100, 
                         glmm.TMB.random.restrict_ground_residence_resource_only100)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.modelsRGRes), sep = " ")

##generate AICc table
aictab(cand.set = Cand.modelsRGRes, modnames = Modnames, sort = TRUE)

##round to 4 digits after decimal point and give log-likelihood
 print(aictab(cand.set = Cand.modelsRGRes, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)
 
 #### models for exploratory ground revisit ####
 Cand.modelsWGRev <- list(glmm.TMB.random.walking_ground_revisit100,
                          glmm.TMB.random.walking_ground_revisit_recursion_only100,
                          glmm.TMB.random.walking_ground_revisit_resource_only100)
 
 ##create a vector of names to trace back models in set
 Modnames <- paste("mod", 1:length(Cand.modelsWGRev), sep = " ")
 
 ##generate AICc table
 aictab(cand.set = Cand.modelsWGRev, modnames = Modnames, sort = TRUE)
 
 ##round to 4 digits after decimal point and give log-likelihood
 print(aictab(cand.set = Cand.modelsWGRev, modnames = Modnames, sort = TRUE),
       digits = 4, LL = TRUE)

 #### models for exploratory ground residence ####
 Cand.modelsWGRes <- list(glmm.TMB.random.walking_ground_residence100, 
                          glmm.TMB.random.walking_ground_residence_recursion_only100, 
                          glmm.TMB.random.walking_ground_residence_resource_only100)
 
 ##create a vector of names to trace back models in set
 Modnames <- paste("mod", 1:length(Cand.modelsWGRes), sep = " ")
 
 ##generate AICc table
 aictab(cand.set = Cand.modelsWGRes, modnames = Modnames, sort = TRUE)
 
 ##round to 4 digits after decimal point and give log-likelihood
 print(aictab(cand.set = Cand.modelsWGRes, modnames = Modnames, sort = TRUE),
       digits = 4, LL = TRUE)
 
 #### models for restricted tree revisit ####
 Cand.modelsRTRev <- list(glmm.TMB.random.restrict_tree_revisit100,
                          glmm.TMB.random.restrict_tree_revisit_recursion_only100,
                          glmm.TMB.random.restrict_tree_revisit_resource_only100)
 
 ##create a vector of names to trace back models in set
 Modnames <- paste("mod", 1:length(Cand.modelsRTRev), sep = " ")
 
 ##generate AICc table
 aictab(cand.set = Cand.modelsRTRev, modnames = Modnames, sort = TRUE)
 
 ##round to 4 digits after decimal point and give log-likelihood
 print(aictab(cand.set = Cand.modelsRTRev, modnames = Modnames, sort = TRUE),
       digits = 4, LL = TRUE)
 
 #### models for restricted Tree residence ####
 Cand.modelsRTRes <- list(glmm.TMB.random.restrict_tree_residence100, 
                          glmm.TMB.random.restrict_tree_residence_recursion_only100, 
                          glmm.TMB.random.restrict_tree_residence_resource_only100)
 
 ##create a vector of names to trace back models in set
 Modnames <- paste("mod", 1:length(Cand.modelsRTRes), sep = " ")
 
 ##generate AICc table
 aictab(cand.set = Cand.modelsRTRes, modnames = Modnames, sort = TRUE)
 
 ##round to 4 digits after decimal point and give log-likelihood
 print(aictab(cand.set = Cand.modelsRTRes, modnames = Modnames, sort = TRUE),
       digits = 4, LL = TRUE)
 
 #### models for walking tree revisit ####
 Cand.modelsWTRev <- list(glmm.TMB.random.walking_tree_revisit100,
                          glmm.TMB.random.walking_tree_revisit_recursion_only100,
                          glmm.TMB.random.walking_tree_revisit_resource_only100)
 
 ##create a vector of names to trace back models in set
 Modnames <- paste("mod", 1:length(Cand.modelsWTRev), sep = " ")
 
 ##generate AICc table
 aictab(cand.set = Cand.modelsWTRev, modnames = Modnames, sort = TRUE)
 
 ##round to 4 digits after decimal point and give log-likelihood
 print(aictab(cand.set = Cand.modelsWTRev, modnames = Modnames, sort = TRUE),
       digits = 4, LL = TRUE)
 
 #### models for walking Tree residence ####
 Cand.modelsWTRes <- list(glmm.TMB.random.walking_tree_residence100, 
                          glmm.TMB.random.walking_tree_residence_recursion_only100, 
                          glmm.TMB.random.walking_tree_residence_resource_only100)
 
 ##create a vector of names to trace back models in set
 Modnames <- paste("mod", 1:length(Cand.modelsWTRes), sep = " ")
 
 ##generate AICc table
 aictab(cand.set = Cand.modelsWTRes, modnames = Modnames, sort = TRUE)
 
 ##round to 4 digits after decimal point and give log-likelihood
 print(aictab(cand.set = Cand.modelsWTRes, modnames = Modnames, sort = TRUE),
       digits = 4, LL = TRUE)
 
 # model performance
 library(gmodels)
 library(pROC)

 #Calc AUC
 # Residence Ground
 pred_resource_gr <- predict(glmm.TMB.random.restrict_ground_residence_resource_only100, grestrict, type="response",allow.new.levels=T)
 
 AUC_resource_gr <- ci(roc(grestrict$case_, pred_resource_gr)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_resource_gr
 
 pred_global_gr <- predict(glmm.TMB.random.restrict_ground_residence100, grestrict, type="response",allow.new.levels=T)
 
 AUC_global_gr <- ci(roc(grestrict$case_, pred_global_gr)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_global_gr
 
 pred_residence_gr <- predict(glmm.TMB.random.restrict_ground_residence_recursion_only100, grestrict, type="response",allow.new.levels=T)
 
 AUC_residence_gr <- ci(roc(grestrict$case_, pred_residence_gr)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_residence_gr
 
 # revisit ground roosting restrict
 pred_resource_gr_rev <- predict(glmm.TMB.random.restrict_ground_revisit_resource_only100, grestrict, type="response",allow.new.levels=T)
 
 AUC_resource_gr_rev <- ci(roc(grestrict$case_, pred_resource_gr_rev)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_resource_gr_rev
 
 pred_global_gr_rev <- predict(glmm.TMB.random.restrict_ground_revisit100, grestrict, type="response",allow.new.levels=T)
 
 AUC_global_gr_rev <- ci(roc(grestrict$case_, pred_global_gr_rev)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_global_gr_rev
 
 pred_residence_gr_rev <- predict(glmm.TMB.random.restrict_ground_revisit_recursion_only100, grestrict, type="response",allow.new.levels=T)
 
 AUC_revisit_gr_rev <- ci(roc(grestrict$case_, pred_residence_gr_rev)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_revisit_gr_rev
 
 # revisit tree roosting restrict
 pred_resource_tr_rev <- predict(glmm.TMB.random.restrict_tree_revisit_resource_only100, trestrict, type="response",allow.new.levels=T)
 
 AUC_resource_tr_rev <- ci(roc(trestrict$case_, pred_resource_tr_rev)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_resource_tr_rev
 
 pred_global_tr_rev <- predict(glmm.TMB.random.restrict_tree_revisit100, trestrict, type="response",allow.new.levels=T)
 
 AUC_global_tr_rev <- ci(roc(trestrict$case_, pred_global_tr_rev)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_global_tr_rev
 
 pred_revisit_tr_rev <- predict(glmm.TMB.random.restrict_tree_revisit_recursion_only100, trestrict, type="response",allow.new.levels=T)
 
 AUC_revisit_tr_rev <- ci(roc(trestrict$case_, pred_residence_tr_rev)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_revisit_tr_rev
 
 # revisit tree roosting residence
 pred_resource_tr <- predict(glmm.TMB.random.restrict_tree_residence_resource_only100, trestrict, type="response",allow.new.levels=T)
 
 AUC_resource_tr <- ci(roc(trestrict$case_, pred_resource_tr)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_resource_tr_rev
 
 pred_global_tr <- predict(glmm.TMB.random.restrict_tree_residence100, trestrict, type="response",allow.new.levels=T)
 
 AUC_global_tr <- ci(roc(trestrict$case_, pred_global_tr)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_global_tr
 
 pred_residence_tr <- predict(glmm.TMB.random.restrict_tree_residence_recursion_only100, trestrict, type="response",allow.new.levels=T)
 
 AUC_residence_tr <- ci(roc(trestrict$case_, pred_residence_tr)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_residence_tr
 #################################### Walking #####################################
 # Residence Ground mobile
 pred_resource_gw <- predict(glmm.TMB.random.walking_ground_residence_resource_only100, gw, type="response",allow.new.levels=T)
 
 AUC_resource_gw <- ci(roc(gw$case_, pred_resource_gw)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_resource_gw
 
 pred_global_gw <- predict(glmm.TMB.random.walking_ground_residence100, gw, type="response",allow.new.levels=T)
 
 AUC_global_gw <- ci(roc(gw$case_, pred_global_gw)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_global_gw
 
 pred_residence_gw <- predict(glmm.TMB.random.walking_ground_residence_recursion_only100, gw, type="response",allow.new.levels=T)
 
 AUC_residence_gw <- ci(roc(gw$case_, pred_residence_gw)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_residence_gw
 
 # revisit ground roosting restrict
 pred_resource_gw_rev <- predict(glmm.TMB.random.walking_ground_revisit_resource_only100, gw, type="response",allow.new.levels=T)
 
 AUC_resource_gw_rev <- ci(roc(gw$case_, pred_resource_gw_rev)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_resource_gw_rev
 
 pred_global_gw_rev <- predict(glmm.TMB.random.walking_ground_revisit100, gw, type="response",allow.new.levels=T)
 
 AUC_global_gw_rev <- ci(roc(gw$case_, pred_global_gw_rev)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_global_gw_rev
 
 pred_residence_gw_rev <- predict(glmm.TMB.random.walking_ground_revisit_recursion_only100, gw, type="response",allow.new.levels=T)
 
 AUC_revisit_gw_rev <- ci(roc(gw$case_, pred_residence_gw_rev)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_revisit_gw_rev
 
 # revisit tree roosting restrict
 pred_resource_tw_rev <- predict(glmm.TMB.random.walking_tree_revisit_resource_only100, tw, type="response",allow.new.levels=T)
 
 AUC_resource_tw_rev <- ci(roc(tw$case_, pred_resource_tw_rev)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_resource_tw_rev
 
 pred_global_tw_rev <- predict(glmm.TMB.random.walking_tree_revisit100, tw, type="response",allow.new.levels=T)
 
 AUC_global_tw_rev <- ci(roc(tw$case_, pred_global_tw_rev)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_global_tw_rev
 
 pred_revisit_tw_rev <- predict(glmm.TMB.random.walking_tree_revisit_recursion_only100, tw, type="response",allow.new.levels=T)
 
 AUC_revisit_tw_rev <- ci(roc(tw$case_, pred_revisit_tw_rev)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_revisit_tw_rev
 
 # revisit tree roosting residence
 pred_resource_tw <- predict(glmm.TMB.random.walking_tree_residence_resource_only100, tw, type="response",allow.new.levels=T)
 
 AUC_resource_tw <- ci(roc(tw$case_, pred_resource_tw)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_resource_tw_rev
 
 pred_global_tw <- predict(glmm.TMB.random.walking_tree_residence100, tw, type="response",allow.new.levels=T)
 
 AUC_global_tw <- ci(roc(tw$case_, pred_global_tw)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_global_tw
 
 pred_residence_tw <- predict(glmm.TMB.random.walking_tree_residence_recursion_only100, tw, type="response",allow.new.levels=T)
 
 AUC_residence_tw <- ci(roc(tw$case_, pred_residence_tw)) %>% as.numeric %>%
   round(digits = 2)
 
 AUC_residence_tw
 
library(performance)
model_performance(glmm.TMB.random.walking_tree_residence_recursion_only100)
# Coefficient Plots
# remove objects from environment
rm(trestrict)
library(broom)
library(dotwhisker)
library(ggpubr)
library(ggplot2)
library(glmmTMB)
library(tidyverse)
library(scales)
## basic coefficient plot
residence_restricted <- dotwhisker::dwplot(list("Ground roosting" = glmm.TMB.random.restrict_ground_residence100,
            "Tree roosting" = glmm.TMB.random.restrict_tree_residence100), effects="fixed",
            dot_args = list(size = 3),
            whisker_args = list(size = 0.75)) %>% 
            relabel_predictors(
              c(sroads.s = "Secondary roads",
                dist_hw.s = "Hardwoods",
                dist_mixed.s = "Mixed pine-hardwoods",
                dist_open.s = "Open treeless",
                dist_pine.s = "Pine",
                dist_shrub.s = "Shrub/scrub",
                ndvi.s = "Normalized density vegetation index",
                cos_ta_ = "Cosine turn angle",
                log_sl_ = "Log of step length",
                residence_raster_foraging.s = "Residence time",
                "sroads.s:residence_raster_foraging.s" = "Secondary roads:Residence time",
                "dist_hw.s:residence_raster_foraging.s" = "Hardwoods:Residence time",
                "dist_mixed.s:residence_raster_foraging.s" = "Mixed pine hardwoods:Residence time",
                "dist_open.s:residence_raster_foraging.s" = "Open treeless:Residence time",
                "dist_pine.s:residence_raster_foraging.s" = "Pine:Residence time",
                "dist_shrub.s:residence_raster_foraging.s" = "Shrub/scrub:Residence time",
                "ndvi.s:residence_raster_foraging.s" = "Normalized density vegetation index:Residence time"
              )) +
            geom_vline(xintercept=0,lty=2) +           
            scale_x_continuous(limits = c(-2.0,2.0)) + xlab("Coeffcient estimate") +
            guides(color = guide_legend(title = "Restricted state"), size = 10) +
            theme_bw() 

residence_restricted <- residence_restricted + theme(text = element_text(size = 10))
residence_restricted


revisit_restricted <- dotwhisker::dwplot(list("Ground roosting revisit restricted" = glmm.TMB.random.restrict_ground_revisit100,
    "Tree roosting revisit restricted" = glmm.TMB.random.restrict_tree_revisit100), effects="fixed",
    dot_args = list(size = 3),
    whisker_args = list(size = 0.75)) %>% 
    relabel_predictors(
    c(sroads.s = "Secondary roads",
      dist_hw.s = "Hardwoods",
      dist_mixed.s = "Mixed pine-hardwoods",
      dist_open.s = "Open treeless",
      dist_pine.s = "Pine",
      dist_shrub.s = "Shrub/scrub",
      ndvi.s = "Normalized density vegetation index",
      cos_ta_ = "Cosine turn angle",
      log_sl_ = "Log of step length",
      revisit_raster_foraging.s = "Revisit",
      "sroads.s:revisit_raster_foraging.s" = "Secondary roads:Revisit",
      "dist_hw.s:revisit_raster_foraging.s" = "Hardwoods:Revisit",
      "dist_mixed.s:revisit_raster_foraging.s" = "Mixed pine-hardwoods:Revisit",
      "dist_open.s:revisit_raster_foraging.s" = "Open treeless:Revisit",
      "dist_pine.s:revisit_raster_foraging.s" = "Pine:Revisit",
      "dist_shrub.s:revisit_raster_foraging.s" = "Shrub/scrub:Revisit",
      "ndvi.s:revisit_raster_foraging.s" = "Normalized density vegetation index:Revisit")) +
  geom_vline(xintercept=0,lty=2) + ylab("Covariate") +
  scale_x_continuous(limits = c(-2.0,2.0)) + xlab("Coeffcient estimate") +
  guides(color = guide_legend(title = "Behavioral state",labels = label_wrap(4)), size = 10) +
  theme_bw()

revisit_restricted <- revisit_restricted + theme(text = element_text(size = 10),legend.position = "none")


revisit_restricted

restricted_figures <- ggarrange(revisit_restricted, residence_restricted, ncol = 2, widths = c(1.45, 1.75))
restricted_figures

#wgres <- dotwhisker::dwplot(glmm.TMB.random.walking_ground_residence100, effects="fixed") + geom_vline(xintercept = 0, lty = 2)
#wgrev <- dotwhisker::dwplot(glmm.TMB.random.walking_ground_revisit100, effects="fixed") + geom_vline(xintercept = 0, lty = 2)
#wtres <- dotwhisker::dwplot(glmm.TMB.random.walking_tree_residence100, effects="fixed") + geom_vline(xintercept = 0, lty = 2)
#wtrev <- dotwhisker::dwplot(glmm.TMB.random.walking_tree_revisit100, effects="fixed") + geom_vline(xintercept = 0, lty = 2)

residence_exploratory <- dotwhisker::dwplot(list("Ground roosting" = glmm.TMB.random.walking_ground_residence100,
  "Tree roosting" = glmm.TMB.random.walking_tree_residence100), effects="fixed",
  dot_args = list(size = 3),
  whisker_args = list(size = 0.75)) %>% 
  relabel_predictors(
    c(sroads.s = "Secondary roads",
      dist_hw.s = "Hardwoods",
      dist_mixed.s = "Mixed pine-hardwoods",
      dist_open.s = "Open treeless",
      dist_pine.s = "Pine",
      dist_shrub.s = "Shrub/scrub",
      ndvi.s = "Normalized density vegetation index",
      cos_ta_ = "Cosine turn angle",
      log_sl_ = "Log of step length",
      residence_raster_walking.s = "Residence time",
      "sroads.s:residence_raster_walking.s" = "Secondary roads:Residence time",
      "dist_hw.s:residence_raster_walking.s" = "Hardwoods:Residence time",
      "dist_mixed.s:residence_raster_walking.s" = "Mixed pine-hardwoods:Residence time",
      "dist_open.s:residence_raster_walking.s" = "Open treeless:Residence time",
      "dist_pine.s:residence_raster_walking.s" = "Pine:Residence time",
      "dist_shrub.s:residence_raster_walking.s" = "Shrub/scrub:Residence time",
      "ndvi.s:residence_raster_walking.s" = "Normalized density vegetation index:Residence time"
    )) +
  geom_vline(xintercept=0,lty=2) +           
  scale_x_continuous(limits = c(-3.00,3.00)) + xlab("Coeffcient estimate") +
  guides(color = guide_legend(title = "Mobile state")) +
  theme_bw()

residence_exploratory <- residence_exploratory + theme(text = element_text(size = 10))


revisit_exploratory <- dotwhisker::dwplot(list("ground roosting" = glmm.TMB.random.walking_ground_revisit100,
  "tree roosting" = glmm.TMB.random.walking_tree_revisit100), effects="fixed",
  dot_args = list(size = 3),
  whisker_args = list(size = 0.75)) %>% 
  relabel_predictors(
    c(sroads.s = "Secondary roads",
      dist_hw.s = "Hardwoods",
      dist_mixed.s = "Mixed pine-hardwoods",
      dist_open.s = "Open treeless",
      dist_pine.s = "Pine",
      dist_shrub.s = "Shrub/scrub",
      ndvi.s = "Normalized density vegetation index",
      cos_ta_ = "Cosine turn angle",
      log_sl_ = "Log of step length",
      revisit_raster_walking.s = "Revisit",
      "sroads.s:revisit_raster_walking.s" = "Secondary roads:Revisit",
      "dist_hw.s:revisit_raster_walking.s" = "Hardwoods:Revisit",
      "dist_mixed.s:revisit_raster_walking.s" = "Mixed pine-hardwoods:Revisit",
      "dist_open.s:revisit_raster_walking.s" = "Open treeless:Revisit",
      "dist_pine.s:revisit_raster_walking.s" = "Pine:Revisit",
      "dist_shrub.s:revisit_raster_walking.s" = "Shrub/scrub:Revisit",
      "ndvi.s:revisit_raster_walking.s" = "Normalized density vegetation index:Revisit")) +
  geom_vline(xintercept=0,lty=2) + ylab("Covariate") +
  scale_x_continuous(limits = c(-3.00,3.00)) + xlab("Coeffcient estimate") +
  theme_bw()

revisit_exploratory <- revisit_exploratory + theme(text = element_text(size = 10),legend.position = "none")

Exploratory_figures <- ggarrange(revisit_exploratory, residence_exploratory, ncol = 2, widths = c(1.45, 1.75))
Exploratory_figures
