library(MODIStsp)
library(rgdal)
library(raster)
library(sp)
library(sf)
library(tidyverse)

setwd("~/GitHub/desert-fires-impact-biodiversity")

##DOWNLOADING NDVI##

# download MODIS monthly NDVI data for Southern California from 2000-01-01 to 2021-12-31
# 1 km resolution MOD13A3 (Vegetation Indices Monthly L3 Global 1km)
# also dl MODIS NDII7 (normalized burn ratio) 
MODIStsp()

##PROCESSING NDVI RASTERS##

# Crop NDVI
rasterfiles <- list.files("D:/MODIS_data/VI_Monthly_1km_v6/NDVI/", pattern = "*.tif", full.name = TRUE)
reps <- (1:length(rasterfiles))
dir.create("raster/MODIS_ndvi")
dir.create("raster/MODIS_ndvi/cropped/")

for (rep in reps) {
  
  #make new name for each file
  newfileName <- paste(sub("D:/MODIS_data/VI_Monthly_1km_v6/NDVI/MOD13A3_", "", rasterfiles[rep]), sep = "")
  
  year <- sub("NDVI_", "", newfileName) 
  year <- strsplit(year, "_")[[1]][1]

  
  #read raster layer
  raster <- raster(rasterfiles[rep])
    
  
  #create crop extent from deserts and reproject to raster crs
  crop_extent <- readOGR("shapefiles/deserts/deserts.shp") %>% 
    spTransform(crs(raster))
  
  
  #crop
  raster_crop <- crop(raster, extent(crop_extent)) #manually change to (-121.4500, -114.0000, 32.60000, 37.60000) because it too short on one side
  
  #create new folder
  dir.create(paste("raster/MODIS_ndvi/cropped/", year, sep = ""), showWarnings = FALSE)
  
  
  #save new cropped raster
  writeRaster(raster_crop, filename = paste("raster/MODIS_ndvi/cropped", year, newfileName, sep = "/"))
  
}
#confirming extent; extent of raster should be within extent of crop_extent 
  #(i.e. xmin_raster < xmin_crop_extent, xmax_raster > xmin_crop_extent, ymin_raster < ymin_crop_extent, ymax_raster > ymin_crop_extent

test <- raster(raster_crop)
extent(test)
extent(crop_extent)



#check raster files
env_data<- list.files(paste("raster/MODIS_ndvi/cropped/", reps, sep = ""), pattern = "*.tif", full.name = TRUE)

r.env <- list()
for(i in env_data) {
  r <- try( raster(i) )
  r.env[i] <- data.frame(name = i, nrow=nrow(r), ncol=ncol(r), res=res(r)[1],               
                         proj=proj4string(r), xmin=extent(r)[1],             
                         xmax=extent(r)[2], ymin=extent(r)[3],                     
                         ymax=extent(r)[4])
}

do.call("rbind", r.env) #NULL means no errors



# average monthly NDVI values to yearly averages
reps <- (2000:2021)
dir.create("raster/MODIS_ndvi/ndvi_yrly_avg/")
for (rep in reps) {
  croppedRasters <- list.files(paste("raster/MODIS_ndvi/cropped/", rep, sep = ""), pattern = "*.tif", full.name = TRUE)
  
  #new file name
  newtotalName <- paste("ndvi", rep, sep = "_")
  
  #create multi-layer raster of all months in one year
  rasterstack <- stack(croppedRasters)
  
  #get average ndvi for the year
  ndvi_avg <- calc(rasterstack, fun = mean, na.rm = TRUE)/10000 #divide by 10000 to convert values to decimal
  
  writeRaster(ndvi_avg, filename = paste("raster/MODIS_ndvi/ndvi_yrly_avg/", newtotalName, ".tif", sep = "" ))
  
}


###EXTRACT NDVI VALUES TO Polygons###

##Extract never-burned NDVI

data <- st_read("shapefiles/deserts/minus-fires/deserts-minus-fires.shp") %>% #polygons for never-burned desert sites
        st_transform(crs = 4326) #transform data to wgs84 to match raster
rasterlist <- list.files("raster/MODIS_ndvi/ndvi_yrly_avg/", pattern = "*.tif", full.name = TRUE)
data_lst <- list()
df <- data.frame() %>% 
  mutate(year= as.numeric(), desert = as.character(), ndvi = as.numeric())

#subset data by desert and store in list
for (i in 1:length(data)){
  data_lst[[i]] <- data[i,]
}


for (i in 1:length(rasterlist)) {
  
  #select year
  year <- sub("raster/MODIS_ndvi/ndvi_yrly_avg/ndvi_", "", rasterlist[i]) 
  year <- as.numeric(strsplit(year, ".tif")[[1]][1])
  
  
  #load raster
  raster <- raster(rasterlist[i])
  
  #create data frame (year, desert, ndvi) by extracting raster values using FUNCTION = mean to desert polygon
  for (j in 1:length(data_lst)) {
   desert <- data_lst[[j]]$Desert
   ndvi <- raster::extract(raster, data_lst[[j]], fun = mean, na.rm = TRUE) 
   df[nrow(df) + 1,] = c(year, desert, ndvi[1,1])
   
   }

}

df$year <- as.numeric(df$year)
df$ndvi <- as.numeric(df$ndvi)

#save :]

write_csv(df, "data/mean_ndvi.csv")


###Extract burned NDVI
#obtain fireids from observations in order to filter fire shapefiles
obs_data <- read_csv("data/animals_2022_v4.csv") 
fireid <- unique(obs_data$fireID) 
fireid <- fireid[2:67]

data <- st_read("shapefiles/fires/fires2000.shp") %>% #polygons for never-burned desert sites
        st_transform(crs = 4326) %>% #transform data to wgs84 to match raster
        filter(fireID %in% fireid)

rasterlist <- list.files("raster/MODIS_ndvi/ndvi_yrly_avg/", pattern = "*.tif", full.name = TRUE)
data_lst <- list()
df <- data.frame() %>% 
  mutate(year= as.numeric(), desert = as.character(), ndvi = as.numeric())

#subset data by desert and store in list
data_lst[[1]] <- data %>% filter(desert == "San Joaquin")
data_lst[[2]] <- data %>% filter(desert == "Mojave")
data_lst[[3]] <- data %>% filter(desert == "Sonoran")


for (i in 1:length(rasterlist)) {
  
  #select year
  year <- sub("raster/MODIS_ndvi/ndvi_yrly_avg/ndvi_", "", rasterlist[i]) 
  year <- as.numeric(strsplit(year, ".tif")[[1]][1])
  
  
  #load raster
  raster <- raster(rasterlist[i])
  
  #create data frame (year, desert, ndvi) by extracting raster values using FUNCTION = mean to desert polygon
  for (j in 1:length(data_lst)) {
    desert <- data_lst[[j]]$desert
    ndvi <- raster::extract(raster, data_lst[[j]], fun = mean, na.rm = TRUE) 
    df[nrow(df) + 1,] = c(year, desert[1], mean(ndvi))
    
  }
  
}

df$year <- as.numeric(df$year)
df$ndvi <- as.numeric(df$ndvi)

#save :]

write_csv(df, "data/mean_ndvi_burned.csv")





###EXTRACT NDVI VALUES TO POINTS###

data <- st_read("shapefiles/animals_v3/animals_v3.shp") 
rasterlist <- list.files("raster/MODIS_ndvi/ndvi_yrly_avg/", pattern = "*.tif", full.name = TRUE)
head(data)
lst <- list()


for (i in 1:length(rasterlist)) {
  
  #select year
  year <- sub("raster/MODIS_ndvi/ndvi_yrly_avg/ndvi_", "", rasterlist[i]) 
  year <- as.numeric(strsplit(year, ".tif")[[1]][1])
  
  #load raster
  raster <- raster(rasterlist[i])
  
  #transform & subset data. Add subset to a list
  lst[[i]] <- data[data$year==year,]
  
  #extract raster values to points
  lst[[i]]$ndvi <- raster::extract(raster, lst[[i]])
  
}

#add data that does not have ndvi values
lst[[23]] <- data[data$year < 2000,]
lst[[23]]$ndvi <- NA

#concatenate list elements into one shapefile
data_ndvi <- do.call(rbind, lst)

#save :]
dir.create("shapefiles/animals_v4")
st_write(data_ndvi, "shapefiles/animals_v4/animals_v4.shp")
