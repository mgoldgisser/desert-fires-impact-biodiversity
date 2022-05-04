#Climate Data from WorldClim.org

library(raster)
library(sp)
library(rgdal)
library(sf)
library(tidyverse)

setwd("~/GitHub/desert-fires-impact-biodiversity")


#download files
download.file(url = "https://biogeo.ucdavis.edu/data/worldclim/v2.1/hist/wc2.1_2.5m_prec_1990-1999.zip", 
              destfile = "~/GitHub/desert-fires-impact-biodiversity")
download.file(url = "https://biogeo.ucdavis.edu/data/worldclim/v2.1/hist/wc2.1_2.5m_prec_2000-2009.zip", 
              destfile = "~/GitHub/desert-fires-impact-biodiversity")
download.file(url = "https://biogeo.ucdavis.edu/data/worldclim/v2.1/hist/wc2.1_2.5m_prec_2010-2018.zip", 
              destfile = "~/GitHub/desert-fires-impact-biodiversity")
download.file(url = "https://biogeo.ucdavis.edu/data/worldclim/v2.1/hist/wc2.1_2.5m_tmin_1990-1999.zip", 
              destfile = "~/GitHub/desert-fires-impact-biodiversity")
download.file(url = "https://biogeo.ucdavis.edu/data/worldclim/v2.1/hist/wc2.1_2.5m_tmin_2000-2009.zip", 
              destfile = "~/GitHub/desert-fires-impact-biodiversity")
download.file(url = "https://biogeo.ucdavis.edu/data/worldclim/v2.1/hist/wc2.1_2.5m_tmin_2010-2018.zip", 
              destfile = "~/GitHub/desert-fires-impact-biodiversity")
download.file(url = "https://biogeo.ucdavis.edu/data/worldclim/v2.1/hist/wc2.1_2.5m_tmax_1990-1999.zip", 
              destfile = "~/GitHub/desert-fires-impact-biodiversity")
download.file(url = "https://biogeo.ucdavis.edu/data/worldclim/v2.1/hist/wc2.1_2.5m_tmax_2000-2009.zip", 
              destfile = "~/GitHub/desert-fires-impact-biodiversity")
download.file(url = "https://biogeo.ucdavis.edu/data/worldclim/v2.1/hist/wc2.1_2.5m_prec_2010-2018.zip", 
              destfile = "~/GitHub/desert-fires-impact-biodiversity")

#unzip, run for each file
system("unzip https://biogeo.ucdavis.edu/data/worldclim/v2.1/hist/wc2.1_2.5m_prec_1990-1999.zip")
  
####TEST CODE###
# #open raster layer
# 
# prec19 <- raster("raster/worldclim/prec/wc2.1_2.5m_prec_2010-01.tif")
# 
# plot(prec19)
# projection(prec19)
# 
# #create crop extent from deserts and reproject to raster crs
# crop_extent <- readOGR("shapefiles/deserts/deserts.shp") %>% 
#   spTransform(crs(prec19))
# 
# plot(crop_extent, axes = TRUE, add = TRUE)
# summary(crop_extent)
# 
# #crop
# prec19_crop <- crop(prec19, crop_extent)
# 
# plot(prec19_crop)
# 
# #save new cropped raster
# writeRaster(prec19_crop, "raster/worldclim/prec_2010-2018/cropped/prec_2019-12_crop.tif")

# string <- "wc2.1_2.5m_prec_2010-01.tif" 
# string <- sub("wc2.1_2.5m_prec_", "", string) 
# string <- strsplit(string,"-")[[1]][1]


#AUTOMATE Precipitation!!

rasterfiles <- list.files("raster/worldclim/prec/original/wc2.1_2.5m_prec_2000-2009/", pattern = "*.tif", full.name = TRUE)
reps <- (1:length(rasterfiles))

for (rep in reps) {
  
  #make new name for each file
  newfileName <- paste(sub("raster/worldclim/prec/original/wc2.1_2.5m_prec_2000-2009/wc2.1_2.5m_", "", rasterfiles[rep]), sep = "")

  year <- sub("prec_", "", newfileName) 
  year <- strsplit(year, "-")[[1]][1]
  
  #read raster layer
  raster <- raster(rasterfiles[rep])

  #create crop extent from deserts and reproject to raster crs
  crop_extent <- readOGR("shapefiles/deserts/deserts.shp") %>% 
    spTransform(crs(raster))
  
  #crop
  raster_crop <- crop(raster, crop_extent)
  
  #create new folder
  dir.create(paste("raster/worldclim/prec/cropped/", year, sep = ""))
  
  
  #save new cropped raster
  writeRaster(raster_crop, filename = paste("raster/worldclim/prec/cropped", year, newfileName, sep = "/"))
  
  }

#check raster files
env_data<- list.files(paste("raster/worldclim/tempmin/cropped/", rep, sep = ""), pattern = "*.tif", full.name = TRUE)

r.env <- list()
for(i in env_data) {
  r <- try( raster(i) )
  r.env[i] <- data.frame(name = i, nrow=nrow(r), ncol=ncol(r), res=res(r)[1],               
                         proj=proj4string(r), xmin=extent(r)[1],             
                         xmax=extent(r)[2], ymin=extent(r)[3],                     
                         ymax=extent(r)[4])
}

do.call("rbind", r.env)

##get total yearly rainfall

reps <- (2000:2009)
for (rep in reps) {
  croppedRasters <- list.files(paste("raster/worldclim/prec/cropped/", rep, sep = ""), pattern = "*.tif", full.name = TRUE)
  
  #new file name
  newtotalName <- paste("precip", rep, sep = "_")

  #create multi-layer raster of all months in one year
  rasterstack <- stack(croppedRasters)
  
  #sum up precipitation for the year
  total_precip <- calc(rasterstack, sum)
  
  writeRaster(total_precip, filename = paste("raster/worldclim/prec/totals/", newtotalName, ".tif", sep = "" ))

}


#AUTOMATE Temperature!!
rasterfiles <- list.files("raster/worldclim/tempmax/original/wc2.1_2.5m_tmax/", pattern = "*.tif", full.name = TRUE)
reps <- (1:length(rasterfiles))

for (rep in reps) {
  
  #make new name for each file
  newfileName <- paste(sub("raster/worldclim/tempmax/original/wc2.1_2.5m_tmax/wc2.1_2.5m_", "", rasterfiles[rep]), sep = "")
  
  year <- sub("tmax_", "", newfileName) 
  year <- strsplit(year, "-")[[1]][1]
  
  #read raster layer
  raster <- raster(rasterfiles[rep])
  
  #create crop extent from deserts and reproject to raster crs
  crop_extent <- readOGR("shapefiles/deserts/deserts.shp") %>% 
    spTransform(crs(raster))
  
  #crop
  raster_crop <- crop(raster, crop_extent)
  
  #create new folder
  dir.create(paste("raster/worldclim/tempmax/cropped/", year, sep = ""))
  
  
  #save new cropped raster
  writeRaster(raster_crop, filename = paste("raster/worldclim/tempmax/cropped", year, newfileName, sep = "/"))
  
}

#check raster files
env_data<- list.files(paste("raster/worldclim/tempmin/cropped/", rep, sep = ""), pattern = "*.tif", full.name = TRUE)

r.env <- list()
for(i in env_data) {
  r <- try( raster(i) )
  r.env[i] <- data.frame(name = i, nrow=nrow(r), ncol=ncol(r), res=res(r)[1],               
                         proj=proj4string(r), xmin=extent(r)[1],             
                         xmax=extent(r)[2], ymin=extent(r)[3],                     
                         ymax=extent(r)[4])
}

do.call("rbind", r.env)

reps <- (2010:2018)
for (rep in reps) {
  croppedRasters <- list.files(paste("raster/worldclim/tempmax/cropped/", rep, sep = ""), pattern = "*.tif", full.name = TRUE)
  
  #new file name
  newtotalName <- paste("tmax", rep, sep = "_")
  
  #create multi-layer raster of all months in one year
  rasterstack <- stack(croppedRasters)
  
  #get min temp for the year
  t_max_year <- max(rasterstack, na.rm = TRUE)
  
  #get average min temp for the year
  t_max_avg <- calc(rasterstack, fun = mean, na.rm = TRUE)
  
  
  writeRaster(t_max_avg, filename = paste("raster/worldclim/tempmax/average/", newtotalName, ".tif", sep = "" ))
  writeRaster(t_max_year, filename = paste("raster/worldclim/tempmax/max/", newtotalName, ".tif", sep = "" ))
  
}


#check raster files
env_data<- list.files(paste("raster/worldclim/tempmin/cropped/", rep, sep = ""), pattern = "*.tif", full.name = TRUE)

r.env <- list()
for(i in env_data) {
  r <- try( raster(i) )
  r.env[i] <- data.frame(name = i, nrow=nrow(r), ncol=ncol(r), res=res(r)[1],               
                         proj=proj4string(r), xmin=extent(r)[1],             
                         xmax=extent(r)[2], ymin=extent(r)[3],                     
                         ymax=extent(r)[4])
}

do.call("rbind", r.env)


###EXTRACT VALUES TO POINTS###

data <- st_read("shapefiles/animals_v4/animals_v4.shp") 
rasterlist <- list.files("raster/worldclim/prec/totals", pattern = "*.tif", full.name = TRUE)
head(data)
lst <- list()


for (i in 1:length(rasterlist)) {
  
  #select year
  year <- sub("raster/worldclim/prec/totals/precip_", "", rasterlist[i]) 
  year <- as.numeric(strsplit(year, ".tif")[[1]][1])
  
  #load raster
  raster <- raster(rasterlist[i])
  
  #transform & subset data. Add subset to a list
  lst[[i]] <- data[data$year==year,]
  
  #extract raster values to points
  lst[[i]]$precip_tot_mm <- raster::extract(raster, lst[[i]])
  
}

#add data that does not have ndvi values
lst[[25]] <- data[data$year > 2018,]
lst[[25]]$precip_tot_mm <- NA

#concatenate list elements into one shapefile
data_new <- do.call(rbind, lst)

data <- data_new

#save :]
dir.create("shapefiles/animals_v5")
st_write(data, "shapefiles/animals_v5/animals_v5.shp")

write.csv(data.frame(data), "data/animals_2022_v4.csv")
