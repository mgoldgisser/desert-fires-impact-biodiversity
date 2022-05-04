#urban areas in CA

library(tidyverse)
library(sf)
library(sp)
library(rgdal)
library(raster)

setwd("~/GitHub/desert-fires-impact-biodiversity")

##download the following links for road data

# download.file("http://www2.census.gov/geo/tiger/TIGER2013/PRISECROADS/tl_2013_06_prisecroads.zip", "C:/Users/marin/Downloads/ca-rds")
# download.file("http://www2.census.gov/geo/tiger/TIGER2013/PRISECROADS/tl_2013_32_prisecroads.zip", "C:/Users/marin/Downloads/nv-rds")
# download.file("http://www2.census.gov/geo/tiger/TIGER2013/PRISECROADS/tl_2013_04_prisecroads.zip", "C:/Users/marin/Downloads/az-rds")

##download the following links for urban area
# download.file(https://www2.census.gov/geo/tiger/TIGER2018/UAC/tl_2018_us_uac10.zip)

baseline <- st_read("shapefiles/animals-baseline-2022-v2")
plot(baseline)

#attach yearly precipitation to obs

rlist <- list.files(path = "raster/worldclim/prec/totals", pattern = "*.tif", full.names = TRUE)
r <- raster(rlist)

rtest <- raster("raster/worldclim/prec/totals/precip_1995.tif")
x <- extract(rtest, baseline)

pointstest <- cbind(baseline, x)

