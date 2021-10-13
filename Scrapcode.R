
#scraps of code for small fixes


#set working directory
setwd("Github/desert-fires-impact-biodiversity/")

get_obs_ID <- read_csv("gold2.csv")
firemonth <- read_csv("bronze.csv") %>% 
  select(c("fireID", "fireendmonth"))
data <- left_join(data, firemonth, by = "fireID")

write_csv(data, "goldcounts.csv")

#rename wildhorse fire to Hackberry Complex 1
data <- data %>% 
  select(-1)
data$firename[data$firename == "WILDHORSE"] = "HACKBERRY COMPLEX 1"

firemonth <- sites %>% select(c("fireID", "fire_month"))
  
#count data
sites <- read_csv("bronze.csv") 
firemonth <- sites %>% select(c("fireID", "fire_month"))
data <- read_csv("goldcounts.csv") %>% 
  left_join(., firemonth, by = "fireID")

data <- left_join(data, firemonth, by = "fireID")
#need to convert month (chr) to month(date)

# convert xlsx to csv
library("rio")
xls <- dir(pattern = "xls") #need to make the data folder your working directory
created <- mapply(convert, xls, gsub("xls", "csv", xls))
unlink(xls) # delete xlsx files


#combine all buffered data into one CSV
one <- read_csv("data/bronze-buffer-1km.csv")
two <- read_csv("data/bronze-buffer-2km_.csv")
three <- read_csv("data/bronze-buffer-3km_.csv")
five <- read_csv("data/bronze-buffer-5km_.csv")

buffer <- rbind(one, two, three, five)

buffer$buffSize_km <- buffer$buffSize_km/1000
head(buffer)
write_csv(buffer, "data/bronze-buffer.csv")

one <- read_csv("data/gold-buffer-1km.csv") %>% 
  select(-c("FID_"))
two <- read_csv("data/gold-buffer-2km_.csv")
three <- read_csv("data/gold-buffer-3km_.csv")
five <- read_csv("data/gold-buffer-5km_.csv")

write_csv(gold, "data/gold-buffer.csv")


bronze <- read_csv("data/bronze-buffer.csv")
gold <- read_csv("data/gold-buffer.csv") 

gold$buffSize_km <- gold$buffSize_km/1000

gold$buffSize_km

head(gold)
head(bronze)

#add obs ID
buffer <- buffer %>% 
  left_join(gold, by = c("species", "year", "month", "day", "decimalLongitude", "decimalLatitude", "individualCount", "basisOfRecord"))

id <- rownames(gbif_sf) 
gbif_sf <- cbind(obsID=id, gbif_sf)

head(buffer)
unique(buffer$obsID)

#add months to buffer data
bronze <- read_csv("data/bronze.csv") %>% 
  select(c(fireID,firemonth))

buffermonth <- left_join(buffer, bronze, by = "fireID")
buffermonth$firemonth
head(buffermonth)

buffer <- read_csv("data/bronze-buffer.csv")



