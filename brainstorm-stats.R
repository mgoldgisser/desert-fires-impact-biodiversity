#stats brainstorm
library(tidyverse)

sites <- read_csv("data/bronze.csv") %>% 
  select(c("fireID", "firename", "firemonth", "firesize", "n_observations"))
buffer <- read_csv("data/bronze-buffer.csv")
datafire <- read_csv("data/gold2.csv") %>% 
  select(c("obsID", "species", "year", "month", "indvdlC", "fireID", "fireyear", "firesize")) %>% 
  rename(obsYear = year)
databuff <- read_csv("data/gold-buffer.csv")%>% 
  select(c("obsID", "species", "year", "month", "indvdlC", "fireID", "fireyear", "firemonth", "firesize", "buffSize_km", "buffArea_acres")) %>% 
  rename(buffarea = buffArea_acres, buffdist = buffSize_km, obsYear = year)

#add fire month to datafire
fmonth <- sites %>% 
  select(c("fireID", "firemonth"))

datafire <- left_join(datafire, fmonth, by = "fireID")

#remove NA for individual counts
datafire$indvdlC[is.na(datafire$indvdlC)] = 0
databuff$indvdlC[is.na(databuff$indvdlC)] = 0

head(datafire)
head(databuff)

#add column comparing observation year to fire year, column = "postingvalue" 
#and whether obs occurred before or after fire = pre_or_post; SY means same year
#n_years = absolute value of years between observation and fire
datafire <- mutate(datafire, postingvalue = obsYear - fireyear, 
                   pre_or_post = case_when(postingvalue > 0 ~ "post", 
                                           postingvalue < 0 ~ "pre", 
                                           postingvalue == 0 ~ "SY"), 
                   n_years = abs(postingvalue))


#find pre- post- for observations in the same year and replace SY in datafire with pre or post
SY <- filter(datafire, postingvalue == 0) %>% 
  mutate(postingmonth = month - firemonth,
    pre_or_post = case_when(postingmonth >0 ~ "post",
                            postingmonth < 0 ~ "pre",
                            postingvalue == 0 ~ "SM"))

datafire$pre_or_post <- ifelse(datafire$pre_or_post == "SY", SY$pre_or_post, datafire$pre_or_post)

datafire$pre_or_post

#do the same for buffered area
databuff <- mutate(databuff, postingvalue = obsYear - fireyear, 
                   pre_or_post = case_when(postingvalue > 0 ~ "post", 
                                           postingvalue < 0 ~ "pre", 
                                           postingvalue == 0 ~ "SY"), 
                   n_years = abs(postingvalue))
  

#find pre- post- for observations in the same year and replace SY in datafire with pre or post

id <- rownames(databuff) 
databuff <- cbind(id=id, databuff) %>% 
  na.exclude

SY <- filter(databuff, postingvalue == 0) %>% 
  mutate(postingmonth = month - firemonth,
         pre_or_post = case_when(postingmonth >0 ~ "post",
                                 postingmonth < 0 ~ "pre",
                                 postingmonth == 0 ~ "SM")) %>% 
  na.exclude()


#retired code - databuff$pre_or_post <- ifelse(databuff$pre_or_post == "SY", SY$pre_or_post, databuff$pre_or_post)

databuff$pre_or_post[match(databuff$id, SY$id)] <- SY$pre_or_post

#some data had observations within the same month of the fire, I manually added pre/post to those
#code below to prepare the csv and get the data rejoined

SM <- filter(SY, pre_or_post == "SM") %>% 
  na.exclude()

test <- databuff %>% 
  select(c("id", "obsID", "species", "fireID", "obsYear", "fireyear", "month", "day")) %>% 
  filter(obsID %in% SM$obsID) %>% 
  na.exclude()

write_csv(test, "data/temp.csv")
##edit csv file manually here

test <- read_csv("data/temp.csv")

SM$pre_or_post <- ifelse(SM$pre_or_post == "SM", test$pre_or_post, NA)

#make some data viz to explore

ggplot(datafire, aes(po))
