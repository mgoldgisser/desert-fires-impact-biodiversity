#stats brainstorm

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

mutate(offset = obsYear - fireyear,
       posting = case_when(offset > 0 ~ "after", offset < 0 ~ "before", offset == 0 ~ "SY"), n_years = abs(offset)) %>% 
  na.exclude() %>% 
  distinct(fireyear, .keep_all= TRUE)

#find pre- post- for observations in the same year and replace SY in datafire with pre or post
SY <- filter(datafire, postingvalue == 0) %>% 
  mutate(postingmonth = month - firemonth,
    pre_or_post = case_when(postingmonth >0 ~ "post",
                            postingmonth < 0 ~ "pre",
                            postingvalue == 0 ~ "SM"))

datafire$pre_or_post <- ifelse(datafire$pre_or_post == "SY", SY$pre_or_post, datafire$pre_or_post)

datafire$posting
