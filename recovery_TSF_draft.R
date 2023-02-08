## Similarity index to compare recovery::TSF

### Use packages and data from desert-fires-impact-biodiversity github index.Rmd

burned <- animals %>% 
  filter(treatmentGroup == "burned") %>% 
  merge(.,sites[,c("fireID", "fireMonth")], by = "fireID") %>% #NEED to add fire months
  mutate(yearsSinceFire = obsYear - fireYear, 
         pre_or_post = case_when(yearsSinceFire > 0 ~ "post",
                                 yearsSinceFire < 0 ~ "pre",
                                 yearsSinceFire == 0 ~ "SY"),
         n_years = abs(yearsSinceFire))

#find pre- post- for observations in the same year and replace SY in datafire with pre or post
SY <- filter(burned, yearsSinceFire == 0) %>% 
  mutate(postingmonth = obsMonth - fireMonth,
         pre_or_post = case_when(postingmonth >0 ~ "post",
                                 postingmonth < 0 ~ "pre",
                                 yearsSinceFire == 0 ~ "SM"))

burned$pre_or_post <- ifelse(burned$pre_or_post == "SY", SY$pre_or_post, burned$pre_or_post)
unique(burned$pre_or_post) #check if any SM (same month); would need to further categorize

burned <- burned %>% 
  filter(pre_or_post == "post") %>% 
  dplyr::select(c("desert", "obsYear", "fireYear", "species",
                  "yearsSinceFire")) %>% 
  mutate(count = 1) %>% 
  mutate(yearsSinceFire = ifelse(yearsSinceFire == 0, 1, yearsSinceFire))

### get fire years to make control
fireYear <- burned %>% 
  group_by(fireYear, desert) %>% 
  summarize(total = n()) %>% 
  dplyr::select(-total) 

### make control
control_test <- burned %>% group_by(desert, fireYear) %>% 
  summarize(maxTime = max(yearsSinceFire)) %>% 
  filter(maxTime != 0)

### get species info for control
control_species <- animals %>% 
  filter(treatmentGroup == "never burned", obsYear > 1999) %>% 
  group_by(desert, obsYear, species) %>%
  summarize(count = n()) 

### how many species each year in control?

control_sp_total <- control_species %>% 
  group_by(desert, obsYear) %>%
  summarize(count = n())


  
## vectors for control df
yearsSinceFire <- c()
desert <- c()
fireYear <- c()


for (i in 1:nrow(control)) {
  yearsSinceFire <-  c(yearsSinceFire, rep(1:control$maxTime[i], 1))
}

for (i in 1:nrow(control)) {
  desert <-  c(desert, rep(control$desert[i], each = control$maxTime[i]))
}

for (i in 1:nrow(control)) {
  fireYear <-  c(fireYear, rep(control$fireYear[i], each = control$maxTime[i]))
}


control <- data.frame(desert, fireYear, yearsSinceFire) %>% 
  mutate(obsYear = fireYear + yearsSinceFire) %>% 
  filter(obsYear < 2021) %>% 
  left_join(., dplyr::select(control_species, -count), by = c("obsYear", "desert"))

## make df with just desert, TSF, and species for both burned and control sites

burned <- burned %>% group_by(desert, yearsSinceFire, species) %>% 
  summarize(count = n()) %>% 
  dplyr::select(-count)

control <- control %>% group_by(desert, yearsSinceFire, species) %>% 
  summarize(count = n()) %>% 
  dplyr::select(-count)

unique_b <- anti_join(burned, control) %>% 
  group_by(desert, yearsSinceFire) %>% 
  summarize(unique_b = n())

unique_c <- anti_join(control, burned) %>% 
  group_by(desert, yearsSinceFire) %>% 
  summarize(unique_c = n())

shared <- anti_join(burned, anti_join(burned, control)) %>% 
  group_by(desert, yearsSinceFire) %>% 
  summarize(shared = n())

tsf_df <- full_join(shared, unique_b, by = c("desert", "yearsSinceFire")) %>% 
  full_join(., unique_c, by = c("desert", "yearsSinceFire")) %>% 
  mutate(unique_b = ifelse(is.na(unique_b), 0, unique_b),
         shared = ifelse(is.na(shared), 0, shared),
         unique_c = ifelse(is.na(unique_c), 0, unique_c),
         si = (2*shared)/(unique_b + unique_c))

ggplot(data = tsf_df, aes(x = yearsSinceFire, y = si, color = desert)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)


