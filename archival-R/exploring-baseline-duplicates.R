#issue with 2021 pull and 2022 pull
#2021 has a lot more records than 2022

#repeat occurrences in both 2021 and 2022
library(tidyverse)
library(arsenal)

setwd("~/GitHub/desert-fires-impact-biodiversity")

#most current (2022)
df2 <- read_csv("data/animals-baseline-2022.csv") %>% 
  select(-c(obsID, kingdom, family, phylum, order, 
            genus, infraspecificEpithet, scientificName, 
            eventDate, basisOfRecord, institutionCode)) %>% 
  rename(day = obsDay, month = obsMonth, year = obsYear, 
         dcmlLng = decimalLongitude, dcmlLtt = decimalLatitude, indvdlC = individualCount) %>% 
  filter(year > 1994, year < 2021) 

#pulled 2021
df1 <- bind_rows((read_csv("C:/Users/marin/OneDrive/Documents/GitHub/desert-fires-impact-biodiversity/data/archive/data/animals-noburn_MOJ.csv")),
                 (read_csv("C:/Users/marin/OneDrive/Documents/GitHub/desert-fires-impact-biodiversity/data/archive/data/animals-noburn_SJ.csv")),
                 (read_csv("C:/Users/marin/OneDrive/Documents/GitHub/desert-fires-impact-biodiversity/data/archive/data/animals-noburn_SON.csv"))) %>% 
  select(-c(FID, obsID, class, order, bssOfRc)) %>% 
  filter(year > 1994, year < 2021)


## address duplicates
dup_df2 <- df2
dup_df2$is_duplicated <- duplicated(df2)
dup_df2 <- dup_df2 %>%   
  filter(is_duplicated == TRUE)



#compare dfs
comparedf(df1, df2)
# Shared: 8 non-by variables and 37284 observations.
# Not shared: 0 variables and 1092 observations.

#filters rows in df1 that are also in df2
shared <- semi_join(df1, df2) #df1 has 2329 rows that df2 also has

#filter rows in df1 that are also in df 2
shared <- semi_join(df2, df1) #df2 has 4853 rows that df1 also has

#filter rows in df1 that are not found in df2
not_shared <- setdiff(df1,df2)

#filter rows in df2 that are not found in df1
not_shared <- setdiff(df2,df1)


uniq_df1 <- anti_join(df1, df2, by = colnames(df1))
uniq_df2 <- anti_join(df2, df1, by = colnames(df1))

all <- bind_rows(df1, df2)
distinct_all <- distinct(all)

#return all rows in 'all' df that do not have a matching row in df1
uniq_df2 <- anti_join(distinct_all, df1, by = colnames(all))
#return all rows in 'all' df that do not have a matching row in df2
uniq_df1 <- anti_join(distinct_all, df2, by = colnames(all))

no_match <- anti_join(df2, df1, by = colnames(all))


## address duplicates
dup_df2$is_duplicated <- duplicated(df2)
dup_df2 <- dup_df2 %>%   
  filter(is_duplicated == TRUE)

test_df <- data.frame(row1 = c(1, 1, 1, 2, 2, 2, 3, 3, 3), row2 = c("a", "a", "b", "a", "b", "c", "c", "c", "c"))
test_df$is_dup <- duplicated(test_df) 
test_df <- test_df %>% filter(is_dup == TRUE) %>% distinct(test_df)
