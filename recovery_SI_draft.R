Sim <- SimilarityPair(recoverytsf_list[[1]][[7]],"incidence_freq", nboot = 200)

# #estimated Sorensens
# Sim[["estimated_richness"]][[1,1]]
# #estimated Sorensens lower limit
# Sim[["estimated_richness"]][[1,3]]
# #estimated Sorensens upper limit
# Sim[["estimated_richness"]][[1,4]]


#estimated ChaoSorensens-abundance
Sim[["estimated_relative"]][[5]]
#estimated Sorensens lower limit
Sim[["estimated_relative"]][[5, 3]]
#estimated Sorensens upper limit
Sim[["estimated_relative"]][[5, 4]]

#Number of observed species in Community 1
Sim[["info"]][["D1"]]
#Number of observed species in Community 2
Sim[["info"]][["D2"]]

## create df of species richness, SI, CI


desert1 <- c()
yearsSinceFire <- c()
SpRich_burned <- c()
SpRich_control <- c()
SI <- c()
SI_CIlower <- c()
SI_CIupper <- c()
error_desert <- c()
error_tsf <- c()

for (d in 1:3) {
  
  for (t in 1:length(recoverytsf_list[[d]])){
    tryCatch(
      #try to do this
      {
    
      Sim <- SimilarityPair(recoverytsf_list[[d]][[t]],"incidence_freq", nboot = 200)
      
      desert1 <- append(desert1, names(recoverytsf_list[d]))
      yearsSinceFire <- append(yearsSinceFire, as.integer(names(recoverytsf_list[[d]][t])))
      SpRich_burned <- append(SpRich_burned, Sim[["info"]][["D1"]])
      SpRich_control <- append(SpRich_control, Sim[["info"]][["D2"]])
      SI <- append(SI, Sim[["estimated_relative"]][[5]])
      SI_CIlower <- append(SI_CIlower, Sim[["estimated_relative"]][[5, 3]])
      SI_CIupper <- append(SI_CIupper, Sim[["estimated_relative"]][[5, 4]])
      
      
      },
      #if an error occurs, say error
      
      error=function(e) {
        message(paste('An Error Occurred for:', desert[d], '| TSF:', as.integer(names(recoverytsf_list[[d]][t]))))
        print(e)
        error_desert <- append(error_desert, desert[d])
        error_tsf <- append(error_tsf, as.integer(names(recoverytsf_list[[d]][tsf])))
        }
      
    )
  }
    
  }
df_tsf <- data.frame(desert1, yearsSinceFire, SpRich_burned, SpRich_control, SI, SI_CIlower, SI_CIupper) %>% 
  rename(desert = desert1)

df_error <- data.frame(error_desert, error_tsf)

#45 errors

ggplot(df_tsf, aes(x=yearsSinceFire, y=SI, color = desert)) +
  geom_point() +
  # geom_errorbar(aes(ymin=SI_CIlower, ymax=SI_CIupper), color = "black", width=.2) +
  geom_smooth(method = "lm", se = FALSE)

