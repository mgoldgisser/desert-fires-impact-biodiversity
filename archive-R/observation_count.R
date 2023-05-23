
animals$institutionalCode <- factor(animals$institutionalCode, levels = c("CLO-eBird",  
                                                                          "iNaturalist", 
                                                                          "PRBO", 
                                                                          "CLO-GBBC", 
                                                                          "Xeno-Canto", 
                                                                          "eButterfly", 
                                                                          "CLO-ML", 
                                                                          "Observation.org", 
                                                                          "LEPSOC",
                                                                          "naturgucker"))

ggplot(filter(animals, institutionalCode != "CLO-eBird"), aes(x=obsYear, fill = institutionalCode)) + 
  geom_bar() +
  scale_fill_manual(values = c("#a6cee3", 
                             "#1f78b4", 
                             "#fecc5c",
                             "#33a02c",
                             "#636363",
                             "#e41a1c",
                             "#756bb1",
                             "#feebe2",
                             "#b3e2cd")) +
  ylab("Count") +
  xlab("Year") +
  labs(fill = "Dataset") +
  theme_bw() +
  facet_grid(vars(treatmentGroup), scales = "free_y")
        

ggplot(animals, aes(x=obsYear, fill = institutionalCode)) + 
  geom_bar() +
  scale_fill_manual(values = c("#8856a7",
                               "#a6cee3", 
                               "#1f78b4", 
                               "#fecc5c",
                               "#33a02c",
                               "#636363",
                               "#e41a1c",
                               "#756bb1",
                               "#feebe2",
                               "#b3e2cd")) +
  ylab("Count") +
  xlab("Year") +
  labs(fill = "Dataset") +
  theme_bw()+
  facet_grid(vars(treatmentGroup), scales = "free_y")


unique(animals$obsYear)
