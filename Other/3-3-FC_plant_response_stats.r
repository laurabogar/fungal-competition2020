# 3-3-FC_plant_response stats

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(tidyverse)
library(stargazer)

alldata = read_csv("processeddata/percent_col_and_mass_data_by_plant.csv")


#### PLANT RESPONSE TO COLONIZATION ####

noNMplants = subset(alldata, Fungi != "None/None")

#### Bringing together colonization and biomass in linear model ####

responselm = lm(plant_response ~ N_level * Fungi * percent_col, data = noNMplants)
summary(responselm)
responseanova = anova(responselm)

# Saving biomass table

sink("stats_tables/plantgrowthresponse_lm_table.html")

stargazer(responselm, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          no.space = TRUE)

sink()

# Saving anova table

sink("stats_tables/plantgrowthresponse_anova_table.html")

stargazer(responseanova, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()

### t tests: Do these differ from zero? ###

# High N:
hightx = subset(noNMplants, N_level == "High")
t.test(hightx$plant_response[hightx$Fungi == "Sp/None"]) # not enough observations
t.test(hightx$plant_response[hightx$Fungi == "Sp/Sp"]) # p = 0.1185
t.test(hightx$plant_response[hightx$Fungi == "Tt/Sp"]) # p = 0.02562
t.test(hightx$plant_response[hightx$Fungi == "Tt/None"]) # p = 0.1547
t.test(hightx$plant_response[hightx$Fungi == "Tt/Tt"]) # p = 7.478e-11

# Low N:
lowtx = subset(noNMplants, N_level == "Low")
t.test(lowtx$plant_response[lowtx$Fungi == "Sp/None"]) # p = 0.1654
t.test(lowtx$plant_response[lowtx$Fungi == "Sp/Sp"]) # p = 0.6937
t.test(lowtx$plant_response[lowtx$Fungi == "Tt/Sp"]) # p = 0.2881
t.test(lowtx$plant_response[lowtx$Fungi == "Tt/None"]) # p = 0.399
t.test(lowtx$plant_response[lowtx$Fungi == "Tt/Tt"]) # p = 0.1192

# Alpha should be 0.05/10 tests = 0.005, per Bonferroni correction. So the only significant 
# non-zero growth effect difference is Tt/Tt in high N




