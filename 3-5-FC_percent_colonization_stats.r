# 3-5 -- FC percent colonization analysis

# Analyzing percent colonization as a function of fungal treatment, N addition, and their interactions
# lm approach

setwd("~/Documents/Fungal competition project/fungal-competition2020/")


library(tidyverse)
# library(apaTables) # I haven't been able to install this successfully
library(stargazer)

alldata = read_csv("processeddata/percent_col_and_mass_data_by_plant.csv")
nonm = subset(alldata, Fungi != "None/None")
onlywithfungi = alldata[-grep("None", alldata$Fungi),]

colonizationlm = lm(percent_col ~ N_level * Fungi, data = onlywithfungi)
summary(colonizationlm)
colanova = anova(colonizationlm)
summary(colanova)

colaov = aov(percent_col ~ N_level * Fungi, data = onlywithfungi)
TukeyHSD(colaov) # Only works if you use aov() above

# apa.reg.table(biomasslm, filename = "stats_tables/Biomass_lm_APA.doc", table.number = 2)

# Saving biomass table

sink("stats_tables/colonization_lm_table.html")

stargazer(colonizationlm, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

sink()

# Saving anova table

sink("stats_tables/colonization_anova_table.html")

stargazer(colanova, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE)

sink()
