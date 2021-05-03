# 3.5 -- FC biomass analysis

# Analyzing plant biomass as a function of fungal treatment, N addition, percent colonization, and their interactions
# lm approach

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(tidyverse)
# library(apaTables) # I haven't been able to install this successfully
library(stargazer)

alldata = read_csv("processeddata/percent_col_and_mass_data_by_plant.csv")

biomasslm = lm(total_biomass ~ N_level * Fungi * percent_col, data = alldata)
summary(biomasslm)
bioanova = anova(biomasslm)


# apa.reg.table(biomasslm, filename = "stats_tables/Biomass_lm_APA.doc", table.number = 2)

# Saving biomass table

sink("stats_tables/biomass_lm_table.html")

stargazer(biomasslm, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          no.space = TRUE)

sink()

# Saving anova table

sink("stats_tables/biomass_anova_table.html")

stargazer(bioanova, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()
