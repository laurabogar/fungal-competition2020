# 3.5 -- FC biomass analysis

# Analyzing plant biomass as a function of fungal treatment, N addition, percent colonization, and their interactions
# lm approach

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(tidyverse)

alldata = read_csv("processeddata/percent_col_and_mass_data_by_plant.csv")

biomasslm = lm(total_biomass ~ N_level * Fungi * percent_col, data = alldata)
summary(biomasslm)
bioanova = anova(biomasslm)
