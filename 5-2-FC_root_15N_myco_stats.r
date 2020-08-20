#5-2-FC_root_15N_myco_stats

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(tidyverse)
library(stargazer)

nitrogeninfo = read_csv("processeddata/isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates.csv")

# STATS FOR PAPER: Root N correlates with mycorrhiza N
rootNformycoN_linear = lm(uncolN15ppmexcess ~ mycoN15ppmexcess, data = nitrogeninfo)
plot(rootNformycoN_linear)
summary(rootNformycoN_linear)

sink("stats_tables/rootNformycoN_lm_results.html")

stargazer(rootNformycoN_linear, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE)

sink()
