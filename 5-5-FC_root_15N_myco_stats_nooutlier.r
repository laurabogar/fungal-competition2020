#5-2-FC_root_15N_myco_stats

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(tidyverse)

nitrogeninfo = read_csv("processeddata/isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates.csv")

nitrogeninfo_nooutlier = subset(nitrogeninfo, Plant != 6041)

# STATS FOR PAPER: Root N correlates with mycorrhiza N
rootNformycoN_linear = lm(nmN15ppmexcess ~ mycoN15ppmexcess, data = nitrogeninfo_nooutlier)
plot(rootNformycoN_linear) # Looks kinda weird, we'll see how this goes.
summary(rootNformycoN_linear)

rootNformycoN_log_linear = lm(nmN15ppmexcess ~ log(forced.mycorrhizas.N15ppmexcess), data = nitrogeninfo_nooutlier)
plot(rootNformycoN_log_linear) 
summary(rootNformycoN_linear)

rootNformycoN_loglog_linear = lm(log(forced.uncolonized.N15ppmexcess) ~ log(forced.mycorrhizas.N15ppmexcess), data = nitrogeninfo_nooutlier)
plot(rootNformycoN_loglog_linear) # data point 45 looks really weird here.
summary(rootNformycoN_linear)

sink("stats_tables/root15N_vs_myco15N_lm_results_NOOUTLIER.html")

stargazer(rootNformycoN_loglog_linear, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()
