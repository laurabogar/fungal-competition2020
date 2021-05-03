# 5-1-FC_hyphal_vs_myco_13C_withoutlier_statsandplot

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(tidyverse)
library(stargazer)

carboninfo = read_csv("processeddata/data_for_carbon_only_analyses.csv")

#### How well did hyphal C track myco C? ####

myhyphaelm_withoutlier = lm(log(hyphae.ppm13Cexcess) ~ log(mycoC13ppmexcess), data = carboninfo)
plot(myhyphaelm_withoutlier)
summary(myhyphaelm_withoutlier)

# LME approach: NOT using this because pairing hyphae
# with mycos in the same compartment controls for batch and
# plant effects in itself -- there's no need to also
# include this as a random effect separately.

# hyphaelme_withoutlier = lmer(hyphae.ppm13Cexcess ~ mycoC13ppmexcess + (1|Batch), data = carboninfo)
# # boundary (singular) fit. Maybe not a problem?
# plot(hyphaelme_withoutlier) #outlier is very clear
# summary(hyphaelme_withoutlier)

sink("stats_tables/hyphae13C_vs_myco13C_lm_results.html")

stargazer(myhyphaelm_withoutlier, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()
