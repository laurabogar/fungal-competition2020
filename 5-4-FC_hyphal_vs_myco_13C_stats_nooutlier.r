#5-4-FC_hyphal_vs_myco_13C_without_outlier_stats

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(stargazer)

carboninfo = read_csv("processeddata/data_for_carbon_only_analyses.csv")
carboninfo_nooutlier = carboninfo[!carboninfo$hyphae.APE13C == max(carboninfo$hyphae.APE13C),] # omit outlier 6024b

#### How well did hyphal C track myco C? ####

myhyphaelm_nooutlier = lm((hyphae.ppm13Cexcess) ~ mycoC13ppmexcess, data = carboninfo_nooutlier)
plot(myhyphaelm_nooutlier)
summary(myhyphaelm_nooutlier)

myhyphaelm_log_nooutlier = lm((hyphae.ppm13Cexcess) ~ log(mycoC13ppmexcess), data = carboninfo_nooutlier)
plot(myhyphaelm_log_nooutlier)
summary(myhyphaelm_log_nooutlier)

myhyphaelm_loglog_nooutlier = lm(log(hyphae.ppm13Cexcess) ~ log(mycoC13ppmexcess), data = carboninfo_nooutlier)
plot(myhyphaelm_loglog_nooutlier) # Log-log plots look the best.
summary(myhyphaelm_loglog_nooutlier)

sink("stats_tables/hyphae13C_vs_myco13C_comparingthreemodels_nooutlier.html")

stargazer(myhyphaelm_nooutlier, myhyphaelm_log_nooutlier, myhyphaelm_loglog_nooutlier, 
          type = "html",
          align = TRUE,
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          no.space = TRUE)

sink()

# LME approach: NOT using this because pairing hyphae
# with mycos in the same compartment controls for batch and
# plant effects in itself -- there's no need to also
# include this as a random effect separately.

sink("stats_tables/hyphae13C_vs_myco13C_lm_results_NOOUTLIER.html")

stargazer(myhyphaelm_nooutlier, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE)

sink()
