# FC - hyphal 13C tracks myco 13C

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(lme4)
library(lmerTest)

together = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly.csv")


#### How well did hyphal C track myco C? ####

carboninfo_nooutlier = carboninfo[!carboninfo$hyphae.APE13C == max(carboninfo$hyphae.APE13C),] # omit outlier 6024b, which is 1216.8 ppm excess 13C

hyphalCformycoC_plot = ggplot(data = carboninfo_nooutlier) +
  geom_point(aes(x = mycoC13ppmexcess,
                 y = hyphae.ppm13Cexcess, 
                 color = N_level,
                 shape =compartment_fungus)) +
  geom_smooth(method = "lm", aes(x = mycoC13ppmexcess,
                                 y = hyphae.ppm13Cexcess),
              color = "black",
              size = 0.5) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  ylab(expression("Hyphal "^13*"C (ppm excess)")) +
  xlab(expression("Mycorrhizal "^13*"C (ppm excess)")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

# STATS for supplement: hyphae C tracks myco C, excluding outlier  
myhyphaelm = lm((hyphae.ppm13Cexcess) ~ mycoC13ppmexcess, data = carboninfo_nooutlier)
plot(myhyphaelm)
summary(myhyphaelm)

save_plot("plots/Regression_hyphal_C_for_myco_C.pdf", 
          hyphalCformycoC_plot,
          ncol = 1,
          base_aspect_ratio = 1.4)

# Doing this with just Tt yields almost exactly the same result.

# INCLUDING OUTLIER

hyphalCformycoC_plot_withoutlier = ggplot(data = carboninfo) +
  geom_point(aes(x = mycorrhizas.APE13C,
                 y = hyphae.APE13C, 
                 color = N_level,
                 shape =compartment_fungus)) +
  geom_smooth(method = "lm", aes(x = mycorrhizas.APE13C,
                                 y = hyphae.APE13C),
              color = "black",
              size = 0.5) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  ylab(expression("Hyphal "^13*"C (atom percent excess)")) +
  xlab(expression("Mycorrhizal "^13*"C (atom percent excess)")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

# Stats with outlier
myhyphaelm_withoutlier = lm((hyphae.ppm13Cexcess) ~ mycoC13ppmexcess, data = carboninfo)
plot(myhyphaelm_withoutlier)
summary(myhyphaelm_withoutlier)

# hyphaelme_withoutlier = lmer(hyphae.APE13C ~ mycorrhizas.APE13C + (1|Batch), data = carboninfo)
# # boundary (singular) fit. Maybe not a problem?
# plot(hyphaelme_withoutlier) #outlier is very clear
# summary(hyphaelme_withoutlier)

sink("stats_tables/hyphae13C_vs_myco13C_withoutlier_lm_results.html")

stargazer(myhyphaelm_withoutlier, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE)

sink()



save_plot("plots/Regression_hyphal_C_for_myco_C_with_outlier.pdf", 
          hyphalCformycoC_plot_withoutlier,
          ncol = 1,
          base_aspect_ratio = 1.4)
