# FC - hyphal 13C tracks myco 13C

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

require(tidyverse)
require(cowplot)

together = read_csv("FCdata/isotope_and_plant_metadata_with_competition_coded_clearly.csv")


#### How well did hyphal C track myco C? ####
carboninfo = subset(together, compartment_fungus != "None" &
                      compartment_fungus != "MIXED" &
                      compartment_fungus != "OTHER")

carboninfo = carboninfo[!is.na(carboninfo$hyphae.APE13C),]
carboninfo = carboninfo[!is.na(carboninfo$mycorrhizas.APE13C),]
# For any analysis involving ONLY C, it makes sense to include
# 1) carbon data from root compartments that
# didn't necessarily receive N label (20 of these, 22 labeled compartments here with hyphal C info), and
# 2) data from plants with failed splits (two of these, 40 successfully split)



justTt = subset(carboninfo, compartment_fungus == "Tt")

carboninfo_nooutlier = carboninfo[!carboninfo$hyphae.APE13C == max(carboninfo$hyphae.APE13C),] # omit outlier 6024b


hyphalCformycoC_plot = ggplot(data = carboninfo_nooutlier) +
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

# STATS for supplement: hyphae C tracks myco C, excluding outlier  
myhyphaelm = lm((hyphae.APE13C) ~ mycorrhizas.APE13C, data = carboninfo_nooutlier)
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


myhyphaelm_withoutlier = lm((hyphae.APE13C) ~ mycorrhizas.APE13C, data = carboninfo)
plot(myhyphaelm_withoutlier)
summary(myhyphaelm_withoutlier)

save_plot("plots/Regression_hyphal_C_for_myco_C_with_outlier.pdf", 
          hyphalCformycoC_plot_withoutlier,
          ncol = 1,
          base_aspect_ratio = 1.4)
