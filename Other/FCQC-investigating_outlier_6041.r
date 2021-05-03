# FCQC-Investigating_outlier_6041

setwd("~/Documents/Fungal competition project/fungal-competition2020/")


require(tidyverse)
require(cowplot)

nitrogeninfo = read_csv("FCdata/isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates.csv")
percent_col = read_csv("FCdata/percent_colonization_and_mass_data_by_compartment.csv")


#### Investigating the outlier ####

# Let's look at compartment 6041b, which has higher C 
# and N enrichment than the others:
# My harvest note for this plant
# indicate that it was supposed to be
# SUIPU/NM, and ended up SUIPU/THETE.
# Furthermore, on the THETE side, for which
# we have the labeling data, my notes say
# "not much fungus; looks young/new."
# I could make an argument that this
# sample appears to represent a very different,
# metabolically active developmental stage,
# and exclude it from analysis for now.

outlier = percent_col[percent_col$Plant == 6041,]
# Very low colonization on the side receiving the label,
# but took up a ton.
# I think this fungus was ACTIVELY growing and courting the plant.
# It was new on that side (should have been NM).

# Including it could make it a lot harder to understand anything about
# what is happening with the rest of the plants.
# On the other hand, I think these are real data that just
# represent a different developmental stage.

# I think it's worth including in the main paper,
# since I think it provides real biological insight,
# and the qualitative patterns are the same with or
# without it. I'll present plots excluding this point
# in the supplement, though.

rootNformycoN_plot = ggplot(data = nitrogeninfo) +
  geom_point(aes(x = mycorrhizas.APE15N,
                 y = uncolonized_roots.APE15N, 
                 color = N_level,
                 shape = mycofungus)) +
  geom_smooth(method = "lm",
              aes(x = mycorrhizas.APE15N,
                  y = uncolonized_roots.APE15N),
              color = "black",
              size = 0.5) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  ylab(bquote(atop("Uncolonized roots "^15*N, "(atom percent excess)"))) +
  xlab(expression("Mycorrhizal "^15*"N (atom percent excess)")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

nitrogeninfo_nooutlier = subset(nitrogeninfo,
                                Plant != 6041)

rootNformycoN_plot_nooutlier = ggplot(data = nitrogeninfo_nooutlier) +
  geom_point(aes(x = mycorrhizas.APE15N,
                 y = uncolonized_roots.APE15N, 
                 color = N_level,
                 shape = mycofungus)) +
  geom_smooth(method = "lm",
              aes(x = mycorrhizas.APE15N,
                  y = uncolonized_roots.APE15N),
              color = "black",
              size = 0.5) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  ylab(bquote(atop("Uncolonized roots "^15*N, "(atom percent excess)"))) +
  xlab(expression("Mycorrhizal "^15*"N (atom percent excess)")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

rootNformycoN_plot_nooutlier = ggplot(data = nitrogeninfo_nooutlier) +
  geom_point(aes(x = mycorrhizas.APE15N,
                 y = uncolonized_roots.APE15N, 
                 color = N_level,
                 shape = mycofungus)) +
  geom_smooth(method = "lm",
              aes(x = mycorrhizas.APE15N,
                  y = uncolonized_roots.APE15N),
              color = "black",
              size = 0.5) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  ylab(bquote(atop("Uncolonized roots "^15*N, "(atom percent excess)"))) +
  xlab(expression("Mycorrhizal "^15*"N (atom percent excess)")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")
