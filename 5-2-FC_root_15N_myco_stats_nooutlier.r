#5-2-FC_root_15N_myco_stats

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(tidyverse)

nitrogeninfo = read_csv("processeddata/isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates.csv")

nitrogeninfo_nooutlier = subset(nitrogeninfo,
                                Plant != 6041)


# STATS FOR PAPER: Root N correlates with mycorrhiza N
rootNformycoN_linear = lm(uncolonized_roots.APE15N ~ mycorrhizas.APE15N, data = nitrogeninfo_nooutlier)
plot(rootNformycoN_linear)
summary(rootNformycoN_linear)

save_plot("plots/Regression_root_N_for_myco_N.pdf", 
          rootNformycoN_plot, 
          base_aspect_ratio = 1.4)

nitrogeninfo[nitrogeninfo$mycorrhizas.APE15N == max(nitrogeninfo$mycorrhizas.APE15N),]
# Plant 6041b is being weird AGAIN. Got a TON of C
# in the hyphae, and a TON of N at the myco.

rootNformycoN_log = glm(forced.uncolonized.APE15N ~ log(forced.mycorrhizas.APE15N), data = nitrogeninfo)
plot(rootNformycoN_log) # very hump-shaped residuals

# INCLUDING OUTLIER:

rootNformycoN_plot_withoutlier = ggplot(data = nitrogeninfo) +
  geom_point(aes(x = mycorrhizas.APE15N,
                 y = uncolonized_roots.APE15N, 
                 color = N_level,
                 shape =compartment_fungus)) +
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

# STATS FOR PAPER: root N corresponds to myco N including outlier

rootNformycoN_linear_withoutlier = lm(uncolonized_roots.APE15N ~ mycorrhizas.APE15N, data = nitrogeninfo)
plot(rootNformycoN_linear_withoutlier)
summary(rootNformycoN_linear_withoutlier)

rootNformycoN_log_withoutlier = lm(uncolonized_roots.APE15N ~ log(forced.mycorrhizas.APE15N), data = nitrogeninfo)
plot(rootNformycoN_log_withoutlier) # not as good as linear
summary(rootNformycoN_log_withoutlier) # also not as good as linear

save_plot("plots/Regression_root_N_for_myco_N_with_outlier.pdf", 
          rootNformycoN_plot_withoutlier, 
          base_aspect_ratio = 1.4)

