#5-3-FC_hyphal_13C_root_15N_myco_comparison_figure

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(cowplot)
library(tidyverse)

carboninfo = read_csv("processeddata/data_for_carbon_only_analyses.csv")
nitrogeninfo = read_csv("processeddata/isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates.csv")


#### Carbon panel: How well did hyphal C track myco C? ####

hyphalCformycoC_plot = ggplot(data = carboninfo) +
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

hyphalCformycoC_plot_nolegend = hyphalCformycoC_plot +
  theme(legend.position = "none")
# save_plot("plots/Regression_hyphal_C_for_myco_C_with_outlier.pdf", 
#           hyphalCformycoC_plot_withoutlier,
#           ncol = 1,
#           base_aspect_ratio = 1.4)
# 
# ggsave("plots/Regression_hyphal_C_for_myco_C_with_outlier.jpeg", 
#        plot = hyphalCformycoC_plot_withoutlier,
#        device = "jpeg",
#        width = 7, height = 5, units = "in")

### Nitrogen panel ####
rootNformycoN_plot = ggplot(data = nitrogeninfo) +
  geom_point(aes(x = mycoN15ppmexcess,
                 y = uncolN15ppmexcess, 
                 color = N_level,
                 shape = mycofungus)) +
  geom_smooth(method = "lm",
              aes(x = mycoN15ppmexcess,
                  y = uncolN15ppmexcess),
              color = "black",
              size = 0.5) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  ylab(bquote(atop("Uncolonized roots "^15*N, "(ppm excess)"))) +
  xlab(expression("Mycorrhizal "^15*"N (ppm excess)")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

rootNformycoN_plot_nolegend = rootNformycoN_plot +
  theme(legend.position = "none")

NforNandCforC = plot_grid(rootNformycoN_plot_nolegend, hyphalCformycoC_plot,
                          labels = c("A", "B"),
                          align = "h",
                          rel_widths = c(1, 1.15))

# save_plot("plots/Multipanel_regressions_NforN_and_CforC.pdf",
#           NforNandCforC,
#           ncol = 2,
#           base_aspect_ratio = 1.4)
# 
# save_plot("plots/Multipanel_regressions_NforN_and_CforC.jpeg",
#           NforNandCforC,
#           ncol = 2,
#           base_aspect_ratio = 1.4)

# When considering N and C together, I need the dataset 
# that only includes plants/compartments that received nitrogen label.

mycoCforN = ggplot(data = nitrogeninfo) +
  geom_point(aes(x = mycoN15ppmexcess,
                 y = mycoC13ppmexcess, 
                 color = N_level,
                 shape = compartment_fungus)) +
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              aes(x = mycoN15ppmexcess,
                  y = mycoC13ppmexcess),
              color = "black",
              size = 0.5) +
  ylab(bquote(atop("Mycorrhizal "^13*"C", "(ppm excess)"))) +
  xlab(bquote(atop("Mycorrhizal "^15*"N (ppm excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

threepanels = plot_grid(rootNformycoN_plot_nolegend, 
                        hyphalCformycoC_plot_nolegend,
                        mycoCforN,
                          labels = c("A", "B", "C"),
                          align = "h",
                          rel_widths = c(1, 1, 1.15))

save_plot("plots/Multipanel_regressions_myco_N_and_C.jpeg",
          threepanels,
          ncol = 3,
          base_aspect_ratio = 1.4)

save_plot("plots/Multipanel_regressions_myco_N_and_C.pdf",
          threepanels,
          ncol = 3)
