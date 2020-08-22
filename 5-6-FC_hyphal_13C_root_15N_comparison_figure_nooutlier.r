#5-6-FC_hyphal_13C_root_15N_comparison_figure_nooutlier.r

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(cowplot)
library(tidyverse)

carboninfo = read_csv("processeddata/data_for_carbon_only_analyses.csv")
nitrogeninfo = read_csv("processeddata/isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates.csv")

nitrogeninfo_nooutlier = subset(nitrogeninfo, Plant != 6041)
carboninfo_nooutlier = carboninfo[!carboninfo$hyphae.APE13C == max(carboninfo$hyphae.APE13C),] # omit outlier 6024b

#### Carbon panel: How well did hyphal C track myco C? ####

hyphalCformycoC_plot_nooutlier = ggplot(data = carboninfo_nooutlier) +
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

hyphalCformycoC_plot_nolegend = hyphalCformycoC_plot_nooutlier +
  theme(legend.position = "none")

### Nitrogen panel ####
rootNformycoN_plot_nooutlier = ggplot(data = nitrogeninfo_nooutlier) +
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

rootNformycoN_plot_nolegend = rootNformycoN_plot_nooutlier +
  theme(legend.position = "none")

# When considering N and C together, I need the dataset 
# that only includes plants/compartments that received nitrogen label.

mycoCforN = ggplot(data = nitrogeninfo_nooutlier) +
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
                        nrow = 1,
                        ncol = 3,
                        rel_widths = c(1, 1, 1.15))

save_plot("plots/Multipanel_regressions_myco_N_and_C_NOOUTLIER.jpeg",
          threepanels,
          base_aspect_ratio = 3.6)

save_plot("plots/Multipanel_regressions_myco_N_and_C_NOOUTLIER.pdf",
          threepanels,
          base_aspect_ratio = 3.6)

CforNlm_nooutlier = lm((mycoC13ppmexcess) ~ mycoN15ppmexcess, data = nitrogeninfo_nooutlier)
plot(CforNlm_nooutlier)
summary(CforNlm_nooutlier)



sink("stats_tables/myco13C_vs_myco_15N_lm_results_NOOUTLIER.html")

stargazer(CforNlm_nooutlier, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE)

sink()
