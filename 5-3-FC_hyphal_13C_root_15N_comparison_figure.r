#5-3-FC_hyphal_13C_root_15N_myco_comparison_figure

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(cowplot)
library(tidyverse)
library(stargazer)

carboninfo = read_csv("processeddata/data_for_carbon_only_analyses.csv")
nitrogeninfo = read_csv("processeddata/isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates.csv")


#### Carbon panel: How well did hyphal C track myco C? ####

hyphalCformycoC_plot = ggplot(data = carboninfo) +
  geom_point(aes(x = mycologC13,
                 y = log(hyphae.ppm13Cexcess), 
                 color = N_level,
                 shape =compartment_fungus)) +
  geom_smooth(method = "lm", aes(x = mycologC13,
                                 y = log(hyphae.ppm13Cexcess)),
              color = "black",
              size = 0.5) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  ylab(expression("Hyphal "^13*"C (ln ppm excess)")) +
  xlab(expression("Mycorrhizal "^13*"C (ln ppm excess)")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
  # geom_abline(intercept = 0, slope = 1, linetype = "dashed")

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
  geom_point(aes(x = mycologN15,
                 y = nmlogN15, 
                 color = N_level,
                 shape = mycofungus)) +
  geom_smooth(method = "lm",
              aes(x = mycologN15,
                  y = nmlogN15),
              color = "black",
              size = 0.5) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  ylab(bquote(atop("Uncolonized roots "^15*N, "(ln ppm excess)"))) +
  xlab(expression("Mycorrhizal "^15*"N (ln ppm excess)")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
  # geom_abline(intercept = 0, slope = 1, linetype = "dashed")

rootNformycoN_plot_nolegend = rootNformycoN_plot +
  theme(legend.position = "none")

# When considering N and C together, I need the dataset 
# that only includes plants/compartments that received nitrogen label.

mycoCforN = ggplot(data = nitrogeninfo) +
  geom_point(aes(x = mycologN15,
                 y = mycologC13, 
                 color = N_level,
                 shape = compartment_fungus)) +
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              aes(x = mycologN15,
                  y = mycologC13),
              color = "black",
              size = 0.5) +
  ylab(bquote(atop("Mycorrhizal "^13*"C", "(ln ppm excess)"))) +
  xlab(bquote(atop("Mycorrhizal "^15*"N (ln ppm excess)"))) +
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

save_plot("plots/Multipanel_regressions_myco_N_and_C.jpeg",
          threepanels,
          base_aspect_ratio = 3.6)

save_plot("plots/Multipanel_regressions_myco_N_and_C.pdf",
          threepanels,
          base_aspect_ratio = 3.6)

### Stats ###
# CforNlm = lm((mycoC13ppmexcess) ~ mycoN15ppmexcess, data = nitrogeninfo)
# plot(CforNlm) # seems okay
# summary(CforNlm)

# Use N15 values forced positive with linear transformation
# so you can try a log fit.
# CforNlm_log = lm((mycoC13ppmexcess) ~ mycologN15, data = nitrogeninfo)
# plot(CforNlm_log) # Kinda weird, esp point 25.
# summary(CforNlm_log) # Adj R^2 = 0.4287

CforNlm_loglog = lm(mycologC13 ~ mycologN15, data = nitrogeninfo)
plot(CforNlm_loglog) # Better
summary(CforNlm_loglog) # Adj R^2 = 0.4644


sink("stats_tables/myco13C_vs_myco_15N_lm_loglogresults.html")

stargazer(CforNlm_loglog, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()

# sink("stats_tables/myco13C_vs_myco_15N_comparingthreemodels.html")
# 
# stargazer(CforNlm, CforNlm_log, CforNlm_loglog, 
#           type = "html",
#           align = TRUE,
#           digits = 3,
#           star.cutoffs = c(0.05, 0.01, 0.001),
#           digit.separator = "",
#           no.space = TRUE)
# 
# sink()
