#5-4-FC_hyphal_vs_myco_13C_without_outlier_stats

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(cowplot)
library(tidyverse)
library(lme4)
library(lmerTest)

carboninfo = read_csv("processeddata/data_for_carbon_only_analyses.csv")

#### How well did hyphal C track myco C? ####

hyphalCformycoC_plot_withoutlier = ggplot(data = carboninfo) +
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

save_plot("plots/Regression_hyphal_C_for_myco_C_with_outlier.pdf", 
          hyphalCformycoC_plot_withoutlier,
          ncol = 1,
          base_aspect_ratio = 1.4)

ggsave("plots/Regression_hyphal_C_for_myco_C_with_outlier.jpeg", 
       plot = hyphalCformycoC_plot_withoutlier,
       device = "jpeg",
       width = 7, height = 5, units = "in")

### Stats with outlier ###
myhyphaelm_withoutlier = lm((hyphae.ppm13Cexcess) ~ mycoC13ppmexcess, data = carboninfo)
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

sink("stats_tables/hyphae13C_vs_myco13C_withoutlier_lm_results.html")

stargazer(myhyphaelm_withoutlier, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE)

sink()