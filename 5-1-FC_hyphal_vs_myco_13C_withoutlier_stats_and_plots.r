# 5-1-FC_hyphal_vs_myco_13C_withoutlier_statsandplot

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(cowplot)
library(tidyverse)
library(lme4)
library(lmerTest)

together = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly.csv")


#### How well did hyphal C track myco C? ####
# For any analysis involving ONLY C, it makes sense to include
# 1) carbon data from root compartments that
# didn't necessarily receive N label (20 of these, 22 labeled compartments here with hyphal C info), and
# 2) data from plants with failed splits (two of these, 40 successfully split)

carboninfo = subset(together, compartment_fungus != "None" &
                      compartment_fungus != "MIXED" &
                      compartment_fungus != "OTHER")

carboninfo = carboninfo[!is.na(carboninfo$hyphae.APE13C),]
carboninfo = carboninfo[!is.na(carboninfo$mycorrhizas.APE13C),]
carboninfo$hyphae.ppm13Cexcessexcess = carboninfo$hyphae.APE13C*(10^4)

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