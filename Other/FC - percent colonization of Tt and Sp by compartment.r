# FC -- analyzing percent colonization of each fungus by compartment

library(tidyverse)
library(cowplot)
library(lme4)
library(lmerTest)

granular_data_bycompt = read_csv("processeddata/granular_mass_and_colonization_data_by_compartment.csv")
granular_data_byplant = read_csv("processeddata/granular_mass_and_colonization_data_by_plant.csv")
totalmass = granular_data_byplant %>%
  group_by(Plant) %>%
  summarize(total_biomass = sum(total_root_mass_wholesystem, shoot_biomass))

granular_data_byplant = left_join(granular_data_byplant, totalmass)

granular_data_byplant = subset(granular_data_byplant, N_level != "None")

granular_data_byplant$bettercomp = fct_relevel(granular_data_byplant$competitors_attempted,
                         "None/None",
                         "Sp/None",
                         "Sp/Sp",
                         "Tt/Sp",
                         "Tt/None",
                         "Tt/Tt")
## With facets
labels = c(High = "High N", Low = "Low N")

Tt_vs_Sp = ggplot(data = granular_data_byplant) +
  geom_point(aes(x = percent_Tt,
                 y = percent_Sp,
                 fill = bettercomp),
             position = position_jitter(width = 0.2),
             pch = 21,
             size = 2) +
  ylab("Percent Sp colonization at harvest") +
  xlab("Percent Tt colonization at harvest") +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  scale_fill_manual(values = c("white", "lightgoldenrod1", "darkgoldenrod1", "tan4", "plum2", "darkorchid4"),
                     name = "Fungi applied") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

save_plot("plots/percent_colonization_Tt_vs_Sp.pdf",
          Tt_vs_Sp,
          base_aspect_ratio = 2)

## No facets
Tt_vs_Sp = ggplot(data = granular_data_byplant) +
  geom_point(aes(x = percent_Tt,
                 y = percent_Sp,
                 color = N_level,
                 shape = bettercomp)) +
  ylab("Percent Sp at harvest") +
  xlab("Percent Tt  at harvest") +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(18, 2, 17, 22, 1, 19),
                     name = "Fungi applied") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

save_plot("plots/percent_colonization_Tt_vs_Sp.pdf",
          Tt_vs_Sp,
          base_aspect_ratio = 1.4)

### Experimental plot ####
granular_data_byplant$Tt_to_Sp = log((granular_data_byplant$percent_Tt+1)/(granular_data_byplant$percent_Sp +1))
granular_data_byplant$percent_col_overall = (granular_data_byplant$Tt_myco_mass_wholesystem + granular_data_byplant$Sp_myco_mass_wholesystem)/granular_data_byplant$percent_col_denominator
relative_dominance_colonization = ggplot(data = granular_data_byplant) +
  geom_jitter(aes(x = Tt_to_Sp,
                 y = percent_col_overall,
                 fill = bettercomp),
             # position = position_jitter(width = 0.2),
             pch = 21,
             size = 2) +
  ylab("Percent colonization at harvest") +
  xlab("Tt to Sp log ratio") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  scale_fill_manual(values = c("white", "lightgoldenrod1", "darkgoldenrod1", "tan4", "plum2", "darkorchid4"),
                    name = "Fungi applied") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

save_plot("plots/percent_colonization_Tt_vs_Sp_relative_dominance.pdf",
          relative_dominance_colonization,
          base_aspect_ratio = 2)

relative_dominance_biomass = ggplot(data = granular_data_byplant) +
  geom_jitter(aes(x = Tt_to_Sp,
                  y = total_biomass,
                  fill = bettercomp),
              # position = position_jitter(width = 0.2),
              pch = 21,
              size = 2) +
  ylab("Total biomass at harvest") +
  xlab("Tt to Sp log ratio") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  scale_fill_manual(values = c("white", "lightgoldenrod1", "darkgoldenrod1", "tan4", "plum2", "darkorchid4"),
                    name = "Fungi applied") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

save_plot("plots/total_biomass_Tt_vs_Sp_relative_dominance.pdf",
          relative_dominance_biomass,
          base_aspect_ratio = 2)

granular_data_byplant$Tt_to_Sp_subtract = granular_data_byplant$percent_Tt - granular_data_byplant$percent_Sp
relative_dominance_colonization_linear = ggplot(data = granular_data_byplant) +
  geom_jitter(aes(x = Tt_to_Sp_subtract,
                  y = percent_col_overall,
                  fill = bettercomp),
              # position = position_jitter(width = 0.2),
              pch = 21,
              size = 2) +
  ylab("Percent colonization at harvest") +
  xlab("Percent Tt minus percent Sp") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  scale_fill_manual(values = c("white", "lightgoldenrod1", "darkgoldenrod1", "tan4", "plum2", "darkorchid4"),
                    name = "Fungi applied") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

save_plot("plots/percent_colonization_Tt_vs_Sp_relative_dominance_linearx.pdf",
          relative_dominance_colonization_linear,
          base_aspect_ratio = 2)

relative_dominance_biomass_linear = ggplot(data = granular_data_byplant) +
  geom_jitter(aes(x = Tt_to_Sp_subtract,
                  y = total_biomass,
                  fill = bettercomp),
              # position = position_jitter(width = 0.2),
              pch = 21,
              size = 2) +
  ylab("Total biomass at harvest") +
  xlab("Percent Tt minus percent Sp") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  scale_fill_manual(values = c("white", "lightgoldenrod1", "darkgoldenrod1", "tan4", "plum2", "darkorchid4"),
                    name = "Fungi applied") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

save_plot("plots/total_biomass_Tt_vs_Sp_relative_dominance_linear.pdf",
          relative_dominance_biomass_linear,
          base_aspect_ratio = 2)
