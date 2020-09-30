# FC -- analyzing percent colonization of each fungus by compartment

library(tidyverse)
library(cowplot)
library(lme4)
library(lmerTest)

granular_data = read_csv("processeddata/granular_mass_and_colonization_data_by_compartment.csv")

ggplot(data = granular_data) +
  geom_point(aes(x = percent_Tt_mycos,
                 y = percent_Sp_mycos,
                 color = N_level,
                 shape = compartment_fungus_attempted))

granular_data$compartment_fungus == "Mixed"
