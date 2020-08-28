# Looking into mixed compartments

# I've just realized that I have (limited) evidence for complimentarity
# between Sp and Tt: The plant got more N in mixed colonization compartments.

# What happens with the carbon?


setwd("~/Documents/Fungal competition project/fungal-competition2020/")

# Libraries needed:
library(cowplot)
library(tidyverse)
library(lme4)
library(lmerTest)
library(stargazer)

mixeddata = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly_INCLUDING_MIXED.csv")

summary(as.factor(mixeddata$compartment_fungus)) # I have 20 mixed compartments

ndata = subset(mixeddata, received15N == "Y")
ndata = subset(ndata, compartment_fungus != "None")

summary(as.factor(ndata$compartment_fungus)) # 10 mixed with N label. Is this really true?

# In mixed compartments, did Sp get more C than Tt?

onlymixed = subset(mixeddata, compartment_fungus == "MIXED")
# Where are 6056a NM roots?


labels = c(High = "High N", Low = "Low N")
carbon_mycos = ggplot(data = onlymixed) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = mycofungus, 
                   y = mycologC13)) +
  geom_jitter(width = 0.20,
              aes(x = mycofungus, 
                  y = mycologC13,
                  shape = mycofungus,
                  color = as.factor(Plant))) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Labeled C in mycorrhizas\n(ln ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus")

nitrogen_mycos = ggplot(data = onlymixed) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = mycofungus, 
                   y = mycologN15)) +
  geom_jitter(width = 0.20,
              aes(x = mycofungus, 
                  y = mycologN15,
                  shape = mycofungus,
                  color = as.factor(Plant))) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Labeled N in mycorrhizas\n(ln ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus")
