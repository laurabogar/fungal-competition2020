# Re-doing C for N figures
# 28 October 2020

library(cowplot)
library(tidyverse)
library(stargazer)
library(lme4)
library(lmerTest)
library(emmeans)

carboninfo = read_csv("processeddata/data_for_carbon_only_analyses_withpctC.csv.csv")
nitrogeninfo = read_csv("processeddata/isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates_withpctN.csv.csv")

carboninfo$hyphalog13C = log(carboninfo$hyphae.ppm13Cexcess)
nitrogeninfo$hyphae.ppm13Cexcess =  nitrogeninfo$hyphae.APE13C*10^4
nitrogeninfo$hyphae.ppm15Nexcess =  nitrogeninfo$hyphae.APE15N*10^4


nitrogeninfo_nooutlier = subset(nitrogeninfo, Plant != 6041)
carboninfo_nooutlier = carboninfo[!carboninfo$hyphae.APE13C == max(carboninfo$hyphae.APE13C),] # omit outlier 6024b

#### EXCLUDING EXTREME VALUE ####

### 1: Using all mycos, N15 in mycos predicts C13 in mycos ####

mycoCforN = ggplot(data = nitrogeninfo_nooutlier) +
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

lmer_model_mycoCforN = lmer(mycologC13 ~ mycologN15*compartment_fungus*N_level + (1|Batch),
                            data = nitrogeninfo_nooutlier)
summary(lmer_model_mycoCforN) #NS
anova(lmer_model_mycoCforN) # mycologN15 significant

### 2: Separating mycos out by species, Tt mycos get more C when nearby roots have 15N. Sp, not so. ####
# labels = c(Tt = "High N", Low = "Low N")

ByfungusCforNmyco = ggplot(data = nitrogeninfo_nooutlier) +
  geom_point(aes(x = mycologN15,
                 y = mycologC13, 
                 color = N_level)) +
  # geom_smooth(method = "lm", 
  #             formula = y ~ x, 
  #             data = subset(nitrogeninfo_nooutlier, compartment_fungus == "Tt"),
  #             aes(x = mycologN15,
  #                 y = mycologC13),
  #             color = "black",
  #             size = 0.5) +
  ylab(bquote(atop("Mycorrhizal "^13*"C", "(ln ppm excess)"))) +
  xlab(bquote(atop("Mycorrhizal "^15*"N (ln ppm excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  facet_grid(. ~ compartment_fungus) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

ByfungusCforNroots = ggplot(data = nitrogeninfo_nooutlier) +
  geom_point(aes(x = nmlogN15,
                 y = mycologC13, 
                 color = N_level)) +
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              aes(x = nmlogN15,
                  y = mycologC13),
              color = "black",
              size = 0.5) +
  ylab(bquote(atop("Mycorrhizal "^13*"C", "(ln ppm excess)"))) +
  xlab(bquote(atop("Root "^15*"N (ln ppm excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  facet_grid(. ~ compartment_fungus) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

ByfungusNforNroots = ggplot(data = nitrogeninfo_nooutlier) +
  geom_point(aes(x = mycologN15,
                 y = nmlogN15, 
                 color = N_level)) +
  # geom_smooth(method = "lm", 
  #             formula = y ~ x, 
  #             aes(x = mycologN15,
  #                 y = nmlogN15),
  #             color = "black",
  #             size = 0.5) +
  xlab(bquote(atop("Mycorrhizal "^15*"N", "(ln ppm excess)"))) +
  ylab(bquote(atop("Root "^15*"N (ln ppm excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  facet_grid(. ~ compartment_fungus) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

ByfungusmycoCforpctN = ggplot(data = nitrogeninfo_nooutlier) +
  geom_point(aes(y = mycologC13,
                 x = mycorrhizas.pctN, 
                 color = N_level)) +
  # geom_smooth(method = "lm", 
  #             formula = y ~ x, 
  #             aes(x = mycologN15,
  #                 y = nmlogN15),
  #             color = "black",
  #             size = 0.5) +
  ylab(bquote(atop("Mycorrhizal "^13*"C", "(ln ppm excess)"))) +
  xlab(bquote(atop("Mycorrhiza percent N"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  facet_grid(. ~ compartment_fungus) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

ByfungusnmCforpctN = ggplot(data = nitrogeninfo_nooutlier) +
  geom_point(aes(y = nmlogC13,
                 x = uncolonized_roots.pctN, 
                 color = N_level)) +
  # geom_smooth(method = "lm", 
  #             formula = y ~ x, 
  #             aes(x = mycologN15,
  #                 y = nmlogN15),
  #             color = "black",
  #             size = 0.5) +
  ylab(bquote(atop("NM "^13*"C", "(ln ppm excess)"))) +
  xlab(bquote(atop("NM percent N"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  facet_grid(. ~ compartment_fungus) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

ByfungusnmCfor15N = ggplot(data = nitrogeninfo_nooutlier) +
  geom_point(aes(y = nmlogC13,
                 x = nmlogN15, 
                 color = N_level)) +
  # geom_smooth(method = "lm", 
  #             formula = y ~ x, 
  #             aes(x = mycologN15,
  #                 y = nmlogN15),
  #             color = "black",
  #             size = 0.5) +
  ylab(bquote(atop("NM "^13*"C", "(ln ppm excess)"))) +
  xlab(bquote(atop("NM 15N (ln ppm excess"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  facet_grid(. ~ compartment_fungus) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

#### INCLUDING EXTREME VALUE ####