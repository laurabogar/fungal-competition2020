# Re-doing C for N figures
# 28 October 2020

library(cowplot)
library(tidyverse)
library(stargazer)
library(lme4)
library(lmerTest)
library(emmeans)
library(MuMIn)

carboninfo = read_csv("processeddata/data_for_carbon_only_analyses_withpctC.csv.csv")
nitrogeninfo = read_csv("processeddata/isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates_withpctN.csv.csv")

carboninfo$hyphalog13C = log(carboninfo$hyphae.ppm13Cexcess)
nitrogeninfo$hyphae.ppm13Cexcess =  nitrogeninfo$hyphae.APE13C*10^4
nitrogeninfo$hyphae.ppm15Nexcess =  nitrogeninfo$hyphae.APE15N*10^4


nitrogeninfo_nooutlier = subset(nitrogeninfo, Plant != 6041)
carboninfo_nooutlier = carboninfo[!carboninfo$hyphae.APE13C == max(carboninfo$hyphae.APE13C),] # omit outlier 6024b

#### EXCLUDING EXTREME VALUE ####

### 1:  N15 in mycos predicts C13 in mycos, but only for Thelephora ####
ByfungusCforNmyco = ggplot(data = nitrogeninfo_nooutlier) +
  geom_point(aes(x = mycologN15,
                 y = mycologC13, 
                 color = N_level)) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              data = subset(nitrogeninfo_nooutlier, compartment_fungus == "Tt"),
              aes(x = mycologN15,
                  y = mycologC13),
              color = "black",
              size = 0.5) +
  ylab(bquote(atop("Mycorrhizal "^13*"C (ln ppm excess)"))) +
  xlab(bquote(atop("Mycorrhizal "^15*"N (ln ppm excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  facet_grid(. ~ compartment_fungus) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

mycoCforN_justTt = lmer(mycologC13 ~ mycologN15*N_level + (1|Batch),
                            data = subset(nitrogeninfo_nooutlier, compartment_fungus == "Tt"))
mycoCforN_justTt_table = summary(mycoCforN_justTt)
class(mycoCforN_justTt) <- "lmerMod"
r.squaredGLMM(mycoCforN_justTt)

sink("stats_tables/mycoCforN_justTt_lmer.html")

stargazer(mycoCforN_justTt, type = "html",
          digits = 3,
          digit.separator = "",
          summary = TRUE)

sink()

mycoCforN_justSp = lmer(mycologC13 ~ mycologN15*N_level + (1|Batch),
                        data = subset(nitrogeninfo_nooutlier, compartment_fungus == "Sp"))
mycoCforN_justSp_table = summary(mycoCforN_justSp)
class(mycoCforN_justSp) <- "lmerMod"
r.squaredGLMM(mycoCforN_justSp)


sink("stats_tables/mycoCforN_justSp_lmer.html")

stargazer(mycoCforN_justSp, type = "html",
          digits = 3,
          digit.separator = "",
          summary = TRUE)

sink()

sink("stats_tables/mycoCforN_TtvsSp_lmer.html")

stargazer(mycoCforN_justTt, mycoCforN_justSp, type = "html",
          digits = 3,
          digit.separator = "",
          summary = TRUE)

sink()

### 2: Fungal N goes to NM roots, but only for Tt ####

ByfungusNforNroots = ggplot(data = nitrogeninfo_nooutlier) +
  geom_point(aes(x = mycologN15,
                 y = nmlogN15, 
                 color = N_level)) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              data = subset(nitrogeninfo_nooutlier, compartment_fungus == "Tt"),
              aes(x = mycologN15,
                  y = nmlogN15),
              color = "black",
              size = 0.5) +
  xlab(bquote(atop("Mycorrhizal "^15*"N (ln ppm excess)"))) +
  ylab(bquote(atop("Root "^15*"N (ln ppm excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  facet_grid(. ~ compartment_fungus) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

NforN_justTt = lmer(nmlogN15 ~ mycologN15*N_level + (1|Batch),
                        data = subset(nitrogeninfo_nooutlier, compartment_fungus == "Tt"))
class(NforN_justTt) <- "lmerMod"
r.squaredGLMM(NforN_justTt)

sink("stats_tables/NforN_justTt_lmer.html")

stargazer(NforN_justTt, type = "html",
          digits = 3,
          digit.separator = "",
          summary = TRUE)

sink()

NforN_justSp = lmer(nmlogN15 ~ mycologN15*N_level + (1|Batch),
                        data = subset(nitrogeninfo_nooutlier, compartment_fungus == "Sp"))
class(NforN_justSp) <- "lmerMod"
r.squaredGLMM(NforN_justSp) #very bad


sink("stats_tables/NforN_justSp_lmer.html")

stargazer(NforN_justSp, type = "html",
          digits = 3,
          digit.separator = "",
          summary = TRUE)

sink()

sink("stats_tables/NforN_TtvsSp_lmer.html")

stargazer(NforN_justTt, NforN_justSp, type = "html",
          digits = 3,
          digit.separator = "",
          summary = TRUE)

sink()

### 3: Plant C goes to hyphae in Tt ####

hyphalCformycoC_plot_nooutlier = ggplot(data = subset(carboninfo_nooutlier, compartment_fungus == "Tt")) +
  geom_point(aes(x = mycologC13,
                 y = log(hyphae.ppm13Cexcess), 
                 color = N_level)) +
  geom_smooth(method = "lm", 
              aes(x = mycologC13,
                  y = log(hyphae.ppm13Cexcess)),
              color = "black",
              size = 0.5) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  ylab(expression("Hyphal "^13*"C in Tt (ln ppm excess)")) +
  xlab(expression("Mycorrhizal "^13*"C in Tt (ln ppm excess)")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) 

CforC_justTt = lmer(hyphalog13C ~ mycologC13*N_level + (1|Batch),
                    data = subset(carboninfo_nooutlier, compartment_fungus == "Tt"))
class(CforC_justTt) <- "lmerMod"
summary(CforC_justTt)
r.squaredGLMM(CforC_justTt)


sink("stats_tables/CforC_justTt_lmer.html")

stargazer(CforC_justTt, type = "html",
          digits = 3,
          digit.separator = "",
          summary = TRUE)

sink()

### Bringing together for three panel figure ####
ByfungusCforNmyco_nolegend = ByfungusCforNmyco +
  theme(legend.position = "none")

ByfungusNforNroots_nolegend = ByfungusNforNroots +
  theme(legend.position = "none")


threepanels_horiz = plot_grid(ByfungusCforNmyco_nolegend, 
                              ByfungusNforNroots_nolegend,
                              hyphalCformycoC_plot_nooutlier,
                              labels = c("A", "B", "C"),
                              align = "h",
                              nrow = 1,
                              ncol = 3,
                              rel_widths = c(1, 1, .9))


save_plot("plots/horizontal_three_panel_plot_CforN_relationships_fig1.pdf", 
          threepanels_horiz,
          base_width = 15)



#### extra stuff ####

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

# lumping Tt and Sp for the sake of illustration
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
#### INCLUDING EXTREME VALUE ####