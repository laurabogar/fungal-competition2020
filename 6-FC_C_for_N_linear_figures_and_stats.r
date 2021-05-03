# Re-doing C for N figures
# 28 October 2020

library(cowplot)
library(tidyverse)
library(stargazer)
library(lme4)
library(lmerTest)
library(emmeans)
library(MuMIn)

carboninfo = read_csv("processeddata/data_for_carbon_only_analyses_withpctC.csv") # from script 4-alternative
nitrogeninfo = read_csv("processeddata/isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates_withpctN.csv")


carboninfo$hyphalog13C = log(carboninfo$hyphae.ppm13Cexcess)
nitrogeninfo$hyphae.ppm13Cexcess =  nitrogeninfo$hyphae.APE13C*10^4
nitrogeninfo$hyphae.ppm15Nexcess =  nitrogeninfo$hyphae.APE15N*10^4

nocontam = read_csv("processeddata/Plants_with_no_Tt_contamination.csv")
nitrogeninfo = mutate(nitrogeninfo, clean = nitrogeninfo$Plant %in% nocontam$Plant)
carboninfo = mutate(carboninfo, clean = carboninfo$Plant %in% nocontam$Plant)


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

# Version for talk:
ByfungusCforNmyco_purple = ggplot(data = nitrogeninfo_nooutlier) +
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
  scale_color_manual(values = c("mediumpurple", "thistle3"),
                     name = "N level") +
  facet_grid(. ~ compartment_fungus) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

save_plot("plots/Myco_C_for_N_purple_for_talk.pdf",
          ByfungusCforNmyco_purple,
          base_aspect_ratio = 1.7)

mycoCforN_justTt = lmer(mycologC13 ~ mycologN15*N_level + (1|Batch),
                            data = subset(nitrogeninfo_nooutlier, compartment_fungus == "Tt"))
mycoCforN_justTt_table = summary(mycoCforN_justTt)
class(mycoCforN_justTt) <- "lmerMod"
r.squaredGLMM(mycoCforN_justTt)
fitTt = as.data.frame(r.squaredGLMM(mycoCforN_justTt))
fitTt_R2 = round(fitTt$R2c, 4)

# Does it change if we only use uncontaminated plants?

mycoCforN_justTt_onlyclean = lmer(mycologC13 ~ mycologN15*N_level + (1|Batch),
                        data = subset(nitrogeninfo_nooutlier, compartment_fungus == "Tt" & clean == TRUE))
mycoCforN_justTt_table = summary(mycoCforN_justTt_onlyclean)
class(mycoCforN_justTt_onlyclean) <- "lmerMod"
r.squaredGLMM(mycoCforN_justTt_onlyclean)
fitTt = as.data.frame(r.squaredGLMM(mycoCforN_justTt_onlyclean))
fitTt_R2 = round(fitTt$R2c, 4) #fit is slightly better with just clean plants, patterns the same

# This allows me to visualize the partial
# residuals of my model -- basically subtracting
# out the influence of N level and its interaction with
# mycorrhiza 15N concentration, I think.
# Ultimately, though, this plot looks very similar
# to the raw data, and is a bit more convoluted to interpret.
library(effects)
est = Effect("mycologN15", partial.residuals = TRUE, mycoCforN_justTt)
plot(est)

# Found remef package hard to use
# library(remef)
# test = remef(mycoCforN_justTt, fix = 2, ran = "all")
# plot(test)

sink("stats_tables/mycoCforN_justTt_lmer.html")

stargazer(mycoCforN_justTt, type = "html",
          digits = 3,
          digit.separator = "",
          ci = TRUE,
          summary = TRUE)

sink()

mycoCforN_justSp = lmer(mycologC13 ~ mycologN15*N_level + (1|Batch),
                        data = subset(nitrogeninfo_nooutlier, compartment_fungus == "Sp"))
mycoCforN_justSp_table = summary(mycoCforN_justSp)
class(mycoCforN_justSp) <- "lmerMod"
r.squaredGLMM(mycoCforN_justSp)
fitSp = as.data.frame(r.squaredGLMM(mycoCforN_justSp))
fitSp_R2 = round(fitSp$R2c, 4)

# Does it look the same if we eliminate plants with any Tt contamination?

mycoCforN_justSp_onlyclean = lmer(mycologC13 ~ mycologN15*N_level + (1|Batch),
                        data = subset(nitrogeninfo_nooutlier, compartment_fungus == "Sp" & clean == TRUE))
mycoCforN_justSp_onlyclean_table = summary(mycoCforN_justSp_onlyclean)
class(mycoCforN_justSp_onlyclean) <- "lmerMod"
r.squaredGLMM(mycoCforN_justSp_onlyclean)
fitSp = as.data.frame(r.squaredGLMM(mycoCforN_justSp_onlyclean))
fitSp_R2 = round(fitSp$R2c, 4)
# Fit gets a little better, but patterns are the same.

sink("stats_tables/mycoCforN_justSp_lmer.html")

stargazer(mycoCforN_justSp, type = "html",
          digits = 3,
          digit.separator = "",
          ci = TRUE,
          summary = TRUE)

sink()

sink("stats_tables/mycoCforN_TtvsSp_lmer.html")

stargazer(mycoCforN_justTt, 
          mycoCforN_justSp, 
          type = "html",
          dep.var.labels = "Mycorrhizal [13C] (ln ppm excess)",
          covariate.labels = c("Mycorrhizal [15N] (ln ppm excess)",
                               "N level",
                               "Mycorrhizal [15N]:N level"),
          column.labels = c("Tt", "Sp"),
          digits = 3,
          digit.separator = "",
          ci = TRUE,
          star.cutoffs = c(0.05, 0.01, 0.001),
          add.lines = list(c("Conditional pseudo-$R2$",
                             fitTt_R2, fitSp_R2)),
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

# with only uncontaminated plants?
NforN_justTt_onlyclean = lmer(nmlogN15 ~ mycologN15*N_level + (1|Batch),
                    data = subset(nitrogeninfo_nooutlier, compartment_fungus == "Tt" & clean == TRUE))
class(NforN_justTt_onlyclean) <- "lmerMod"
r.squaredGLMM(NforN_justTt_onlyclean)
# fit gets a little worse, but pattern is the same

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

stargazer(NforN_justSp, 
          type = "html",
          digits = 3,
          digit.separator = "",
          summary = TRUE)

sink()


fitTt = as.data.frame(r.squaredGLMM(NforN_justTt))
fitTt_R2 = round(fitTt$R2c, 4)
fitSp = as.data.frame(r.squaredGLMM(NforN_justSp))
fitSp_R2 = round(fitSp$R2c, 4)


sink("stats_tables/NforN_TtvsSp_lmer.html")

stargazer(NforN_justTt, 
          NforN_justSp, 
          type = "html",
          dep.var.labels = "Uncolonized root [15N] (ln ppm excess)",
          covariate.labels = c("Mycorrhizal [15N] (ln ppm excess)",
                               "N level",
                               "Mycorrhizal [15N]:N level"),
          column.labels = c("Tt", "Sp"),
          digits = 3,
          digit.separator = "",
          ci = TRUE,
          star.cutoffs = c(0.05, 0.01, 0.001),
          add.lines = list(c("Conditional pseudo-$R2$",
                             fitTt_R2, fitSp_R2)),
          summary = TRUE)

sink()

# ### 3: Plant C goes to hyphae in Tt ####
# 
# hyphalCformycoC_plot = ggplot(data = subset(carboninfo, compartment_fungus == "Tt")) +
#   geom_point(aes(x = mycologC13,
#                  y = log(hyphae.ppm13Cexcess), 
#                  color = N_level)) +
#   geom_smooth(method = "lm", 
#               aes(x = mycologC13,
#                   y = log(hyphae.ppm13Cexcess)),
#               color = "black",
#               size = 0.5) +
#   scale_color_manual(values = c("steelblue4", "steelblue1"),
#                      name = "N level") +
#   scale_shape_manual(values = c(17, 15),
#                      name = "Fungus") +
#   ylab(expression("Hyphal "^13*"C in Tt (ln ppm excess)")) +
#   xlab(expression("Mycorrhizal "^13*"C in Tt (ln ppm excess)")) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm")) 
# test = lm(hyphalog13C ~ mycologC13*N_level, data = carboninfo_nooutlier)
# plot(test)
# CforC_justTt = lmer(hyphalog13C ~ mycologC13*N_level + (1|Batch),
#                     data = subset(carboninfo, compartment_fungus == "Tt"))
# class(CforC_justTt) <- "lmerMod"
# summary(CforC_justTt)
# r.squaredGLMM(CforC_justTt)
# 
# fitTt = as.data.frame(r.squaredGLMM(CforC_justTt))
# fitTt_R2 = round(fitTt$R2c, 4)
# 
# sink("stats_tables/CforC_justTt_lmer.html")
# 
# stargazer(CforC_justTt, type = "html",
#           dep.var.labels = "Hyphal [13C] (ln ppm excess)",
#           covariate.labels = c("Mycorrhizal [13C] (ln ppm excess)",
#                                "N level",
#                                "Mycorrhizal [13C]:N level"),
#           column.labels = c("Tt", "Sp"),
#           digits = 3,
#           digit.separator = "",
#           ci = TRUE,
#           star.cutoffs = c(0.05, 0.01, 0.001),
#           add.lines = list(c("Conditional pseudo-$R2$",
#                              fitTt_R2)),
#           summary = TRUE)
# 
# sink()

# ### Bringing together for three panel figure ####
# ByfungusCforNmyco_nolegend = ByfungusCforNmyco +
#   theme(legend.position = "none")
# 
# ByfungusNforNroots_nolegend = ByfungusNforNroots +
#   theme(legend.position = "none")
# 
# hyphalCformycoC_plot_legendbelow = hyphalCformycoC_plot_nooutlier +
#   theme(legend.position = "bottom")
# 
# threepanels_horiz = plot_grid(ByfungusCforNmyco_nolegend, 
#                               ByfungusNforNroots_nolegend,
#                               hyphalCformycoC_plot,
#                               labels = c("a", "b", "c"),
#                               align = "h",
#                               axis = "b",
#                               nrow = 1,
#                               ncol = 3,
#                               rel_widths = c(1, 1, .9))
# 
# 
# save_plot("plots/horizontal_three_panel_plot_CforN_relationships_fig1.pdf", 
#           threepanels_horiz,
#           base_width = 15)
# 
# threepanels_vertical = plot_grid(ByfungusCforNmyco_nolegend, 
#                               ByfungusNforNroots_nolegend,
#                               hyphalCformycoC_plot_legendbelow,
#                               labels = c("a", "b", "c"),
#                               align = "v",
#                               axis = "l",
#                               nrow = 3,
#                               ncol = 1,
#                               rel_heights = c(1,1,1))
# 
# save_plot("plots/vertical_three_panel_plot_CforN_relationships_fig1.pdf", 
#           threepanels_vertical,
#           base_height = 12,
#           base_width = 5)



### Two panels excluding extreme value ####

ByfungusCforNmyco_nolegend = ByfungusCforNmyco +
  theme(legend.position = "none")

ByfungusNforNroots_legendbelow = ByfungusNforNroots +
  theme(legend.position = "bottom")

twopanels_vertical_noextreme = plot_grid(ByfungusCforNmyco_nolegend, 
                               ByfungusNforNroots_legendbelow,
                               labels = c("a", "b"),
                               align = "v",
                               axis = "l",
                               nrow = 2,
                               ncol = 1,
                               rel_heights = c(1,1))

save_plot("plots/vertical_TWO_panel_plot_CforN_relationships_noextreme_forsupp.pdf", 
          twopanels_vertical_noextreme,
          base_height = 8,
          base_width = 5)

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

nitrogen_nocontam = subset(nitrogeninfo, Plant %in% nocontam$Plant)
# if you do this, you lose 21 samples of 59 total. This is definitely not worth it.

### 1:  N15 in mycos predicts C13 in mycos, but only for Thelephora ####
ByfungusCforNmyco_ie = ggplot(data = nitrogeninfo) +
  geom_point(aes(x = mycologN15,
                 y = mycologC13, 
                 color = N_level)) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              data = subset(nitrogeninfo, compartment_fungus == "Tt"),
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

mycoCforN_justTt_includingextreme = lmer(mycologC13 ~ mycologN15*N_level + (1|Batch),
                        data = subset(nitrogeninfo, compartment_fungus == "Tt"))
mycoCforN_justTt_table_ie = summary(mycoCforN_justTt_includingextreme)
class(mycoCforN_justTt_includingextreme) <- "lmerMod"
fitTt = as.data.frame(r.squaredGLMM(mycoCforN_justTt_includingextreme))
fitTt_R2 = round(fitTt$R2c, 4)


mycoCforN_justSp_ie = lmer(mycologC13 ~ mycologN15*N_level + (1|Batch),
                        data = subset(nitrogeninfo, compartment_fungus == "Sp"))
mycoCforN_justSp_table_ie = summary(mycoCforN_justSp_ie)
class(mycoCforN_justSp_ie) <- "lmerMod"
fitSp = as.data.frame(r.squaredGLMM(mycoCforN_justSp_ie))
fitSp_R2 = round(fitSp$R2c, 4)

# If we use only uncontaminated plants?
mycoCforN_justTt_includingextreme_onlyclean = lmer(mycologC13 ~ mycologN15*N_level + (1|Batch),
                                         data = subset(nitrogeninfo, compartment_fungus == "Tt" & clean == TRUE))
mycoCforN_justTt_table_ie_onlyclean = summary(mycoCforN_justTt_includingextreme_onlyclean)
class(mycoCforN_justTt_includingextreme_onlyclean) <- "lmerMod"
fitTt_onlyclean = as.data.frame(r.squaredGLMM(mycoCforN_justTt_includingextreme_onlyclean))
fitTt_R2_onlyclean = round(fitTt_onlyclean$R2c, 4)
# Fit for Tt is about the same, significance patterns & relationships also the same

mycoCforN_justSp_ie_onlyclean = lmer(mycologC13 ~ mycologN15*N_level + (1|Batch),
                           data = subset(nitrogeninfo, compartment_fungus == "Sp" & clean == TRUE))
mycoCforN_justSp_table_ie_onlyclean = summary(mycoCforN_justSp_ie_onlyclean)
class(mycoCforN_justSp_ie_onlyclean) <- "lmerMod"
fitSp_onlyclean = as.data.frame(r.squaredGLMM(mycoCforN_justSp_ie_onlyclean))
fitSp_R2_onlyclean = round(fitSp_onlyclean$R2c, 4)
# Improves fit a bit but patterns the same for Sp.

# sink("stats_tables/mycoCforN_TtvsSp_lmer_includingextreme.html")

models_CforNmycos_includingextreme = 
  stargazer(mycoCforN_justTt_includingextreme, 
          mycoCforN_justSp_ie, 
          type = "html",
          dep.var.labels = "Mycorrhizal [13C] (ln ppm excess)",
          covariate.labels = c("Mycorrhizal [15N] (ln ppm excess)",
                               "N level",
                               "Mycorrhizal [15N]:N level"),
          column.labels = c("Tt", "Sp"),
          digits = 3,
          digit.separator = "",
          ci = TRUE,
          star.cutoffs = c(0.05, 0.01, 0.001),
          add.lines = list(c("Conditional pseudo-$R2$",
                             fitTt_R2, fitSp_R2)),
          summary = TRUE)

cat(models_CforNmycos_includingextreme, 
    sep = '\n', 
    file = "stats_tables/mycoCforN_TtvsSp_lmer_includingextreme.html")

# no contam
models_CforNmycos_includingextreme_onlyclean = 
  stargazer(mycoCforN_justTt_includingextreme_onlyclean, 
            mycoCforN_justSp_ie_onlyclean, 
            type = "html",
            dep.var.labels = "Mycorrhizal [13C] (ln ppm excess)",
            covariate.labels = c("Mycorrhizal [15N] (ln ppm excess)",
                                 "N level",
                                 "Mycorrhizal [15N]:N level"),
            column.labels = c("Tt", "Sp"),
            digits = 3,
            digit.separator = "",
            ci = TRUE,
            star.cutoffs = c(0.05, 0.01, 0.001),
            add.lines = list(c("Conditional pseudo-$R2$",
                               fitTt_R2_onlyclean, fitSp_R2_onlyclean)),
            summary = TRUE)

cat(models_CforNmycos_includingextreme_onlyclean, 
    sep = '\n', 
    file = "stats_tables/mycoCforN_TtvsSp_lmer_includingextreme_onlyclean.html")
# sink()

### 2: Fungal N goes to NM roots, but only for Tt ####

ByfungusNforNroots_ie = ggplot(data = nitrogeninfo) +
  geom_point(aes(x = mycologN15,
                 y = nmlogN15, 
                 color = N_level)) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              data = subset(nitrogeninfo, compartment_fungus == "Tt"),
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

NforN_justTt_ie = lmer(nmlogN15 ~ mycologN15*N_level + (1|Batch),
                    data = subset(nitrogeninfo, compartment_fungus == "Tt"))
summary(NforN_justTt_ie)
class(NforN_justTt_ie) <- "lmerMod"
r.squaredGLMM(NforN_justTt_ie)
fitTt = as.data.frame(r.squaredGLMM(NforN_justTt_ie))
fitTt_R2 = round(fitTt$R2c, 4)


NforN_justSp_ie = lmer(nmlogN15 ~ mycologN15*N_level + (1|Batch),
                    data = subset(nitrogeninfo, compartment_fungus == "Sp"))
summary(NforN_justSp_ie)
class(NforN_justSp_ie) <- "lmerMod"
r.squaredGLMM(NforN_justSp_ie)
fitSp = as.data.frame(r.squaredGLMM(NforN_justSp_ie))
fitSp_R2 = round(fitSp$R2c, 4)


sink("stats_tables/NforN_TtvsSp_lmer_includingextreme.html")

stargazer(NforN_justTt_ie, 
          NforN_justSp_ie, 
          type = "html",
          dep.var.labels = "Uncolonized root [15N] (ln ppm excess)",
          covariate.labels = c("Mycorrhizal [15N] (ln ppm excess)",
                               "N level",
                               "Mycorrhizal [15N]:N level"),
          column.labels = c("Tt", "Sp"),
          digits = 3,
          digit.separator = "",
          ci = TRUE,
          star.cutoffs = c(0.05, 0.01, 0.001),
          add.lines = list(c("Conditional pseudo-$R2$",
                             fitTt_R2, fitSp_R2)),
          summary = TRUE)

sink()

# What if I only include uncontaminated plants?
NforN_justTt_ie_onlyclean = lmer(nmlogN15 ~ mycologN15*N_level + (1|Batch),
                                 data = subset(nitrogeninfo, compartment_fungus == "Tt" & clean == TRUE))
summary(NforN_justTt_ie_onlyclean)
class(NforN_justTt_ie_onlyclean) <- "lmerMod"
r.squaredGLMM(NforN_justTt_ie_onlyclean)
fitTt_onlyclean = as.data.frame(r.squaredGLMM(NforN_justTt_ie_onlyclean))
fitTt_R2 = round(fitTt_onlyclean$R2c, 4)
# better fit here, only myco N15 is a significant predictor anymore (no marginal effects)


NforN_justSp_ie_onlyclean = lmer(nmlogN15 ~ mycologN15*N_level + (1|Batch),
                                 data = subset(nitrogeninfo, compartment_fungus == "Sp" & clean == TRUE))
summary(NforN_justSp_ie_onlyclean)
class(NforN_justSp_ie_onlyclean) <- "lmerMod"
r.squaredGLMM(NforN_justSp_ie_onlyclean)
fitSp_onlyclean = as.data.frame(r.squaredGLMM(NforN_justSp_ie_onlyclean))
fitSp_R2_onlyclean = round(fitSp_onlyclean$R2c, 4)
# this model is equally worthless to the first one. (rsq is 0.05)

# only clean
sink("stats_tables/NforN_TtvsSp_lmer_includingextreme_onlyclean.html")

stargazer(NforN_justTt_ie_onlyclean, 
          NforN_justSp_ie_onlyclean, 
          type = "html",
          dep.var.labels = "Uncolonized root [15N] (ln ppm excess)",
          covariate.labels = c("Mycorrhizal [15N] (ln ppm excess)",
                               "N level",
                               "Mycorrhizal [15N]:N level"),
          column.labels = c("Tt", "Sp"),
          digits = 3,
          digit.separator = "",
          ci = TRUE,
          star.cutoffs = c(0.05, 0.01, 0.001),
          add.lines = list(c("Conditional pseudo-$R2$",
                             fitTt_R2_onlyclean, fitSp_R2_onlyclean)),
          summary = TRUE)

sink()

### 3: Plant C goes to hyphae in Tt ####

hyphalCformycoC_plot = ggplot(data = subset(carboninfo, compartment_fungus == "Tt")) +
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


test = lm(hyphalog13C ~ mycologC13*N_level, data = carboninfo_nooutlier)
plot(test)


CforC_justTt = lmer(hyphalog13C ~ mycologC13*N_level + (1|Batch),
                    data = subset(carboninfo, compartment_fungus == "Tt"))
summary(CforC_justTt)
class(CforC_justTt) <- "lmerMod"
r.squaredGLMM(CforC_justTt)

fitTt = as.data.frame(r.squaredGLMM(CforC_justTt))
fitTt_R2 = round(fitTt$R2c, 4)



sink("stats_tables/CforC_justTt_lmer2.html")

stargazer(CforC_justTt, type = "html",
          dep.var.labels = "Hyphal [13C] (ln ppm excess)",
          covariate.labels = c("Mycorrhizal [13C] (ln ppm excess)",
                               "N level",
                               "Mycorrhizal [13C]:N level"),
          column.labels = c("Tt", "Sp"),
          digits = 3,
          digit.separator = "",
          ci = TRUE,
          star.cutoffs = c(0.05, 0.01, 0.001),
          add.lines = list(c("Conditional pseudo-$R2$",
                             fitTt_R2)),
          summary = TRUE)

sink()

# what if we use only uncontaminated microcosms?

CforC_justTt_onlyclean = lmer(hyphalog13C ~ mycologC13*N_level + (1|Batch),
                              data = subset(carboninfo, compartment_fungus == "Tt" & clean == TRUE))
summary(CforC_justTt_onlyclean)
class(CforC_justTt_onlyclean) <- "lmerMod"
r.squaredGLMM(CforC_justTt_onlyclean)
# slightly better fit, same overall pattern.

fitTt_onlyclean = as.data.frame(r.squaredGLMM(CforC_justTt_onlyclean))
fitTt_R2_onlyclean = round(fitTt$R2c, 4)

sink("stats_tables/CforC_justTt_lmer_onlyclean.html")

stargazer(CforC_justTt_onlyclean, type = "html",
          dep.var.labels = "Hyphal [13C] (ln ppm excess)",
          covariate.labels = c("Mycorrhizal [13C] (ln ppm excess)",
                               "N level",
                               "Mycorrhizal [13C]:N level"),
          column.labels = c("Tt", "Sp"),
          digits = 3,
          digit.separator = "",
          ci = TRUE,
          star.cutoffs = c(0.05, 0.01, 0.001),
          add.lines = list(c("Conditional pseudo-$R2$",
                             fitTt_R2_onlyclean)),
          summary = TRUE)

sink()

### Bringing together for three panel figure ####
ByfungusCforNmyco_nolegend = ByfungusCforNmyco_ie +
  theme(legend.position = "none")

ByfungusNforNroots_nolegend = ByfungusNforNroots_ie +
  theme(legend.position = "none")

hyphalCformycoC_plot_legendbelow = hyphalCformycoC_plot +
  theme(legend.position = "bottom")

threepanels_horiz = plot_grid(ByfungusCforNmyco_nolegend, 
                              ByfungusNforNroots_nolegend,
                              hyphalCformycoC_plot,
                              labels = c("a", "b", "c"),
                              align = "h",
                              axis = "b",
                              nrow = 1,
                              ncol = 3,
                              rel_widths = c(1, 1, .9))


save_plot("plots/horizontal_three_panel_plot_CforN_relationships_fig1.pdf", 
          threepanels_horiz,
          base_width = 15)

threepanels_vertical = plot_grid(ByfungusCforNmyco_nolegend, 
                                 ByfungusNforNroots_nolegend,
                                 hyphalCformycoC_plot_legendbelow,
                                 labels = c("a", "b", "c"),
                                 align = "v",
                                 axis = "l",
                                 nrow = 3,
                                 ncol = 1,
                                 rel_heights = c(1,1,1))

save_plot("plots/vertical_three_panel_plot_CforN_relationships_fig1.pdf", 
          threepanels_vertical,
          base_height = 12,
          base_width = 5)



### Two panels including extreme value ####

ByfungusCforNmyco_ie = ggplot(data = nitrogeninfo) +
  geom_point(aes(x = mycologN15,
                 y = mycologC13, 
                 color = N_level)) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              data = subset(nitrogeninfo, compartment_fungus == "Tt"),
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

ByfungusNforNroots_ie = ggplot(data = nitrogeninfo) +
  geom_point(aes(x = mycologN15,
                 y = nmlogN15, 
                 color = N_level)) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              data = subset(nitrogeninfo, compartment_fungus == "Tt"),
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

ByfungusCforNmyco_ie_nolegend = ByfungusCforNmyco_ie +
  theme(legend.position = "none")

ByfungusNforNroots_ie_legendbelow = ByfungusNforNroots_ie +
  theme(legend.position = "bottom")

twopanels_vertical = plot_grid(ByfungusCforNmyco_ie_nolegend, 
                                 ByfungusNforNroots_ie_legendbelow,
                                 labels = c("a", "b"),
                                 align = "v",
                                 axis = "l",
                                 nrow = 2,
                                 ncol = 1,
                                 rel_heights = c(1,1))

save_plot("plots/vertical_TWO_panel_plot_CforN_relationships_includingextreme_forsupp.pdf", 
          twopanels_vertical,
          base_height = 8,
          base_width = 5)
