# 14 - Fungal competition: Code and stats associated with a new panel for figure 2,
# looking at shoot 15N concentrations as a function of mycorrhizal 15N concentrations.

# Linear relationships between 13C and 15N in mycorrhizas;
# 15N in mycorrhizas and 15N in uncolonized roots;
# and 13C in mycorrhizas and 13C in hyphae.

library(cowplot)
library(tidyverse)
library(stargazer)
library(lme4)
library(lmerTest)
library(emmeans)
library(MuMIn)

forcefactor = 18.32042 # This is the absolute value of the minimum 15N ppm excess
# in uncolonized fine roots, as determined in script 5. I'll use it to bring 
# my leaf 15N values into line so they can be compared with the other 15N.

carboninfo = read_csv("processeddata/data_for_carbon_only_analyses_withpctC.csv") # from script 4-alternative
nitrogeninfo = read_csv("processeddata/isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates_withpctN.csv")
leafinfo = read_csv("processeddata/Cleaned_processed_FC_leaf_isotopes.csv")
moremetadata = read_csv("processeddata/percent_col_and_mass_data_by_plant.csv")

carboninfo$hyphalog13C = log(carboninfo$hyphae.ppm13Cexcess)
nitrogeninfo$hyphae.ppm13Cexcess =  nitrogeninfo$hyphae.APE13C*10^4
nitrogeninfo$hyphae.ppm15Nexcess =  nitrogeninfo$hyphae.APE15N*10^4

leafinfo$leaf.ppm13Cexcess =  leafinfo$leaf.APE13C*10^4
leafinfo$leaf.ppm15Nexcess =  leafinfo$leaf.APE15N*10^4

nocontam = read_csv("processeddata/Plants_with_no_Tt_contamination.csv")
nitrogeninfo = mutate(nitrogeninfo, clean = nitrogeninfo$Plant %in% nocontam$Plant)
carboninfo = mutate(carboninfo, clean = carboninfo$Plant %in% nocontam$Plant)
leafinfo = mutate(leafinfo, clean = leafinfo$Plant %in% nocontam$Plant)

nitrowithleaves = left_join(leafinfo, nitrogeninfo)

moremetawithleaves = left_join(leafinfo, moremetadata)
moremetawithleaves_clean = filter(moremetawithleaves, clean == TRUE)
mml_enriched = filter(moremetawithleaves, enriched == 1)

mml_enriched$Fungi[mml_enriched$Plant == 6028] = "Mixed/Tt"

# mml_unenriched = filter(moremetawithleaves, enriched == 0)

mml_enriched$leaflog13C = log(mml_enriched$leaf.ppm13Cexcess)


# The 15N values would always have needed to be transformed with the
# forcefactor from the other script (specified above)

mml_enriched$forced.leaves.N15ppmexcess = mml_enriched$leaf.ppm15Nexcess +
  forcefactor + 1

mml_enriched$leaflog15N = log(mml_enriched$forced.leaves.N15ppmexcess)

N_tojoin = select(nitrogeninfo, Plant, Side, mycofungus,
                  mycologC13:nmlogN15)

alltogether = left_join(mml_enriched, N_tojoin)

#### STATS AND PLOTS ####

### 1:  Does N15 in leaves predict C13 in leaves? ####
ByfungusCforN_leaves = ggplot(data = mml_enriched) +
  geom_point(aes(x = leaflog15N,
                 y = leaflog13C, 
                 color = N_level)) +
  # geom_smooth(method = "lm",
  #             formula = y ~ x,
  #             data = subset(nitrogeninfo, compartment_fungus == "Tt"),
  #             aes(x = mycologN15,
  #                 y = mycologC13),
  #             color = "black",
  #             size = 0.5) +
  ylab(bquote(atop("Leaf "^13*"C (ln ppm excess)"))) +
  xlab(bquote(atop("Leaf "^15*"N (ln ppm excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  # facet_grid(. ~ Fungi) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

# Not really.
# I don't think this is fundamentally the interesting question, anyway.

### 2: Does mycorrhizal 15N correlate to leaf 15N? ####

justmycos = filter(alltogether, !is.na(mycofungus))
# build plot, filtering out entries for which you have leaf but not myco isotopes
LeafNformycoN = ggplot(data = justmycos) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              data = subset(justmycos, mycofungus == "Sp"),
              aes(x = mycologN15,
                  y = leaflog15N),
              color = "black",
              size = 0.5) +
  geom_point(aes(x = mycologN15,
                 y = leaflog15N, 
                 color = N_level)) +
  xlab(bquote(atop("Mycorrhizal "^15*"N (ln ppm excess)"))) +
  ylab(bquote(atop("Leaf "^15*"N (ln ppm excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  facet_grid(. ~ mycofungus) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

leafN_mycoN_justTt = lmer(leaflog15N ~ mycologN15*N_level + (1|Batch),
                       data = filter(justmycos, mycofungus == "Tt"))
summary(leafN_mycoN_justTt)
class(leafN_mycoN_justTt) <- "lmerMod"
r.squaredGLMM(leafN_mycoN_justTt)
fitTt = as.data.frame(r.squaredGLMM(leafN_mycoN_justTt))
fitTt_R2 = round(fitTt$R2c, 4)


leafN_mycoN_justSp = lmer(leaflog15N ~ mycologN15*N_level + (1|Batch),
                          data = filter(justmycos, mycofungus == "Sp"))
summary(leafN_mycoN_justSp)
class(leafN_mycoN_justSp) <- "lmerMod"
r.squaredGLMM(leafN_mycoN_justSp)
fitSp = as.data.frame(r.squaredGLMM(leafN_mycoN_justSp))
fitSp_R2 = round(fitSp$R2c, 4)

sink("stats_tables/LeafN_for_mycoN_TtvsSp_lmer.html")

stargazer(leafN_mycoN_justTt, 
          leafN_mycoN_justSp, 
          type = "html",
          dep.var.labels = "Leaf [15N] (ln ppm excess)",
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

### 3: Does uncolonized root 15N correlate to leaf 15N? ####

justnm = filter(alltogether, !is.na(nmlogN15))
# build plot, filtering out entries for which you have leaf but not nm root isotopes
LeafNfornmN = ggplot(data = justnm) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              data = subset(justnm, mycofungus == "Sp"),
              aes(x = nmlogN15,
                  y = leaflog15N),
              color = "black",
              size = 0.5) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              data = subset(justnm, mycofungus == "Tt"),
              aes(x = nmlogN15,
                  y = leaflog15N),
              color = "black",
              size = 0.5) +
  geom_point(aes(x = nmlogN15,
                 y = leaflog15N, 
                 color = N_level)) +
  xlab(bquote(atop("Uncolonized root "^15*"N (ln ppm excess)"))) +
  ylab(bquote(atop("Leaf "^15*"N (ln ppm excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  facet_grid(. ~ mycofungus) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

leafN_nmN_justTt = lmer(leaflog15N ~ nmlogN15*N_level + (1|Batch),
                          data = filter(justnm, mycofungus == "Tt"))
summary(leafN_nmN_justTt)
class(leafN_nmN_justTt) <- "lmerMod"
r.squaredGLMM(leafN_nmN_justTt)
fitTt = as.data.frame(r.squaredGLMM(leafN_nmN_justTt))
fitTt_R2 = round(fitTt$R2c, 4)


leafN_nmN_justSp = lmer(leaflog15N ~ nmlogN15*N_level + (1|Batch),
                          data = filter(justmycos, mycofungus == "Sp"))
summary(leafN_nmN_justSp)
class(leafN_nmN_justSp) <- "lmerMod"
r.squaredGLMM(leafN_nmN_justSp)
fitSp = as.data.frame(r.squaredGLMM(leafN_nmN_justSp))
fitSp_R2 = round(fitSp$R2c, 4)

sink("stats_tables/LeafN_for_nm_root_N_TtvsSp_lmer.html")

stargazer(leafN_nmN_justTt, 
          leafN_nmN_justSp, 
          type = "html",
          dep.var.labels = "Leaf [15N] (ln ppm excess)",
          covariate.labels = c("Uncolonized root [15N] (ln ppm excess)",
                               "N level",
                               "Uncolonized root [15N]:N level"),
          column.labels = c("Tt", "Sp"),
          digits = 3,
          digit.separator = "",
          ci = TRUE,
          star.cutoffs = c(0.05, 0.01, 0.001),
          add.lines = list(c("Conditional pseudo-$R2$",
                             fitTt_R2, fitSp_R2)),
          summary = TRUE)

sink()

### Bringing together for two panel figure ####
LeafNformycoN_nolegend = LeafNformycoN +
  theme(legend.position = "none")

twopanels_horiz = plot_grid(LeafNformycoN_nolegend, 
                            LeafNfornmN,
                              labels = c("a", "b"),
                              align = "h",
                              axis = "b",
                              nrow = 1,
                              ncol = 2,
                              rel_widths = c(1, 1.2))


save_plot("plots/horizontal_two_panel_plot_leaf15N.pdf", 
          twopanels_horiz,
          base_width = 10)

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

### 4: Fungal N goes to leaves? ####

ByfungusNmycoforNleaves_ie = ggplot(data = nitrowithleaves) +
  geom_point(aes(x = mycologN15,
                 y = leaflog15N, 
                 color = N_level)) +
  xlab(bquote(atop("Mycorrhizal "^15*"N (ln ppm excess)"))) +
  ylab(bquote(atop("Leaf "^15*"N (ln ppm excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  facet_grid(. ~ compartment_fungus) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))


mycoNforleafN_justTt = lmer(leaflog15N ~ mycologN15*N_level + (1|Batch),
                       data = subset(nitrowithleaves, compartment_fungus == "Tt"))
summary(mycoNforleafN_justTt)
class(mycoNforleafN_justTt) <- "lmerMod"
r.squaredGLMM(mycoNforleafN_justTt)
fitTt = as.data.frame(r.squaredGLMM(mycoNforleafN_justTt))
fitTt_R2 = round(fitTt$R2c, 4)


mycoNforleafN_justSp = lmer(leaflog15N ~ mycologN15*N_level + (1|Batch),
                       data = subset(nitrowithleaves, compartment_fungus == "Sp"))
summary(mycoNforleafN_justSp) # singular fit
class(mycoNforleafN_justSp) <- "lmerMod"
r.squaredGLMM(mycoNforleafN_justSp)
fitSp = as.data.frame(r.squaredGLMM(mycoNforleafN_justSp))
fitSp_R2 = round(fitSp$R2c, 4)


sink("stats_tables/mycoNforleafN_TtvsSp_lmer.html")

stargazer(mycoNforleafN_justTt, 
          mycoNforleafN_justSp, 
          type = "html",
          dep.var.labels = "Leaf [15N] (ln ppm excess)",
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

ByfungusNrootforNleaves_ie = ggplot(data = nitrowithleaves) +
  geom_point(aes(x = nmlogN15,
                 y = leaflog15N, 
                 color = N_level)) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              data = subset(nitrowithleaves, compartment_fungus == "Tt"),
              aes(x = nmlogN15,
                  y = leaflog15N),
              color = "black",
              size = 0.5) +
  xlab(bquote(atop("Root "^15*"N (ln ppm excess)"))) +
  ylab(bquote(atop("Leaf "^15*"N (ln ppm excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  facet_grid(. ~ compartment_fungus) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))


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
