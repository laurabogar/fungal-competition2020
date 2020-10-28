#5-6-FC_hyphal_13C_root_15N_comparison_figure_nooutlier.r

# setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(cowplot)
library(tidyverse)
library(stargazer)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggiraphExtra) # required for ggpredict() visualization.

carboninfo = read_csv("processeddata/data_for_carbon_only_analyses.csv")
nitrogeninfo = read_csv("processeddata/isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates.csv")

carboninfo$hyphalog13C = log(carboninfo$hyphae.ppm13Cexcess)
nitrogeninfo$hyphae.ppm13Cexcess =  nitrogeninfo$hyphae.APE13C*10^4
nitrogeninfo$hyphae.ppm15Nexcess =  nitrogeninfo$hyphae.APE15N*10^4


nitrogeninfo_nooutlier = subset(nitrogeninfo, Plant != 6041)
carboninfo_nooutlier = carboninfo[!carboninfo$hyphae.APE13C == max(carboninfo$hyphae.APE13C),] # omit outlier 6024b

#### Carbon panel: How well did hyphal C track myco C? ####

hyphalCformycoC_plot_nooutlier = ggplot(data = carboninfo_nooutlier) +
  geom_point(aes(x = mycologC13,
                 y = log(hyphae.ppm13Cexcess), 
                 color = N_level,
                 shape =compartment_fungus)) +
  geom_smooth(method = "lm", 
              aes(x = mycologC13,
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

hyphalCformycoC_plot_nolegend = hyphalCformycoC_plot_nooutlier +
  theme(legend.position = "none")

hyphaCmycoC_model = lm(hyphalog13C ~ mycologC13*compartment_fungus*N_level, data = carboninfo)
summary(hyphaCmycoC_model)
anova(hyphaCmycoC_model)
ggPredict(hyphaCmycoC_model)

hyphaCmycoC_justTt = lm(hyphalog13C ~ mycologC13*N_level, 
                        data = subset(carboninfo_nooutlier, compartment_fungus == "Tt"))
summary(hyphaCmycoC_justTt) # much clearer when you remove undersampled Sp data
ggPredict(hyphaCmycoC_justTt)
plot(hyphaCmycoC_justTt)

hyphaCmycoC_justTt_lme = lmer(hyphalog13C ~ mycologC13*N_level + (1|Batch), 
                        data = subset(carboninfo_nooutlier, compartment_fungus == "Tt"))
summary(hyphaCmycoC_justTt_lme) # mycologC13 is significant
ggPredict(hyphaCmycoC_justTt)
plot(hyphaCmycoC_justTt)

nolog = lm(hyphae.ppm13Cexcess ~ mycoC13ppmexcess*N_level, 
           data = subset(carboninfo_nooutlier, compartment_fungus == "Tt"))
plot(nolog)
summary(nolog)
ggPredict(nolog)

# When I examine histograms of each data set (hyphae and mycos C13,
# logged or not), I see long right tails on both distributions.
# The log-transformed hyphae data look sort of uniform, while the myco
# data look basically normal. Log transformation is probably the way to go.

### Nitrogen panel ####
rootNformycoN_plot_nooutlier = ggplot(data = nitrogeninfo_nooutlier) +
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

rootNformycoN_plot_nolegend = rootNformycoN_plot_nooutlier +
  theme(legend.position = "none")

# When considering N and C together, I need the dataset 
# that only includes plants/compartments that received nitrogen label.

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


model_mycoCforN = lm(mycologC13 ~ mycologN15, data = nitrogeninfo_nooutlier)
summary(model_mycoCforN)

lmer_model_mycoCforN = lmer(mycologC13 ~ mycologN15*compartment_fungus*N_level + (1|Batch),
                            data = nitrogeninfo_nooutlier)
summary(lmer_model_mycoCforN) #NS
anova(lmer_model_mycoCforN) # mycologN15 significant

lmer_model_mycoCfornmN = lmer(mycologC13 ~ nmlogN15*compartment_fungus*N_level + (1|Batch),
                            data = nitrogeninfo)
summary(lmer_model_mycoCfornmN) #NS
anova(lmer_model_mycoCfornmN) #NS

onlyTt = lmer(mycologC13 ~ nmlogN15*N_level + (1|Batch),
                              data = subset(nitrogeninfo, compartment_fungus == "Tt"))
summary(onlyTt) #Intercept AND NM N15 AND interaction with N level very significant.
anova(onlyTt)

onlySp = lmer(mycologC13 ~ nmlogN15*N_level + (1|Batch),
              data = subset(nitrogeninfo, compartment_fungus == "Sp"))
summary(onlySp) #NS.

threepanels = plot_grid(rootNformycoN_plot_nolegend, 
                        hyphalCformycoC_plot_nolegend,
                        mycoCforN,
                        labels = c("A", "B", "C"),
                        align = "h",
                        nrow = 1,
                        ncol = 3,
                        rel_widths = c(1, 1, 1.15))

save_plot("plots/Multipanel_regressions_myco_N_and_C_NOOUTLIER.jpeg",
          threepanels,
          base_aspect_ratio = 3.6)

save_plot("plots/Multipanel_regressions_myco_N_and_C_NOOUTLIER.pdf",
          threepanels,
          base_aspect_ratio = 3.6)

## Adding fourth panel? ####
data_for_hypharootplot = nitrogeninfo_nooutlier %>% drop_na(hyphae.ppm13Cexcess)
data_for_hypharootplot$hyphalogC13 = log(data_for_hypharootplot$hyphae.ppm13Cexcess)

data_for_hypharootplot$hyphalogN15 = log(data_for_hypharootplot$hyphae.ppm15Nexcess+abs(min(data_for_hypharootplot$hyphae.ppm15Nexcess))+1)
data_for_hypharootplot_justTt = subset(data_for_hypharootplot,
                                       compartment_fungus == "Tt")


# I have literally one Suillus compartment in here.
# Can't analyze effect of fungus, so best to exclude it.

# data_for_hypharootplot = subset(data_for_hypharootplot, compartment_fungus != "Sp")


hyphaCforrootN = ggplot(data = data_for_hypharootplot_justTt) +
  geom_point(aes(y = hyphalogC13,
                 x = nmlogN15, 
                 color = N_level)) +
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              aes(y = hyphalogC13,
                  x = nmlogN15),
              color = "black",
              size = 0.5) +
  ylab(bquote(atop("Hyphal "^13*"C", "(ln ppm excess)"))) +
  xlab(bquote(atop("Root "^15*"N (ln ppm excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  # scale_shape_manual(values = c(17, 15),
  #                    name = "Fungus") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

lmer_model_hyphaCfornmN = lmer(hyphalogC13 ~ nmlogN15*N_level + (1|Batch),
                            data = data_for_hypharootplot_justTt)
anova(lmer_model_hyphaCfornmN) #NS, boo
summary(lmer_model_hyphaCfornmN) # intercept significant

lmer_model_hyphalCformycoN = lmer(hyphalogC13 ~ mycologN15*N_level + (1|Batch),
                                data = data_for_hypharootplot_justTt)
anova(lmer_model_hyphalCformycoN) #NS, boo
summary(lmer_model_hyphalCformycoN) #still NS

lmer_model_mycoCfornmN = lmer(mycologC13 ~ nmlogN15*N_level + (1|Batch),
                               data = data_for_hypharootplot_justTt)
anova(lmer_model_mycoCfornmN) #NS, boo
summary(lmer_model_mycoCfornmN) #NM log N15 is significant, as is intercept

lmer_model_mycoCformycoN = lmer(mycologC13 ~ mycologN15*N_level + (1|Batch),
                              data = data_for_hypharootplot_justTt)
anova(lmer_model_mycoCformycoN) #NS, boo
summary(lmer_model_mycoCformycoN) #mycologN15 is significant.





model_hyphaerootCforN_justTt = lm(hyphalogC13 ~ nmlogN15*N_level, data = data_for_hypharootplot_justTt)
summary(model_hyphaerootCforN_justTt)
ggPredict(model_hyphaerootCforN_justTt)

model_mycorootCforN_justTt = lm(mycologC13 ~ nmlogN15*N_level, data = data_for_hypharootplot_justTt)
summary(model_mycorootCforN_justTt)
ggPredict(model_mycorootCforN_justTt)

model_rootrootCforN_justTt = lm(nmlogC13 ~ nmlogN15*N_level, data = data_for_hypharootplot_justTt)
summary(model_rootrootCforN_justTt)
ggPredict(model_rootrootCforN_justTt)

model_hyphahyphaCforN_justTt = lm(hyphalogC13 ~ hyphalogN15*N_level, data = data_for_hypharootplot_justTt)
ggPredict(model_hyphahyphaCforN_justTt)

model_hyphamycoCforN_justTt = lm(hyphalogC13 ~ mycologN15*N_level, data = data_for_hypharootplot_justTt)
ggPredict(model_hyphamycoCforN_justTt)

model_mycomycoCforN_justTt = lm(mycologC13 ~ mycologN15*N_level, data = data_for_hypharootplot_justTt)
ggPredict(model_mycomycoCforN_justTt)
# This web page is very helpful for thinking about slope differences:
# https://blog.minitab.com/blog/adventures-in-statistics-2/how-to-compare-regression-lines-between-different-models
# Looking at just Tt (and I think when I include the one Sp point),
# the intercept is almost significant (p = 0.06) and
# the slope with nmlogn15 is almost significant (p = 0.0538),
# but N level and its interaction with N15 are not.
# I think this means I should include only the significant slope in my figure?

model_HRCforN_justTt_noN = lm(hyphalog13C ~ nmlogN15, data = data_for_hypharootplot_justTt)
summary(model_HRCforN_justTt_noN) # Huh, it's a much less predictive

anova(model_hyphaerootCforN_nooutlier)

model_hyphaerootCforN = lm(log(hyphae.ppm13Cexcess) ~ nmlogN15*N_level, data = nitrogeninfo)
summary(model_hyphaerootCforN) # Literally the same
# as the outlier removed version. I think the outlier
# didn't have hyphal C and root N data together.

mytest = lm(log(hyphae.ppm13Cexcess) ~ nmlogN15, data = data_for_hypharootplot)
summary(mytest) # definitely worse to exclude N level here.

fullmycomodel = lm(mycologC13 ~ mycologN15*N_level*compartment_fungus, data = nitrogeninfo_nooutlier)
summary(fullmycomodel)
anova(fullmycomodel)
ggPredict(fullmycomodel)

test = lm(mycologC13 ~ mycologN15, data = nitrogeninfo_nooutlier)
summary(test)

justTt = lm(mycologC13 ~ mycologN15*N_level, 
            data = subset(nitrogeninfo_nooutlier, compartment_fungus == "Tt"))
summary(justTt)
anova(justTt)

justSp = lm(mycologC13 ~ mycologN15*N_level, 
            data = subset(nitrogeninfo_nooutlier, compartment_fungus == "Sp"))
summary(justSp)
anova(justSp)

sink("stats_tables/hyphae13C_vs_root_15N_lm_loglogresults.html")

stargazer(model_hyphaerootCforN_nooutlier, type = "html",
          digits = 3,
          star.cutoffs = c(0.058,0.035, 0.01, 0.001),
          digit.separator = "",
          summary = TRUE,
          report = 'vc*s*p')

sink()

# How does the relationship change if we look
# just at mycorrhizas in this same, reduced dataset?

#### C for N changes a lot depending on N level and tissue ####
model_justmycos = lm(mycologC13 ~ mycologN15*N_level, data = data_for_hypharootplot_justTt)
summary(model_justmycos)
plot(model_justmycos)

model_justmycos_nolog = lm(mycoC13ppmexcess ~ mycoN15ppmexcess*N_level, data = data_for_hypharootplot_justTt)
plot(model_justmycos_nolog)

anova(model_justmycos)
ggPredict(model_justmycos, interactive = TRUE)
ggPredict(model_hyphaerootCforN_justTt, interactive = TRUE)

model_mycosC_rootsN = lm(mycologC13 ~ nmlogN15*N_level, data = data_for_hypharootplot_justTt)
ggPredict(model_mycosC_rootsN, interactive = TRUE)



mycoCformycoN_directcomp = ggplot(data = data_for_hypharootplot) +
  geom_point(aes(y = mycologC13,
                 x = mycologN15, 
                 color = N_level,
                 shape = compartment_fungus)) +
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              aes(y = mycologC13,
                  x = mycologN15),
              color = "black",
              size = 0.5) +
  ylab(bquote(atop("Mycorrhiza "^13*"C", "(ln ppm excess)"))) +
  xlab(bquote(atop("Mycorrhiza "^15*"N (ln ppm excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

mycoCforN_justSp = ggplot(data = subset(nitrogeninfo_nooutlier,
                                        compartment_fungus == "Sp")) +
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

mycoCforN_justTt = ggplot(data = subset(nitrogeninfo_nooutlier,
                                        compartment_fungus == "Tt")) +
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


#### Recoloring the myco C for N figure to accommodate different slopes for Sp and Tt ####
summary(fullmycomodel)
# From what I can piece together here https://cran.r-project.org/web/packages/ggiraphExtra/vignettes/ggPredict.html,
# the model summary should allow me to figure out the equations
# for Sp and Tt.

# Tt 13C = -0.8718 + 2.2295 (myco15N) +
# equation1=function(x){coef(fullmycomodel)[2]*x+coef(fit1)[1]}

alternativemodel = lm(mycologC13 ~ mycologN15 * compartment_fungus * N_level,
                      data = nitrogeninfo)

ggPredict(alternativemodel, interactive = TRUE) # This is very cool, but I am not super
# interested
# WOW the interactive model is amazing!!!

ggPredict(alternativemodel)




fullmodel_rootsasN = lm(mycologC13 ~ nmlogN15 * compartment_fungus * N_level,
                        data = nitrogeninfo)
summary(fullmodel_rootsasN)
ggPredict(fullmodel_rootsasN, interactive = TRUE)

fullmodel_flippingaxes = lm(nmlogN15 ~ mycologC13 * compartment_fungus * N_level,
                        data = nitrogeninfo)
ggPredict(fullmodel_flippingaxes)


mycoCforN_recolored = ggplot(data = nitrogeninfo_nooutlier) +
  geom_point(aes(x = mycologN15,
                 y = mycologC13, 
                 shape = N_level,
                 color = compartment_fungus)) +
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              aes(x = mycologN15,
                  y = mycologC13),
              color = "black",
              size = 0.5) +
  geom_smooth(data = subset(nitrogeninfo_nooutlier, 
                            compartment_fungus == "Tt"),
              method = "lm", 
              formula = y ~ x, 
              aes(x = mycologN15,
                  y = mycologC13),
              color = "orchid4",
              size = 0.5) +
  ylab(bquote(atop("Mycorrhizal "^13*"C", "(ln ppm excess)"))) +
  xlab(bquote(atop("Mycorrhizal "^15*"N (ln ppm excess)"))) +
  scale_color_manual(values = c("thistle3", "orchid4"),
                     name = "Fungus") +
  scale_shape_manual(values = c(19, 1),
                     name = "N level") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

# CforNlm_nooutlier = lm((mycoC13ppmexcess) ~ mycoN15ppmexcess, data = nitrogeninfo_nooutlier)
# plot(CforNlm_nooutlier)
# summary(CforNlm_nooutlier)
# 
# 
# sink("stats_tables/myco13C_vs_myco_15N_lm_results_NOOUTLIER.html")
# 
# stargazer(CforNlm_nooutlier, type = "html",
#           digits = 3,
#           star.cutoffs = c(0.05, 0.01, 0.001),
#           digit.separator = "",
#           summary = FALSE)
# 
# sink()

# Use N15 values forced positive with linear transformation
# so you can try a log fit.
# CforNlm_log_nooutlier = lm((mycoC13ppmexcess) ~ mycologN15, data = nitrogeninfo_nooutlier)
# plot(CforNlm_log_nooutlier) # Looks fine
# summary(CforNlm_log_nooutlier) # Adj R^2 = 0.2586
# 
CforNlm_loglog_nooutlier = lm(mycologC13 ~ mycologN15, data = nitrogeninfo_nooutlier)
# plot(CforNlm_loglog_nooutlier) # Better
# summary(CforNlm_loglog_nooutlier) # Adj R^2 = 0.402
# # So a log-log relationship is best when we remove the extreme values.
# 
# 
# summary(CforNlm_nooutlier) # Adj R^2 = 0.1686
# 
# sink("stats_tables/myco13C_vs_myco_15N_lm_logresults_nooutlier.html")
# 
# stargazer(CforNlm_log_nooutlier, type = "html",
#           digits = 3,
#           star.cutoffs = c(0.05, 0.01, 0.001),
#           digit.separator = "",
#           summary = FALSE)
# 
# sink()

sink("stats_tables/myco13C_vs_myco_15N_lm_loglogresults_nooutlier.html")

stargazer(CforNlm_loglog_nooutlier, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE)

sink()

# sink("stats_tables/myco13C_vs_myco_15N_comparingthreemodels_nooutlier.html")
# 
# stargazer(CforNlm_nooutlier, CforNlm_log_nooutlier, CforNlm_loglog_nooutlier, 
#           type = "html",
#           align = TRUE,
#           digits = 3,
#           star.cutoffs = c(0.05, 0.01, 0.001),
#           digit.separator = "",
#           no.space = TRUE)
# 
# sink()
