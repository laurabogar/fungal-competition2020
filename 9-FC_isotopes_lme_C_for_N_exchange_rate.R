# 7 - FC_isotopes_lme_C_for_N_exchange_rate
# June 2020
# Using an LME framework to more rigorously calculate how C for N 
# exchange rates varied in this experiment.

# setwd("~/Documents/Fungal competition project/fungal-competition2020/")

# Load required packages
library(agricolae)
library(cowplot)
library(emmeans)
library(tidyverse)
library(lme4)
library(lmerTest)
library(stargazer)

# Loading required data
# forN = read_csv("processeddata/isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates_withmixed.csv")
# forN = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly_INCLUDING_MIXED.csv")
forN = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly_INCLUDING_MIXED_betterversus.csv")


forN = subset(forN, received15N == "Y")
forN$mycoCforN = forN$mycoC13ppmexcess/forN$mycoN15ppmexcess
forN$improvedCforN = forN$mycoC13ppmexcess/forN$nmN15ppmexcess
# log rule: log(x/y) = log(x) - log(y)
forN$logmycoCforN = forN$mycologC13 - forN$mycologN15
forN$logimprovedCforN = forN$mycologC13 - forN$nmlogN15

justmycos = forN[!is.na(forN$logmycoCforN),]
justmycos = subset(justmycos, mycofungus != "None")

justmycos$mycofungus = factor(justmycos$mycofungus, levels = c("Sp", "Tt"))
justmycos$compartment_fungus = as.factor(justmycos$compartment_fungus)

justmycos_nomixed = justmycos[-grep("MIXED", justmycos$competitors),]

#### Statistical analyses ####

# No competition
# CforN = lmer(logmycoCforN ~ mycofungus * N_level + (1|Batch/Plant),
#              data = justmycos)

CforN = lmer(logmycoCforN ~ mycofungus * N_level + (1|Batch),
             data = justmycos_nomixed)

summary(CforN)

CforN.anova = anova(CforN)

CforNposthoc = emmeans(CforN, list(pairwise ~ mycofungus*N_level), adjust = "tukey")

sink("stats_tables/CforN_exchangerates_anova.html")

stargazer(CforN.anova, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()

sink("stats_tables/CforN_exchangerates_anova_posthoc.txt")

CforNposthoc

sink()

# WITH competition
# CforN = lmer(logmycoCforN ~ mycofungus * N_level * versus3 + (1|Batch/Plant),
#              data = justmycos)

CforN = lmer(logmycoCforN ~ compartment_fungus * N_level * versus + (1|Batch),
             data = justmycos_nomixed)

summary(CforN)

CforN.anova = anova(CforN)

CforNposthoc = emmeans(CforN, list(pairwise ~ compartment_fungus*N_level), adjust = "tukey")

sink("stats_tables/CforN_exchangerates_anova.html")

stargazer(CforN.anova, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()

sink("stats_tables/CforN_exchangerates_anova_posthoc.txt")

CforNposthoc

sink()

# Excluding accidental Tt contamination
nocontam = read_csv("processeddata/Plants_with_no_Tt_contamination.csv")
onlyclean = mutate(justmycos_nomixed, 
                   clean = justmycos_nomixed$Plant %in% nocontam$Plant) %>%
  filter(clean == TRUE)



CforN.onlyclean = lmer(logmycoCforN ~ compartment_fungus * N_level * versus + (1|Batch),
             data = onlyclean) # rank deficient, dropping one column/coefficient

summary(CforN.onlyclean) # here, it's just N level and marginally fungal species, no sig. interaction

CforN.anova = anova(CforN.onlyclean)

CforNposthoc = emmeans(CforN.onlyclean, list(pairwise ~ compartment_fungus*N_level), adjust = "tukey")

sink("stats_tables/CforN_exchangerates_anova_onlyclean.html")

stargazer(CforN.anova, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()

sink("stats_tables/CforN_exchangerates_anova_posthoc_onlyclean.txt")

CforNposthoc

sink()

### Including root N instead of myco N ####

CforN_improved = lmer(logimprovedCforN ~ mycofungus * N_level * versus3 + (1|Batch),
             data = justmycos_nomixed)

CforN_improvedlm = lm(logimprovedCforN ~ N_level * versus3 * mycofungus,
                        data = justmycos_nomixed)
ggPredict(CforN_improvedlm)

summary(CforN_improved)

CforN.anova_improved = anova(CforN_improved)

CforNposthoc = emmeans(CforN_improved, list(pairwise ~ mycofungus*N_level), adjust = "tukey")

sink("stats_tables/CforN_exchangerates_anova.html")

stargazer(CforN.anova, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()

sink("stats_tables/CforN_exchangerates_anova_posthoc.txt")

CforNposthoc

sink()

# tx = with(justmycos, interaction(mycofungus, N_level))
# forlabels = aov(logmycoCforN~tx, data = justmycos)
# mylabels = HSD.test(forlabels, "tx", group = TRUE)
# These match my lme posthoc results above.

### Plot exchange rates (plant C to fungal N) in mycorrhizas ###
labels = c(High = "High N", Low = "Low N")
annotations = data.frame(x = c((1:2), (1:2)),
                         y = c(4.3, 4, 3.3, 3.3),
                         N_level = c(rep("High", 2), rep("Low", 2)),
                         labs = c(paste(c("a", "b")), paste(c("b", "b"))),
                         x1 = c(0.6, 1.6, 0.6, 1.6), 
                         x2 = c(1.4, 2.4, 1.4, 2.4), 
                         y1 = c(4, 3.7, 3, 3), 
                         y2 = c(4, 3.7, 3, 3))


exchangerate_plot = ggplot(data = justmycos_nomixed) +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(.9),
               aes(x = mycofungus, 
                   y = logmycoCforN,
                   fill = versus3)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.9),
             aes(x = mycofungus, 
                 y = logmycoCforN,
                 fill = versus3)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Log exchange rate (plant C to\nfungal N in mycorrhizas)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  scale_fill_manual(values = c("lightgray", "gray42", "white")) +
  labs(fill = "Competitor") +
  geom_segment(data = annotations, aes(x = x1, xend = x2,
                                     y = y1, yend = y2),
               colour = "black") +
  geom_text(data = annotations, aes(x, y, label = labs))


save_plot("plots/CforN_exchangerates_withcomp.pdf", 
          exchangerate_plot,
          ncol = 1,
          base_aspect_ratio = 1.8)

### ROOT N INSTEAD: plot exchange rates (plant C to fungal N) in mycorrhizas ###
labels = c(High = "High N", Low = "Low N")
annotations = data.frame(x = c((1:2), (1:2)),
                         y = c(4.3, 4, 3.3, 3.3),
                         N_level = c(rep("High", 2), rep("Low", 2)),
                         labs = c(paste(c("a", "b")), paste(c("b", "b"))),
                         x1 = c(0.6, 1.6, 0.6, 1.6), 
                         x2 = c(1.4, 2.4, 1.4, 2.4), 
                         y1 = c(4, 3.7, 3, 3), 
                         y2 = c(4, 3.7, 3, 3))


exchangerate_plot_improved = ggplot(data = justmycos_nomixed) +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(.9),
               aes(x = mycofungus, 
                   y = logimprovedCforN,
                   fill = versus3)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15),
             aes(x = mycofungus, 
                 y = logimprovedCforN,
                 fill = versus3)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Log exchange rate (plant C\nin mycorrhizas to fungal N in roots") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  scale_fill_manual(values = c("lightgray", "gray42", "white")) +
  labs(fill = "Competitor") 
  # geom_segment(data = annotations, aes(x = x1, xend = x2,
  #                                      y = y1, yend = y2),
  #              colour = "black") +
  # geom_text(data = annotations, aes(x, y, label = labs))


save_plot("plots/CforN_exchangerates_withcomp.pdf", 
          exchangerate_plot,
          ncol = 1,
          base_aspect_ratio = 1.8)


#### Building three panel plot ####
together = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly_INCLUDING_MIXED_betterversus.csv")
nonm = together[!is.na(together$mycorrhizas.APE13C),]
nonm = subset(nonm, compartment_fungus != "None")
excluding_mixed = nonm[-grep("MIXED", nonm$competitors),]

nocontam = read_csv("processeddata/Plants_with_no_Tt_contamination.csv")
excluding_mixed_onlyclean = mutate(excluding_mixed, clean = excluding_mixed$Plant %in% nocontam$Plant)
excluding_mixed_onlyclean = subset(excluding_mixed_onlyclean, clean == TRUE)

# nonm$versus2 = nonm$versus
# nonm$versus3 = nonm$versus
# for (i in 1:nrow(nonm)) {
#   if ((nonm$versus[i] == "Mixed")|
#       (nonm$compartment_fungus[i] == "MIXED")) {
#     nonm$versus2[i] = "Other"
#     nonm$versus3[i] = "Other"
#   } else if ((nonm$versus[i] == "Tt" & nonm$compartment_fungus[i] == "Tt")|
#              (nonm$versus[i] == "Sp" & nonm$compartment_fungus[i] == "Sp")) {
#     nonm$versus2[i] = "Self"
#     nonm$versus3[i] = "Self"
#   } else if ((nonm$versus[i] == "Tt" & nonm$compartment_fungus[i] == "Sp")|
#              (nonm$versus[i] == "Sp" & nonm$compartment_fungus[i] == "Tt")) {
#     nonm$versus2[i] = "Other"
#     nonm$versus3[i] = "Other"
#   } else if (nonm$versus[i] == "None") {
#     nonm$versus2[i] = "Other"
#     nonm$versus3[i] = "None"
#     
#   }
# }
# excluding_mixed = nonm[-grep("MIXED", nonm$competitors),]
# 
# ndata = subset(nonm, received15N == "Y")
# ndata = subset(ndata, compartment_fungus != "None")
# justmycos = subset(ndata, mycofungus != "None")

labels = c(High = "High N", Low = "Low N")

### Carbon plot ####
annotations = data.frame(x = c((1:2), (1:2)),
                         y = c(7.3, 7.8, 7.1, 6.3),
                         N_level = c(rep("High", 2), rep("Low", 2)),
                         labs = c(paste(c("a", "a")), paste(c("a", "b"))),
                         x1 = c(0.6, 1.6, 0.6, 1.6), 
                         x2 = c(1.4, 2.4, 1.4, 2.4), 
                         y1 = c(7, 7.5, 6.8, 6))

carboncomparison_withcompetition = ggplot(data = excluding_mixed) +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(.9),
               aes(x = compartment_fungus, 
                   y = mycologC13,
                   fill = versus3)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.9),
             aes(x = compartment_fungus, 
                 y = mycologC13,
                 fill = versus3)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Plant C in mycorrhizas\n(ln ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  # geom_text(data = annotations, aes(x, y, label = labs))
  scale_fill_manual(values = c("lightgray", "gray42", "white")) +
  geom_segment(data = annotations, aes(x = x1, xend = x2,
                                   y = y1, yend = y1),
               colour = "black") +
  geom_text(data = annotations, aes(x, y, label = labs)) +
  labs(fill = "Competitor")


### Nitrogen plot ####
collabels = data.frame(N_level = c("High", "High", "Low", "Low"),
                       x1 = c(0.6, 1.6, 0.6, 1.6), 
                       x2 = c(1.4, 2.4, 1.4, 2.4), 
                       y1 = c(4, 6.5, 5.6, 4), 
                       y2 = c(4, 6.5, 5.6, 4),
                       xstar = c(1, 2, 1, 2), ystar = c(4.3, 6.8, 5.9, 4.3),
                       lab = c("a", "b", "ab", "a"))

# margsig = data.frame(N_level = "High",
#                      x1 = 1,
#                      x2 = 2,
#                      xstar = 1.5,
#                      y1 = 7.2,
#                      y2 = 7.2,
#                      ystar = 7.9,
#                      lab = ".")

labels = c(High = "High N", Low = "Low N")
nitrogencomparison_mycos_withcomp = ggplot(data = justmycos_nomixed) +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(.9),
               aes(x = mycofungus, 
                   y = mycologN15,
                   fill = versus3)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.9),
             aes(x = mycofungus, 
                 y = mycologN15,
                 fill = versus3)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Fungal N in mycorrhizas\n(ln ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  # geom_text(data = annotations, aes(x, y, label = labs)) +
  scale_fill_manual(values = c("lightgray", "gray42", "white")) +
  labs(fill = "Competitor") +
  geom_text(data = collabels, aes(x = xstar,  y = ystar, label = lab)) +
  geom_segment(data = collabels, aes(x = x1, xend = x2,
                                     y = y1, yend = y2),
               colour = "black")
  # geom_text(data = margsig, size = 10, aes(x = xstar,  y = ystar, label = lab)) +
  # geom_segment(data = margsig, aes(x = x1, xend = x2,
  #                                  y = y1, yend = y2),
  #              colour = "black")

carboncomparison_withcompetition_nolegend = carboncomparison_withcompetition +
  theme(legend.position = "none")

nitrogencomparison_mycos_withcomp_nolegend = nitrogencomparison_mycos_withcomp +
  theme(legend.position = "none")

exchangerate_plot_legendbelow = exchangerate_plot +
  theme(legend.position = "bottom")

threepanels_horiz = plot_grid(carboncomparison_withcompetition_nolegend, 
                        nitrogencomparison_mycos_withcomp_nolegend,
                        exchangerate_plot,
                        labels = c("a", "b", "c"),
                        align = "h",
                        # nrow = 1,
                        ncol = 3,
                        rel_widths = c(1, 1, 1.6))


save_plot("plots/horizontal_three_panel_plot_CC_NN_CN_withcompetition_updated.pdf", 
          threepanels_horiz)

threepanels_vertical = plot_grid(carboncomparison_withcompetition_nolegend, 
                              nitrogencomparison_mycos_withcomp_nolegend,
                              exchangerate_plot_legendbelow,
                              labels = c("a", "b", "c"),
                              align = "v",
                              nrow = 3,
                              ncol = 1,
                              rel_heights = c(1, 1, 1.2))



save_plot("plots/vertical_three_panel_plot_CC_NN_CN_withcompetition_updated.pdf", 
          threepanels_vertical,
          base_height = 10,
          base_aspect_ratio = .5)

#### Three panel plot with ONLY clean plants (no Tt contam) ####

together = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly_INCLUDING_MIXED_betterversus.csv")
nonm = together[!is.na(together$mycorrhizas.APE13C),]
nonm = subset(nonm, compartment_fungus != "None")
excluding_mixed = nonm[-grep("MIXED", nonm$competitors),]

nocontam = read_csv("processeddata/Plants_with_no_Tt_contamination.csv")
excluding_mixed_onlyclean = mutate(excluding_mixed, clean = excluding_mixed$Plant %in% nocontam$Plant)
excluding_mixed_onlyclean = subset(excluding_mixed_onlyclean, clean == TRUE)



labels = c(High = "High N", Low = "Low N")

### Carbon plot ####

carboncomparison_withcompetition = ggplot(data = excluding_mixed_onlyclean) +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(.9),
               aes(x = compartment_fungus, 
                   y = mycologC13,
                   fill = versus3)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.9),
             aes(x = compartment_fungus, 
                 y = mycologC13,
                 fill = versus3)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Plant C in mycorrhizas\n(ln ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  scale_fill_manual(values = c("lightgray", "gray42", "white")) +
  labs(fill = "Competitor")


### Nitrogen plot ####

labels = c(High = "High N", Low = "Low N")
nitrogencomparison_mycos_withcomp = ggplot(data = onlyclean) +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(.9),
               aes(x = mycofungus, 
                   y = mycologN15,
                   fill = versus3)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.9),
             aes(x = mycofungus, 
                 y = mycologN15,
                 fill = versus3)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Fungal N in mycorrhizas\n(ln ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  scale_fill_manual(values = c("lightgray", "gray42", "white")) +
  labs(fill = "Competitor")

### Exchange rate plot ###

exchangerate_plot = ggplot(data = onlyclean) +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(.9),
               aes(x = mycofungus, 
                   y = logmycoCforN,
                   fill = versus3)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.9),
             aes(x = mycofungus, 
                 y = logmycoCforN,
                 fill = versus3)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Log exchange rate (plant C to\nfungal N in mycorrhizas)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  scale_fill_manual(values = c("lightgray", "gray42", "white")) +
  labs(fill = "Competitor")


carboncomparison_withcompetition_nolegend = carboncomparison_withcompetition +
  theme(legend.position = "none")

nitrogencomparison_mycos_withcomp_nolegend = nitrogencomparison_mycos_withcomp +
  theme(legend.position = "none")

exchangerate_plot_legendbelow = exchangerate_plot +
  theme(legend.position = "bottom")


threepanels_vertical = plot_grid(carboncomparison_withcompetition_nolegend, 
                                 nitrogencomparison_mycos_withcomp_nolegend,
                                 exchangerate_plot_legendbelow,
                                 labels = c("a", "b", "c"),
                                 align = "v",
                                 nrow = 3,
                                 ncol = 1,
                                 rel_heights = c(1, 1, 1.2))



save_plot("plots/vertical_three_panel_plot_CC_NN_CN_withcompetition_onlyclean.pdf", 
          threepanels_vertical,
          base_height = 10,
          base_aspect_ratio = .5)



#### Looking at percent N ####
# Overall, the patterns in percent N look pretty similar to the patterns
# that we can see with the N15 label. Overall, Tt (whether
# mycos or adjacent fine roots) has the most N (in this case,
# as a percent of tissue contents) in the high N tx, trailed
# by Sp in high N, then Sp in low N and Tt in low N.
# The percent N data indicate possible hoarding by Tt in high N,
# when competing versus Suillus? (lower percent N in fine roots
# than in other competition treatments)




includingall = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly_INCLUDING_MIXED_and_pctCN.csv")
all_nomix = includingall[-grep("MIXED", includingall$competitors),]
all_nomix = all_nomix[!is.na(all_nomix$competitors),]
# problem: this includes each plant twice.

forpct = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly_INCLUDING_MIXED_and_pctCN.csv")

# forpct = read_csv("processeddata/isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates_withpctN.csv")
forpct_nomix = forpct[-grep("MIXED", forpct$competitors),]

granular_data_bycompt = read_csv("processeddata/granular_mass_and_colonization_data_by_compartment.csv")
granular_nomix = granular_data_bycompt[-grep("MIXED", granular_data_bycompt$competitors),]
granular_nomix$Side = tolower(granular_nomix$Side)
granular_nomix = subset(granular_nomix, compartment_fungus != "Failed")
granular_tojoin = select(granular_nomix, Plant, Side, 
                         total_root_biomass_compartment,
                         compartment_fungus)

tojoin = select(forpct_nomix, Plant, Side, uncolonized_roots.pctN, mycorrhizas.pctN,
                mycofungus, competition_treatment, N_level)
scalingNbymass = left_join(tojoin, granular_tojoin) %>% distinct()
scalingNbymass$g_N = scalingNbymass$uncolonized_roots.pctN*scalingNbymass$total_root_biomass_compartment

scalingNbymass = subset(scalingNbymass, 
                        competition_treatment != "Mixed" &
                          competition_treatment != 0)

labels = c(High = "High N", Low = "Low N")
nitrogenbymass = ggplot(data = scalingNbymass) +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(.9),
               aes(x = compartment_fungus, 
                   y = g_N)) +
  geom_jitter(aes(x = compartment_fungus, 
                 y = g_N)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Total g N in roots\n(by compartment)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") 
  # geom_text(data = annotations, aes(x, y, label = labs)) +
  # scale_fill_manual(values = c("lightgray", "gray42", "white")) +
  # labs(fill = "Competitor")

labels = c(High = "High N", Low = "Low N")
nitrogen_withcomp = ggplot(data = forpct_nomix) +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(.9),
               aes(x = mycofungus, 
                   y = mycorrhizas.pctN,
                   fill = versus)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.9),
             aes(x = mycofungus, 
                 y = mycorrhizas.pctN,
                 fill = versus)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Percent N in mycorrhizas") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  # geom_text(data = annotations, aes(x, y, label = labs)) +
  scale_fill_manual(values = c("lightgray", "gray42", "white")) +
  labs(fill = "Competitor")

nitrogenhoarding_mycos_withcomp = ggplot(data = forpct_nomix) +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(.9),
               aes(x = mycofungus, 
                   y = log(mycorrhizas.pctN/uncolonized_roots.pctN),
                   fill = versus)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.9),
             aes(x = mycofungus, 
                 y = log(mycorrhizas.pctN/uncolonized_roots.pctN),
                 fill = versus)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Log ratio of N in mycorrhizas\nto N in adjacent roots") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  # geom_text(data = annotations, aes(x, y, label = labs)) +
  scale_fill_manual(values = c("lightgray", "gray42", "white")) +
  labs(fill = "Competitor")

percentN_roots_withcomp = ggplot(data = forpct_nomix) +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(.9),
               aes(x = mycofungus, 
                   y = uncolonized_roots.pctN,
                   fill = versus)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.9),
             aes(x = mycofungus, 
                 y = uncolonized_roots.pctN,
                 fill = versus)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Percent N in fine roots\nadjacent to mycorrhizas") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  # geom_text(data = annotations, aes(x, y, label = labs)) +
  scale_fill_manual(values = c("lightgray", "gray42", "white")) +
  labs(fill = "Competitor")

# Can I look at total N in the root system?
tojoin = select(newcomp, Plant, Side, comp_root_mass)
tojoin$Side = tolower(tojoin$Side)
Nwithmass = left_join(forpct, tojoin)
Ntouse = Nwithmass %>% distinct(Plant, Side, 
                                mycorrhizas.pctN, 
                                uncolonized_roots.pctN,
                                comp_root_mass)


labels = c(High = "High N", Low = "Low N")
nitrogenpctcomparison_mycos_withcomp = ggplot(data = Nwithmass) +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(.9),
               aes(x = mycofungus, 
                   y = mycorrhizas.pctN,
                   fill = versus)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.9),
             aes(x = mycofungus, 
                 y = mycorrhizas.pctN,
                 fill = versus)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Percent N in mycorrhizas") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  # geom_text(data = annotations, aes(x, y, label = labs)) +
  scale_fill_manual(values = c("lightgray", "gray42", "white")) +
  labs(fill = "Competitor")

#### JUST ROOT N: Nitrogen plot ####

collabels = data.frame(N_level = c("High", "High", "Low", "Low"),
                       x1 = c(0.6, 1.6, 0.6, 1.6), 
                       x2 = c(1.4, 2.4, 1.4, 2.4), 
                       y1 = c(4, 6.5, 5.6, 4), 
                       y2 = c(4, 6.5, 5.6, 4),
                       xstar = c(1, 2, 1, 2), ystar = c(4.3, 6.8, 5.9, 4.3),
                       lab = c("ab", "b", "ab", "a"))

margsig = data.frame(N_level = "High",
                     x1 = 1,
                     x2 = 2,
                     xstar = 1.5,
                     y1 = 7.2,
                     y2 = 7.2,
                     ystar = 7.9,
                     lab = ".")

labels = c(High = "High N", Low = "Low N")
nitrogencomparison_mycos_withcomp_improved = ggplot(data = justmycos_nomixed) +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(.9),
               aes(x = mycofungus, 
                   y = nmlogN15,
                   fill = versus3)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15),
             aes(x = mycofungus, 
                 y = nmlogN15,
                 fill = versus3)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Labeled N in uncolonized roots\n(ln ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  # geom_text(data = annotations, aes(x, y, label = labs)) +
  scale_fill_manual(values = c("lightgray", "gray42", "white")) +
  labs(fill = "Competitor") 
  # geom_text(data = collabels, aes(x = xstar,  y = ystar, label = lab)) +
  # geom_segment(data = collabels, aes(x = x1, xend = x2,
  #                                    y = y1, yend = y2),
  #              colour = "black") +
  # geom_text(data = margsig, size = 10, aes(x = xstar,  y = ystar, label = lab)) +
  # geom_segment(data = margsig, aes(x = x1, xend = x2,
  #                                  y = y1, yend = y2),
  #              colour = "black")

#### HOARDING ? ####

myhoardmodel = lm(nmlogN15 ~ mycologN15*mycofungus*N_level, data = justmycos_nomixed)
ggPredict(myhoardmodel)

#### OLD STUFF ####

annotations = data.frame(x = c((1:2), (1:2)),
                         y = c(7.2, 7.2, 7.2, 6),
                         N_level = c(rep("High", 2), rep("Low", 2)),
                         labs = c(paste(c("a", "a")), paste(c("a", "b"))))


carboncomparison = ggplot(data = excluding_mixed) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = compartment_fungus, y = mycologC13)) +
  geom_jitter(width = 0.20,
              aes(x = compartment_fungus, 
                  y = mycologC13)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Plant C in mycorrhizas\n(ln ppm excess)") +
  # theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  geom_text(data = annotations, aes(x, y, label = labs))

carboncomparison_nolegend = carboncomparison +
  theme(legend.position = "none")

labels = c(High = "High N", Low = "Low N")
annotations = data.frame(x = c((1:2), (1:2)),
                         y = c(4.5, 6.6, 5.6, 4.5),
                         N_level = c(rep("High", 2), rep("Low", 2)),
                         labs = c(paste(c("a", "b")), paste(c("ab", "a"))))

nitrogencomparison_mycos = ggplot(data = justmycos) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = mycofungus, 
                   y = mycologN15)) +
  geom_jitter(width = 0.20,
              aes(x = mycofungus, 
                  y = mycologN15)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Labeled N in mycorrhizas\n(ln ppm excess)") +
  # theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  geom_text(data = annotations, aes(x, y, label = labs))

nitrogencomparison_mycos_nolegend = nitrogencomparison_mycos +
  theme(legend.position = "none")

threepanels = plot_grid(carboncomparison_nolegend, 
                        nitrogencomparison_mycos_nolegend,
                        exchangerate_plot,
                        labels = c("A", "B", "C"),
                        align = "v",
                        nrow = 3)

save_plot("plots/vertical_three_panel_plot_CC_NN_CN.pdf", 
          threepanels,
          base_height = 10,
          base_aspect_ratio = .5)

#### OLD UNNECESSARY STUFF FROM HERE ON ####

# forN$mycoCforN = numeric(nrow(forN))
# for (i in 1:nrow(forN)) {
#    for (j in 1:nrow(forN)) {
#      if (forN$Plant[i] == forN$Plant[j] &
#          forN$Side[i] == forN$Side[j]) {
#        if (forN$mycofungus[i] == "None" &
#            forN$mycofungus[j] != "None") {
#          forN$mycoCforN = forN$mycoC13ppmexcess[j]/forN$mycoN15ppmexcess
#        }
#      }
# }
# together = read_csv("./FCdata/isotope_and_plant_metadata_with_competition_coded_clearly.csv")
# exchangerates = subset(together, received15N == "Y" & Batch != "NA" & mycorrhizas.APE13C != "NA" & mycorrhizas.APE15N != "NA")
# 
# # Constants
# 
# forcefactor = -min(exchangerates$mycorrhizas.APE15N) + 0.000001 # need to coerce small negative N enrichment values to a near-zero value compatible with ratios and log transformation
# myoutliers = c(6032, 6090) # just outlandishly high C for N ratios; >= 1 order of magnitude greater than other values
# 
# ### Preparing data and calculations ###
# 
# exchangerates$forced.mycorrhizas.APE15N = exchangerates$mycorrhizas.APE15N + forcefactor
# exchangerates$forced.mycorrhizas.ppm15N = exchangerates$forced.mycorrhizas.APE15N * 10^4
#  
# exchangerates$forced.mycoC13forN15 = exchangerates$mycorrhizas.APE13C/exchangerates$forced.mycorrhizas.APE15N

# toexamine = exchangerates[exchangerates$Plant %in% myoutliers,]
# 
# toexamine$mycorrhizas.APE13C
# toexamine$mycorrhizas.APE15N # Okay, both of these
# had the same APE value, because (looking back at
# original UCSC data) they had the same d15N: 3.86.

# exchangerates_withoutliers = exchangerates

# exchangerates_nooutliers = exchangerates[!exchangerates$Plant %in% myoutliers,]

# mymax = max(exchangerates_nooutliers$forced.mycoC13forN15)
# outlierfudgefactor = mymax

# I know the fungi were getting a lot of 13C for very little 15N
# in these outliers.
# The ratio makes them artificially high values.
# Let's coerce them to the maximum value observed in other samples,
# since I think their value ("high") is still scientifically meaningful, even
# though the absolute ratio is too large to analyze.

# exchangerates$forced.mycoC13forN15[exchangerates$Plant %in% myoutliers] = outlierfudgefactor
# exchangerates$Batch = as.factor(exchangerates$Batch)

#### Did the fungi get different exchange rates at diff N levels? ####

CforN.full = lmer(log(mycoC13forN15) ~ compartment_fungus*versus*N_level + (1|Batch), 
                  data = exchangerates) # fixed-effect model matrix is rank deficient so dropping 3 columns/coefficients

anova(CforN.full) # Compartment fungus is marginally significant (p = 0.077),
# while the interaction between that and the N level is highly significant.

# Note: I tried this same test using "competition treatment"
# instead of "versus," reasoning that maybe a fungus has a 
# fundamentally different experience if it is competing
# with a heterospecific competitor... but indeed, that version
# (with competition treatment instead of straight up competitor identity)
# was qualitatively the same -- competition treatment may even be less significant than versus.

### Plot exchange rates (plant C to fungal N) in mycorrhizas ###
labels = c(High = "High N", Low = "Low N")

ggplot(data = exchangerates) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = compartment_fungus, y = log(mycoC13forN15))) +
  geom_jitter(aes(x = compartment_fungus, y = log(mycoC13forN15)),
              width = 0.25) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  ylab("Exchange rate\n(plant C to fungal N in mycorrhizas)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus")


### How does N provisioning vary by fungus and N level? ###

justN.full = lmer(log(mycoN15ppmexcess) ~ compartment_fungus*versus*N_level + (1|Batch), 
                  data = exchangerates)

anova(justN.full) # only compartment_fungus:N_level is significant.
# Nothing else even marginal.

labels = c(High = "High N", Low = "Low N")

ggplot(data = exchangerates) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = compartment_fungus, y = (mycoN15ppmexcess))) +
  geom_jitter(aes(x = compartment_fungus, y = (mycoN15ppmexcess)),
              width = 0.25) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  ylab("Fungal N in mycorrhizas (ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus")

### Old stuff ###
# Oh man, these only are true if I exclude the outliers:
annotations = data.frame(x = c((1:2), (1:2)),
                         y = c(5.5, 4, 4, 5.5),
                         N_level = c(rep("High", 2), rep("Low", 2)),
                         labs = c("a", "b", "b", "ab"))

exchangerate_plot = ggplot(data = exchangerates) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = compartment_fungus, 
                   y = log(mycoC13forN15))) + 
  geom_jitter(aes(x = compartment_fungus,
                  y = log(mycoC13forN15)),
              width = 0.25) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  xlab("Fungus receiving nitrogen label") +
  ylab(bquote(atop("Log exchange rate", "("^13*"C to "^15*"N excess in mycorrhizas)"))) +
  geom_text(data = annotations, aes(x, y, label = labs))


exchangerate_test = aov(log(mycoC13forN15) ~ compartment_fungus*N_level, data = exchangerates)
plot(exchangerate_test) #blech

# 25 and 15 look like real outliers.
# These are 6055 and 6031.
# Okay, but removing those, now 18, 22, and 24 are problems
# These are 6043, 6050, 6064 # Haha okay now all the suillus high N are gone.
# Sticking with only removing the crazy outliers
# that were >2 orders of magnitude removed from
# the median rest of the data. (e.g. >1000 ratio values)

# Definitely not quite normal, but anova is pretty robust to that.
exchangerate_test_lm = lm(log(mycoC13forN15) ~ compartment_fungus*N_level, data = exchangerates)
exchangerate_test_orderflipped = aov(log(mycoC13forN15) ~ N_level*compartment_fungus, data = exchangerates)

summary(exchangerate_test) # significant interaction.
summary(exchangerate_test_orderflipped) # Same result as when I put comp fungus first. Interaction sinificant.
summary(exchangerate_test_lm) # everything's significant if I just use lm

TukeyHSD(exchangerate_test)

median(exchangerates$mycoC13forN15[exchangerates$compartment_fungus == "Sp" & exchangerates$N_level == "High"])
median(exchangerates$mycoC13forN15[exchangerates$compartment_fungus == "Tt" & exchangerates$N_level == "High"])


exchangerate_tukey = TukeyHSD(exchangerate_test)
write.csv(exchangerate_tukey$N_level, "Statistical_tables/exchangerate_Tukey_output_Nlevel.csv")
write.csv(exchangerate_tukey$compartment_fungus, "Statistical_tables/exchangerate_Tukey_output_fungus.csv")
write.csv(exchangerate_tukey$`compartment_fungus:N_level`, "Statistical_tables/exchangerate_Tukey_output_Nlevel-Fungus.csv")

# pdf("plots/exchange_rate_boxplot_nooutliers.pdf", width = 7, height = 5)
# exchangerate_plot
# dev.off()

save_plot("plots/MAIN_Exchange_rate_boxplot_by_fungus.pdf", exchangerate_plot)


pdf("plots/exchange_rate_boxplot_fortalk.pdf", width = 8, height = 4)
exchangerate_plot
dev.off()

save_plot("plots/MAIN_Exchange_rate_boxplot_by_fungus.pdf", 
          exchangerate_plot,
          base_aspect_ratio = 1.4)


tx = with(exchangerates, interaction(N_level, compartment_fungus))

forlabels = aov(log(mycoC13forN15) ~ tx, data = exchangerates)
mylabels = HSD.test(forlabels, "tx", group = TRUE)

#### Exchange rates by competition treatment ####

exchangerates$competition_treatment = numeric(nrow(exchangerates))
for (i in 1:nrow(exchangerates)) {
  if (exchangerates$Fungi[i] == "Sp/None" |
      exchangerates$Fungi[i] == "Tt/None" |
      exchangerates$Fungi[i] == "THETE") {
    exchangerates$competition_treatment[i] = "None"
  } else if (exchangerates$Fungi[i] == "Sp/Sp" | exchangerates$Fungi[i] == "Tt/Tt") {
    exchangerates$competition_treatment[i] = "Self"
  } else if (exchangerates$Fungi[i] == "MIXED/THETE" | 
             exchangerates$Fungi[i] == "MIXED/SUIPU" |
             exchangerates$Fungi[i] == "Tt/Sp") {
    exchangerates$competition_treatment[i] = "Other"
  }
}

onlyTt = subset(exchangerates, compartment_fungus == "Tt")

competition_exchangerate_plot = ggplot(data = exchangerates) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competition_treatment,
                   y = log(mycoC13forN15))) +
  geom_jitter(width = 0.20,
              aes(x = competition_treatment,
                  y = log(mycoC13forN15),
                  shape = compartment_fungus)) +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  xlab("Competition treatment") +
  ylab(bquote(atop("Log exchange rate", "("^13*"C to "^15*"N excess in mycorrhizas)")))

save_plot("plots/MAIN_Competition_and_N_treatment_exchange_rate_plot.pdf",
          competition_exchangerate_plot,
          base_aspect_ratio = 1.4)

exchangeratebycomp = aov(log(mycoC13forN15) ~ competition_treatment*N_level, data = exchangerates)
plot(exchangeratebycomp)
summary(exchangeratebycomp)

allfactorsexchangerate = lm(log(mycoC13forN15) ~ competition_treatment*N_level*compartment_fungus, data = exchangerates)
plot(allfactorsexchangerate)
summary(allfactorsexchangerate)

justhighN = subset(exchangerates, N_level == "High")
Nanova = aov(log(mycoC13forN15) ~ competition_treatment, data = justhighN)
plot(Nanova)
summary(Nanova)

onlyTt_selfvsother = subset(onlyTt, competition_treatment != "None")
onlyTt_selfvsother$competition_treatment = as.factor(onlyTt_selfvsother$competition_treatment)

require(ggsignif)

seginfo = data.frame(x = 0.9, 
                     xend = 2.1,
                     y = 5.3, 
                     yend = 5.3, 
                     color = "black",
                     N_level = "High")

starinfo = data.frame(x = c((1.5), (1.5)),
                      y = c(5.7, 5.7),
                      N_level = c(rep("High", 1), rep("Low", 1)),
                      labs = c(".", ""))



onlyTt_selfvsother_plot = ggplot(data = onlyTt_selfvsother) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competition_treatment,
                   y = log(mycoC13forN15))) +
  geom_jitter(width = 0.25,
              aes(x = competition_treatment,
                  y = log(mycoC13forN15))) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  geom_segment(data=seginfo,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE) +
  geom_text(data = starinfo, aes(x, y, label = labs), size = 8) +
  xlab("Competition treatment") +
  theme(plot.margin = unit(c(1,1,0.5,0.5), "cm")) +
  ylab(bquote(atop("Log exchange rate: Tt", "("^13*"C to "^15*"N excess in mycorrhizas)")))

pdf("plots/MAIN_Thelephora_self_vs_other_t_test_plot_fortalk.pdf", width = 8, height = 4)
onlyTt_selfvsother_plot