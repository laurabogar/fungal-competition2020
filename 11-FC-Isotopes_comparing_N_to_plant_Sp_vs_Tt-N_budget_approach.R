# Fungal competition: Code and stats associated with Figure 3b

# Comparing fungal N (15N) in mycorrhizas of Sp vs Tt, by N level and competition treatment


# Libraries needed:
library(cowplot)
library(emmeans)
library(tidyverse)
library(lme4)
library(lmerTest)
library(stargazer)

# Data:

bio_and_col_byplant = read_csv("processeddata/biomass_and_colonization_data_by_plant.csv")
bio_and_col_bycompt = read_csv("processeddata/biomass_and_colonization_data_by_compartment.csv")
isotopes = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly_INCLUDING_MIXED_and_pctCN.csv")
isotopes_nomixed = isotopes[-grep("MIXED", isotopes$Fungi),]


# Creating new data structure

masses = bio_and_col_bycompt %>%
  mutate(mass_Sp = (percent_col_Sp/100)*total_root_biomass_compartment,
         mass_Tt = (percent_col_Tt/100)*total_root_biomass_compartment,
         Side = tolower(Side)) %>%
  mutate(mass_NM_roots = total_root_biomass_compartment - mass_Sp - mass_Tt) %>%
  select(everything(), -percent_col_overall, -percent_col_Sp, -percent_col_Tt)

iso_justroots = isotopes_nomixed %>%
  select(Plant, Side, mycofungus, mycoC13ppmexcess, mycoN15ppmexcess, nmC13ppmexcess, nmN15ppmexcess, Batch)

total_compt_biomass = masses %>%
  select(Plant, Side, total_root_biomass_compartment) %>%
  spread(Side, total_root_biomass_compartment) %>%
  rename(total_root_biomass_compartment.a = a, total_root_biomass_compartment.b = b)

total_Tt_mass = masses %>%
  select(Plant, Side, mass_Tt) %>%
  spread(Side, mass_Tt) %>%
  rename(mass_Tt.a = a, mass_Tt.b = b)

total_Sp_mass = masses %>%
  select(Plant, Side, mass_Sp) %>%
  spread(Side, mass_Sp) %>%
  rename(mass_Sp.a = a, mass_Sp.b = b)

total_NM_mass = masses %>%
  select(Plant, Side, mass_NM_roots) %>%
  spread(Side, mass_NM_roots) %>%
  rename(mass_NM_roots.a = a, mass_NM_roots.b = b)

allmasses = left_join(total_Tt_mass, total_Sp_mass) %>%
  left_join(total_NM_mass)

iso_sp = iso_justroots %>%
  subset(mycofungus == "Sp")

C_sp = iso_sp %>%
  select(Plant, Side, mycoC13ppmexcess) %>%
  spread(Side, mycoC13ppmexcess) %>%
  rename(SpC13ppmexcess.a = a, SpC13ppmexcess.b = b)

N_sp = iso_sp %>%
  select(Plant, Side, mycoN15ppmexcess) %>%
  spread(Side, mycoN15ppmexcess) %>%
  rename(SpN15ppmexcess.a = a, SpN15ppmexcess.b = b)

iso_tt = iso_justroots %>%
  subset(mycofungus == "Tt")

C_tt = iso_tt %>%
  select(Plant, Side, mycoC13ppmexcess) %>%
  spread(Side, mycoC13ppmexcess) %>%
  rename(TtC13ppmexcess.a = a, TtC13ppmexcess.b = b)

N_tt = iso_tt %>%
  select(Plant, Side, mycoN15ppmexcess) %>%
  spread(Side, mycoN15ppmexcess) %>%
  rename(TtN15ppmexcess.a = a, TtN15ppmexcess.b = b)

C_NM = iso_justroots %>%
  select(Plant, Side, nmC13ppmexcess) %>%
  spread(Side, nmC13ppmexcess) %>%
  rename(nmC13ppmexcess.a = a, nmC13ppmexcess.b = b)

N_NM = iso_justroots %>%
  select(Plant, Side, nmN15ppmexcess) %>%
  spread(Side, nmN15ppmexcess) %>%
  rename(nmN15ppmexcess.a = a, nmN15ppmexcess.b = b)

allisos = left_join(C_sp, C_tt) %>%
  left_join(C_NM) %>%
  left_join(N_sp) %>%
  left_join(N_tt) %>%
  left_join(N_NM)





together = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly_INCLUDING_MIXED_and_pctCN.csv")

ndata = subset(together, received15N == "Y")
ndata = subset(ndata, compartment_fungus != "None")

ndata$versus2 = ndata$versus
ndata$versus3 = ndata$versus
for (i in 1:nrow(ndata)) {
  if ((ndata$versus[i] == "Mixed")|
      (ndata$compartment_fungus[i] == "MIXED")) {
    ndata$versus2[i] = "Other"
    ndata$versus3[i] = "Other"
  } else if ((ndata$versus[i] == "Tt" & ndata$compartment_fungus[i] == "Tt")|
             (ndata$versus[i] == "Sp" & ndata$compartment_fungus[i] == "Sp")) {
    ndata$versus2[i] = "Self"
    ndata$versus3[i] = "Self"
  } else if ((ndata$versus[i] == "Tt" & ndata$compartment_fungus[i] == "Sp")|
             (ndata$versus[i] == "Sp" & ndata$compartment_fungus[i] == "Tt")) {
    ndata$versus2[i] = "Other"
    ndata$versus3[i] = "Other"
  } else if (ndata$versus[i] == "None") {
    ndata$versus2[i] = "Other"
    ndata$versus3[i] = "None"
    
  }
}

# smallconstant = abs(min(ndata$nmN15ppmexcess)) + 1 # add 1 ppm plus minimum value
# # to make all n enrichment data log-able.
# 
# ndata$nmlogN15 = log(ndata$nmN15ppmexcess + smallconstant)
# ndata$mycologN15 = log(ndata$mycoN15ppmexcess + smallconstant)

mixedcompartments = ndata[ndata$compartment_fungus == "MIXED",]
allmixedcompartments = together[together$compartment_fungus == "MIXED",]

# nonm = together[!is.na(together$mycorrhizas.APE13C),]
# nonm = subset(nonm, compartment_fungus != "None")

# Exclude COMPARTMENTS with mixed cultures

excluding_mixed = ndata[-grep("Mixed", ndata$compartment_fungus),]

excluding_mixed$versus = as.factor(excluding_mixed$versus)
# relevel(excluding_mixed$versus, levels = c("None", "Sp", "Tt"))

justmycos = subset(ndata, mycofungus != "None")
justmycos_nomixed = justmycos[-grep("MIXED", justmycos$competitors),]

#### N-15 enrichment of mycos by species ####
# Does the N-15 enrichment of mycorrhizas depend on the species
# of fungus forming the mycorrhiza, controlling for competitor
# identity and N addition level?

# n15.myco.full = lmer(mycologN15 ~ mycofungus * N_level + (1|Batch/Plant), 
#                      data = justmycos) # I don't have any random effects here that I don't think I need

n15.myco.full = lmerTest::lmer(mycologN15 ~ compartment_fungus * N_level * versus + (1|Batch), 
                     data = justmycos_nomixed)

# class(n15.myco.full) = "lmerMod"

summary(n15.myco.full)

anovaresults = anova(n15.myco.full)

n15.myco.posthoc = emmeans(n15.myco.full, list(pairwise ~ compartment_fungus*N_level*versus), adjust = "tukey")
n15.myco.posthoc.interpretable = emmeans(n15.myco.full, list(pairwise ~ compartment_fungus*N_level), adjust = "tukey")


sink("stats_tables/N_by_fungus_competition_N_lme_results_mycos.html")
# sink("stats_tables/test.html")

stargazer(anovaresults, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()

sink("stats_tables/N_by_fungus_competition_N_mycos_anova_posthoc.txt")

n15.myco.posthoc.interpretable

sink()

### Including competition again ####
# n15.myco.withcomp = lmer(mycologN15 ~ mycofungus * N_level * versus3 + (1|Batch/Plant), 
#                      data = justmycos) # I don't have any random effects here that I don't think I need

n15.myco.withcomp = lmer(mycologN15 ~ mycofungus * N_level * versus3 + (1|Batch), 
                         data = justmycos_nomixed)

summary(n15.myco.withcomp)

anovaresults = anova(n15.myco.withcomp)

n15.myco.withcomp.posthoc = emmeans(n15.myco.withcomp, list(pairwise ~ mycofungus*N_level), adjust = "tukey")

sink("stats_tables/N_by_fungus_competition_N_lme_results_mycos.html")

stargazer(anovaresults, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()

sink("stats_tables/N_by_fungus_competition_N_mycos_anova_posthoc.txt")

n15.myco.withcomp.posthoc

sink()

# Examining results when you omit accidentally Tt plants
nocontam = read_csv("processeddata/Plants_with_no_Tt_contamination.csv")
justmycos_onlyclean = mutate(justmycos_nomixed, 
                             clean = justmycos_nomixed$Plant %in% nocontam$Plant) %>%
  filter(clean == TRUE)

n15.myco.withcomp.onlyclean = lmer(mycologN15 ~ mycofungus * N_level * versus3 + (1|Batch), 
                         data = justmycos_onlyclean) # rank deficient, dropping 3 columns

summary(n15.myco.withcomp.onlyclean)

anovaresults = anova(n15.myco.withcomp.onlyclean) 
# missing cells for: mycofungusTt:versus3None, N_levelHigh:versus3None, mycofungusSp:N_levelHigh:versus3None, mycofungusTt:N_levelHigh:versus3None, mycofungusTt:N_levelLow:versus3None.  
# Interpret type III hypotheses with care.

n15.myco.withcomp.posthoc.onlyclean = emmeans(n15.myco.withcomp.onlyclean, list(pairwise ~ mycofungus*N_level), adjust = "tukey")
# NOTE: Results may be misleading due to involvement in interactions
sink("stats_tables/N_by_fungus_competition_N_lme_results_mycos_onlyclean.html")

stargazer(anovaresults, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()

sink("stats_tables/N_by_fungus_competition_N_mycos_anova_posthoc_onlyclean.txt")

n15.myco.withcomp.posthoc.onlyclean

sink()




# Does the N-15 enrichment of NM roots depend on the species 
# of fungi in a given root compartment?
# n15.nm.full = lmer(nmlogN15 ~ compartment_fungus * N_level + (1|Batch), 
#                 data = ndata) # I don't have any random effects here that I don't think I need
# 
# summary(n15.nm.full)

# anovaresults = anova(n15.nm.full)
# 
# 
# sink("stats_tables/C_by_fungus_competition_N_lme_results_nmroots.html")
# 
# stargazer(anovaresults, type = "html",
#           digits = 3,
#           star.cutoffs = c(0.05, 0.01, 0.001),
#           digit.separator = "",
#           summary = FALSE,
#           no.space = TRUE)
# 
# sink()


#### Plot ####
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
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Labeled N in mycorrhizas\n(ln ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  geom_text(data = annotations, aes(x, y, label = labs))


nitrogencomparison_nm = ggplot(data = ndata) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = compartment_fungus, 
                   y = nmlogN15)) +
  geom_jitter(width = 0.20,
              aes(x = compartment_fungus, 
                  y = nmlogN15)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Labeled N in uncolonized\nfine roots (ln ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus")

save_plot("plots/N_comparison_byfungus_mycorrhizas.pdf",
          nitrogencomparison_mycos,
          base_aspect_ratio = 1.4)

save_plot("plots/N_comparison_byfungus_mycorrhizas.jpeg",
          nitrogencomparison_mycos,
          base_aspect_ratio = 1.4)

### Plot with competition ###

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
                     ystar = 7.7,
                     lab = ".")

labels = c(High = "High N", Low = "Low N")
nitrogencomparison_mycos_withcomp = ggplot(data = justmycos) +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(.9),
               aes(x = mycofungus, 
                   y = mycologN15,
                   fill = versus3)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15),
             aes(x = mycofungus, 
                 y = mycologN15,
                 fill = versus3)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Labeled N in mycorrhizas\n(ln ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  # geom_text(data = annotations, aes(x, y, label = labs)) +
  scale_fill_manual(values = c("lightgray", "gray42", "white")) +
  labs(fill = "Competitor") +
  geom_text(data = collabels, aes(x = xstar,  y = ystar, label = lab)) +
  geom_segment(data = collabels, aes(x = x1, xend = x2,
                                     y = y1, yend = y2),
               colour = "black") +
  geom_text(data = margsig, size = 10, aes(x = xstar,  y = ystar, label = lab)) +
  geom_segment(data = margsig, aes(x = x1, xend = x2,
                                     y = y1, yend = y2),
               colour = "black")

save_plot("plots/N_comparison_byfungus_mycorrhizas_withcomp.pdf",
          nitrogencomparison_mycos_withcomp,
          base_aspect_ratio = 1.4)
