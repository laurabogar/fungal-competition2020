# 8 - overall N to plant, Tt vs Sp

# Trying to analyze FC isotopes in a more cohesive way
# This is NOT currently one of my central scripts 8/19/2020

# I have cool evidence of functional complementarity here, but
# I need to rearrange my data frame to show this most effectively.
# Probably best to leave this one for last, once the other analyses
# are looking good.

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

# Libraries needed:
library(cowplot)
library(tidyverse)
library(lme4)
library(lmerTest)
library(stargazer)

# Data:
# together = read_csv("./FCdata/isotope_and_plant_metadata_with_competition_coded_clearly.csv")
# together = read_csv("processeddata/isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates.csv")
together = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly.csv")


ndata = subset(together, received15N == "Y")
ndata = subset(ndata, compartment_fungus != "None")

smallconstant = abs(min(ndata$nmN15ppmexcess)) + 1 # add 1 ppm plus minimum value
# to make all n enrichment data log-able.

ndata$nmlogN15 = log(ndata$nmN15ppmexcess + smallconstant)
ndata$mycologN15 = log(ndata$mycoN15ppmexcess + smallconstant)

mixedcompartments = ndata[ndata$compartment_fungus == "MIXED",]
allmixedcompartments = together[together$compartment_fungus == "MIXED",]

# nonm = together[!is.na(together$mycorrhizas.APE13C),]
# nonm = subset(nonm, compartment_fungus != "None")

# Exclude COMPARTMENTS with mixed cultures

excluding_mixed = ndata[-grep("MIXED", ndata$compartment_fungus),]
excluding_mixed$versus = as.factor(excluding_mixed$versus)
# relevel(excluding_mixed$versus, levels = c("None", "Sp", "Tt"))

#### N-15 enrichment of mycos by species ####
# Does the N-15 enrichment of mycorrhizas depend on the species
# of fungus forming the mycorrhiza, controlling for competitor
# identity and N addition level?

n15.nm.full = lmer(nmlogN15 ~ compartment_fungus * N_level + (1|Batch), 
                data = excluding_mixed) # I don't have any random effects here that I don't think I need

summary(n15.nm.full)

anovaresults = anova(n15.nm.full)


sink("stats_tables/C_by_fungus_competition_N_lme_results_nmroots.html")

stargazer(anovaresults, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()

n15.myco.full = lmer(mycologN15 ~ compartment_fungus * N_level + (1|Batch), 
                   data = excluding_mixed) # I don't have any random effects here that I don't think I need

summary(n15.myco.full)

anovaresults = anova(n15.myco.full)


sink("stats_tables/C_by_fungus_competition_N_lme_results_mycos.html")

stargazer(anovaresults, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()


#### Plot ####
labels = c(High = "High N", Low = "Low N")
nitrogencomparison_mycos = ggplot(data = excluding_mixed) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = compartment_fungus, 
                   y = mycologN15)) +
  geom_jitter(width = 0.20,
              aes(x = compartment_fungus, 
                  y = mycologN15,
                  shape = mycofungus)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Labeled N in mycorrhizas\n(ln ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus")

nitrogencomparison_nm = ggplot(data = excluding_mixed) +
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

save_plot("plots/N_comparison_byfungus_nmroots.pdf",
          nitrogencomparison_nm)

# How well does 13C in mycos correspond to 13C in adjacent fine roots?

ggplot(data = excluding_mixed) +
  geom_point(aes(x = nmlogC13, y = mycologC13, shape = compartment_fungus, color = N_level)) +
  ylab("Plant C in mycorrhizas (log ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Plant C in NM roots (log ppm excess)") +
  geom_smooth(method = lm, aes(x = nmlogC13, y = mycologC13))

nmvsmycos = lmer(mycologC13 ~ nmlogC13 * versus * N_level + (1|Batch/Plant), 
                 data = excluding_mixed) # I don't have any random effects here that I don't think I need

summary(nmvsmycos)
anova(nmvsmycos) # I do NOT know how to understand this in
# the context of my plot. Let's go simpler.

nmvsmycos = lm(mycologC13 ~ nmlogC13, data = excluding_mixed) # I don't have any random effects here that I don't think I need
stargazer(nmvsmycos, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

# Relationship: y = 0.741x + 1.497
# Adjusted R^2 = 0.366, p<0.001

# How well does 13C in mycos correspond to % N in those mycos?
# (Like Bogar et al. 2019?)

ggplot(data = excluding_mixed) +
  geom_point(aes(x = pctN, y = mycologC13, shape = compartment_fungus, color = N_level)) +
  ylab("Plant C in mycorrhizas (log ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Percent N") +
  geom_smooth(method = lm, aes(x = pctN, y = mycologC13))

CvspctN = lm(mycologC13 ~ pctN, data = excluding_mixed)
stargazer(CvspctN, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

# I'm honestly not confident that the N label really... worked. These
# values are so close to zero.
ggplot(data = excluding_mixed) +
  geom_point(aes(x = pctN, y = log(mycorrhizas.APE15N*10^4), shape = compartment_fungus, color = N_level)) +
  ylab("Labeled N in mycorrhizas\n(log ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Percent N") +
  geom_smooth(method = lm, aes(x = pctN, y = log(mycorrhizas.APE15N*10^4)))

# What if I included NM plants to see if mycos are different from fine roots?
everything_nomixed = together[-grep("MIXED", together$competitors),]

ggplot(data = everything_nomixed) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = compartment_fungus, y = nmlogC13)) +
  geom_jitter(width = 0.20,
              aes(x = compartment_fungus, y = nmlogC13)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Plant C in NM roots (log ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungal treatment")
