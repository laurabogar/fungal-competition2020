# Trying to analyze FC isotopes in a more cohesive way
# This is NOT currently one of my central scripts 8/19/2020

# Attempting a linear mixed model a la Argüello et al. 2016
# August 2019/June 2020

# This tutorial was quite helpful: https://ourcodingclub.github.io/2017/03/15/mixed-models.html

# THE PROCESS THAT HAS FINALLY WORKED:
# 1) Load lme4 and then lmerTest
# 2) Build full linear mixed effects model
# 3) Call "anova()" on that model to get a normal looking ANOVA table
# with significance labels.

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

# Libraries needed:
library(cowplot)
library(tidyverse)
library(lme4)
library(lmerTest)

# Data:
together = read_csv("./FCdata/isotope_and_plant_metadata_with_competition_coded_clearly.csv")

# Shouldn't model mycorrhiza-specific phenomena with NM plants

nonm = together[!is.na(together$mycorrhizas.APE13C),]
nonm = subset(nonm, compartment_fungus != "None")
# Should I exclude microcosms with mixed cultures?

excluding_mixed = nonm[-grep("MIXED", nonm$competitors),]
# excluding_mixed$versus = relevel(excluding_mixed$versus, levels = c("None", "Sp", "Tt"))

#### C-13 enrichment of mycos by species ####
# Does the C-13 enrichment of mycorrhizas depend on the species
# of fungus forming the mycorrhiza, controlling for competitor
# identity and N addition level?

# Argüello model: 
# log(14C in hyphal compartment A) ~ 
# AMF fungus side A*AMF fungus side B + 
# plant species identity + random pot effect
# with C labeling group "as a covariate." Does this mean as a random effect?

# Plant ID and C-13 labeling batch should  be random effects here.
# Ideally, I'd have a random slope for each,
# but a random intercept is all I can do with my dataset
# and is probably reasonable.


# To be clear about my hypotheses:
# Fungal identity should be significant
# Competitor identity should be significant
# N level should be significant
# *Interaction between fungal identity and N level MIGHT be significant
# *Interaction between fungal identity and competitor identity MIGHT be significant
# *Interaction between competitor identity and N level MIGHT be significant
# * Three way interaction between fungus and competitor and N level MIGHT be significant.



c13.full = lmer(mycologC13 ~ compartment_fungus * versus * N_level + (1|Batch/Plant), 
                data = excluding_mixed) # I don't have any random effects here that I don't think I need

summary(c13.full)

anovaresults = anova(c13.full)


sink("stats_tables/C_by_fungus_competition_N_lme_results.txt")

anovaresults

sink()


# Significant factors: fungal identity (compartment fungus),
# N level, and the interaction between fungal identity and N level.
# Marginal: Interaction between competitor identity (versus) and N level
# Not significant: Competitor identity alone, 
# interaction between fungal identity and competitior identity (interesting!),
# three way interaction between fungal identity and competitor identity and N level.

#### Plot ####
labels = c(High = "High N", Low = "Low N")
ggplot(data = excluding_mixed) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = compartment_fungus, y = transmycoC13)) +
  geom_jitter(width = 0.20,
              aes(x = compartment_fungus, y = transmycoC13)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Plant C in mycorrhizas (log ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungal treatment")

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

library(stargazer)
nmvsmycos = lm(transmycoC13 ~ nmlogC13, data = excluding_mixed) # I don't have any random effects here that I don't think I need
stargazer(nmvsmycos, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

# Relationship: y = 0.741x + 1.497
# Adjusted R^2 = 0.366, p<0.001

# How well does 13C in mycos correspond to % N in those mycos?
# (Like Bogar et al. 2019?)

ggplot(data = excluding_mixed) +
  geom_point(aes(x = pctN, y = transmycoC13, shape = compartment_fungus, color = N_level)) +
  ylab("Plant C in mycorrhizas (log ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Percent N") +
  geom_smooth(method = lm, aes(x = pctN, y = transmycoC13))

CvspctN = lm(transmycoC13 ~ pctN, data = excluding_mixed)
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
