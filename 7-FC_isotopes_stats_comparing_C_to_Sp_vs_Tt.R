# Fungal competition: Code and stats associated with Figure 3a

# Comparing plant C allocated to Sp vs Tt, by N level and competition treatment



# This tutorial was quite helpful: https://ourcodingclub.github.io/2017/03/15/mixed-models.html

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

# Libraries needed:
library(tidyverse)
library(cowplot)
library(lme4)
library(lmerTest)
library(emmeans)
library(stargazer)

# Data:
together = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly_INCLUDING_MIXED_and_pctCN.csv")

# Shouldn't model mycorrhiza-specific phenomena with NM plants

nonm = together[!is.na(together$mycorrhizas.APE13C),]
nonm = subset(nonm, compartment_fungus != "None")
nonm$versus2 = nonm$versus
nonm$versus3 = nonm$versus
for (i in 1:nrow(nonm)) {
  if (nonm$versus[i] == "Mixed") {
    nonm$versus2[i] = "Other"
    nonm$versus3[i] = "Other"
  } else if ((nonm$versus[i] == "Tt" & nonm$mycofungus[i] == "Tt")|
        (nonm$versus[i] == "Sp" & nonm$mycofungus[i] == "Sp")) {
          nonm$versus2[i] = "Self"
          nonm$versus3[i] = "Self"
  } else if ((nonm$versus[i] == "Tt" & nonm$mycofungus[i] == "Sp")|
             (nonm$versus[i] == "Sp" & nonm$mycofungus[i] == "Tt")) {
    nonm$versus2[i] = "Other"
    nonm$versus3[i] = "Other"
  } else if (nonm$versus[i] == "None") {
    nonm$versus2[i] = "Other"
    nonm$versus3[i] = "None"
    
  }
}

write_csv(nonm, "processeddata/isotope_and_plant_metadata_with_competition_coded_clearly_INCLUDING_MIXED_betterversus.csv")

# Should I exclude microcosms with mixed cultures?
nonm = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly_INCLUDING_MIXED_betterversus.csv")
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

# c13.full = lmer(mycologC13 ~ mycofungus * versus2 * N_level + (1|Batch/Plant), 
#                 data = nonm) 


c13.full = lmer(mycologC13 ~ compartment_fungus * versus * N_level + (1|Batch/Plant),
                data = excluding_mixed) # I don't have any random effects here that I don't think I need

summary(c13.full)

anovaresults = anova(c13.full)


sink("stats_tables/C_by_fungus_competition_N_lme_anova.html")

stargazer(anovaresults, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()

Cbyfungusposthoc = emmeans(c13.full, list(pairwise ~ compartment_fungus*N_level), adjust = "tukey")

Cbyfungusposthoc_withversus = emmeans(c13.full, list(pairwise ~ compartment_fungus*N_level*versus), adjust = "tukey")

sink("stats_tables/C_by_fungus_competition_N_lme_anova_posthoc_withversus.txt")

Cbyfungusposthoc_withversus
sink()

sink("stats_tables/C_by_fungus_competition_N_lme_anova_posthoc.txt")

Cbyfungusposthoc
sink()

# Trying this with only plants that were never accidentally inoculated with Tt

nocontam = read_csv("processeddata/Plants_with_no_Tt_contamination.csv")
excluding_mixed_onlyclean = mutate(excluding_mixed, clean = excluding_mixed$Plant %in% nocontam$Plant)
excluding_mixed_onlyclean = filter(excluding_mixed_onlyclean, clean == TRUE)

c13.full.onlyclean = lmer(mycologC13 ~ compartment_fungus * versus * N_level + (1|Batch/Plant),
                data = excluding_mixed_onlyclean) # I don't have any random effects here that I don't think I need

summary(c13.full.onlyclean)

anovaresults = anova(c13.full.onlyclean)


sink("stats_tables/C_by_fungus_competition_N_lme_anova_onlyclean.html")

stargazer(anovaresults, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()

Cbyfungusposthoc = emmeans(c13.full.onlyclean, list(pairwise ~ compartment_fungus*N_level), adjust = "tukey")

Cbyfungusposthoc_withversus = emmeans(c13.full.onlyclean, list(pairwise ~ compartment_fungus*N_level*versus), adjust = "tukey")

sink("stats_tables/C_by_fungus_competition_N_lme_anova_posthoc_withversus_onlyclean.txt")

Cbyfungusposthoc_withversus
sink()

sink("stats_tables/C_by_fungus_competition_N_lme_anova_posthoc_onlyclean.txt")

Cbyfungusposthoc
sink()


# sink("stats_tables/C_by_fungus_competition_N_lme_results_noanova.html")
# 
# stargazer(c13.full, type = "html",
#           digits = 3,
#           star.cutoffs = c(0.05, 0.01, 0.001),
#           digit.separator = "",
#           summary = TRUE,
#           no.space = TRUE)
# 
# sink()


# Significant factors: fungal identity (compartment fungus),
# N level, and the interaction between fungal identity and N level.
# Marginal: Interaction between competitor identity (versus) and N level
# Not significant: Competitor identity alone, 
# interaction between fungal identity and competitior identity (interesting!),
# three way interaction between fungal identity and competitor identity and N level.

#### Plot ####
labels = c(High = "High N", Low = "Low N")

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
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Plant C in mycorrhizas\n(ln ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  geom_text(data = annotations, aes(x, y, label = labs))


save_plot("plots/Sp_gets_more_C_than_Tt_in_low_N.pdf",
          carboncomparison,
          base_aspect_ratio = 1.4)

save_plot("plots/Sp_gets_more_C_than_Tt_in_low_N.jpeg",
          carboncomparison,
          base_aspect_ratio = 1.4)

carboncomparison_withcompetition = ggplot(data = excluding_mixed) +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(.9),
               aes(x = compartment_fungus, 
                   y = mycologC13,
                   fill = versus3)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15),
              aes(x = compartment_fungus, 
                  y = mycologC13,
                  fill = versus3)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Plant C in mycorrhizas\n(ln ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  # geom_text(data = annotations, aes(x, y, label = labs))
  scale_fill_manual(values = c("lightgray", "gray42", "white")) +
  labs(fill = "Competitor")
  

save_plot("plots/Sp_gets_more_C_than_Tt_in_low_N.pdf",
          carboncomparison,
          base_aspect_ratio = 1.4)



carboncomparison_includingmixed = ggplot(data = nonm) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = mycofungus, y = mycologC13)) +
  geom_jitter(width = 0.20,
              aes(x = mycofungus, 
                  y = mycologC13)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Plant C in mycorrhizas\n(ln ppm excess)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungus") +
  geom_text(data = annotations, aes(x, y, label = labs))


save_plot("plots/Sp_gets_more_C_than_Tt_in_low_N.pdf",
          carboncomparison,
          base_aspect_ratio = 1.4)
#### EXTRA ANALYSES NOT PRESENTED IN MAIN TEXT FOLLOW ####
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
