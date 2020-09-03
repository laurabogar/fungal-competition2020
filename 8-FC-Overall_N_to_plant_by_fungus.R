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
# together = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly.csv")
together = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly_INCLUDING_MIXED.csv")

ndata = subset(together, received15N == "Y")
ndata = subset(ndata, compartment_fungus != "None")

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

excluding_mixed = ndata[-grep("MIXED", ndata$compartment_fungus),]
excluding_mixed$versus = as.factor(excluding_mixed$versus)
# relevel(excluding_mixed$versus, levels = c("None", "Sp", "Tt"))

justmycos = subset(ndata, mycofungus != "None")

#### N-15 enrichment of mycos by species ####
# Does the N-15 enrichment of mycorrhizas depend on the species
# of fungus forming the mycorrhiza, controlling for competitor
# identity and N addition level?

n15.myco.full = lmer(mycologN15 ~ mycofungus * N_level + (1|Batch/Plant), 
                     data = justmycos) # I don't have any random effects here that I don't think I need

summary(n15.myco.full)

anovaresults = anova(n15.myco.full)

n15.myco.posthoc = emmeans(n15.myco.full, list(pairwise ~ mycofungus*N_level), adjust = "tukey")

sink("stats_tables/N_by_fungus_competition_N_lme_results_mycos.html")

stargazer(anovaresults, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()

sink("stats_tables/N_by_fungus_competition_N_mycos_anova_posthoc.txt")

n15.myco.posthoc

sink()

tx = with(justmycos, interaction(mycofungus, N_level))
forlabels = aov(mycologN15~tx, data = justmycos)
mylabels = HSD.test(forlabels, "tx", group = TRUE)
# Dang, these don't match perfectly with my more 
# sophisticated post-hoc results. Will have to annotate
# manually.

annotations = data.frame(x = c((1:2), (1:2)),
                         y = c(4.5, 6.6, 5.6, 4.5),
                         N_level = c(rep("High", 2), rep("Low", 2)),
                         labs = c(paste(c("a", "b")), paste(c("ab", "a"))))


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

