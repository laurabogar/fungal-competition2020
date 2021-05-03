# Trying to analyze FC isotopes in a more cohesive way
# Attempting a linear mixed model a la Argüello et al. 2016
# August 2019/June 2020

# This tutorial was quite helpful: https://ourcodingclub.github.io/2017/03/15/mixed-models.html

# THE PROCESS THAT HAS FINALLY WORKED:
# 1) Load lme4 and then lmerTest
# 2) Build full linear mixed effects model
# 3) Call "anova()" on that model to get a normal looking ANOVA table
# with significance labels.

setwd("~/Documents/2018-2019/Fungal competition/fungal-competition2019/")

# Libraries needed:
require(tidyverse)
require(cowplot)
require(lme4)
require(lmerTest)

# Data:
together = read_csv("./FCdata/isotopes_together_as_analyzed.csv")
# metadata = read_csv("./FCdata/Fungal_competition_plant_tracking.csv")
metadata_byplant = read_csv("./FCdata/percent_col_and_mass_data_by_plant.csv")

batchtomerge = select(metadata_byplant, Plant, Batch)

together = left_join(together, batchtomerge)

# Filtering

together$Batch = as.factor(together$Batch)

together = together[-grep("FAILED", together$competitors),]

together$versus = numeric(nrow(together))

for (i in 1:nrow(together)) {
  if (together$compartment_fungus[i] == "Sp") {
    if (together$competitors[i] == "SUIPU/NM") {
      together$versus[i] = "None"
    } else if (together$competitors[i] == "SUIPU/SUIPU") {
      together$versus[i] = "Sp"
    } else if (together$competitors[i] == "THETE/SUIPU") {
      together$versus[i] = "Tt"
    } else if (grepl("MIXED", together$competitors[i])){
      together$versus[i] = "Mixed"
    }
  } else if (together$compartment_fungus[i] == "Tt") {
    if (together$competitors[i] == "THETE/NM") {
      together$versus[i] = "None"
    } else if (together$competitors[i] == "THETE/THETE") {
      together$versus[i] = "Tt"
    } else if (together$competitors[i] == "THETE/SUIPU") {
      together$versus[i] = "Sp"
    } else if (grepl("MIXED", together$competitors[i])) {
      together$versus[i] = "Mixed"
    }
  } else if (together$compartment_fungus[i] == "None") {
    if (together$competitors[i] == "SUIPU/NM") {
      together$versus[i] = "Sp"
    } else if (together$competitors[i] == "NM/NM") {
      together$versus[i] = "None"
    } else if (together$competitors[i] == "THETE/NM") {
      together$versus[i] = "Tt"
    }
  } else if (together$compartment_fungus[i] == "MIXED") {
      if (together$competitors[i] == "MIXED/SUIPU") {
        together$versus[i] = "Sp"
      } else if (together$competitors[i] == "MIXED/THETE") {
        together$versus[i] = "Tt"
     }
  }
}

together$mycoC13ppmexcess = together$mycorrhizas.APE13C * (10^4)

together = together[together$enriched != 0,]

min(together$mycoC13ppmexcess[!is.na(together$mycoC13ppmexcess)])

together$transmycoC13 = (log(together$mycoC13ppmexcess))

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



c13.full = lmer(transmycoC13 ~ compartment_fungus * versus * N_level + (1|Batch/Plant), 
                data = excluding_mixed) # I don't have any random effects here that I don't think I need


anova(c13.full)

# OH THANK GOD This is the tool I've been looking for this whole time.
# Some p values next to my factors.
# Done!!!

# Significant factors: fungal identity (compartment fungus),
# N level, and the interaction between fungal identity and N level.
# Marginal: Interaction between competitor identity (versus) and N level
# Not significant: Competitor identity alone, 
# interaction between fungal identity and competitior identity (interesting!),
# three way interaction between fungal identity and competitor identity and N level.


labels = c(High = "High N", Low = "Low N")
ggplot(data = excluding_mixed) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = compartment_fungus, y = transmycoC13)) +
  geom_jitter(width = 0.20,
              aes(x = compartment_fungus, y = transmycoC13)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Plant C in mycorrhizas (log ppm excess)") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungal treatment")

### Carbon allocation ratio (Tt only) ###

justTt = subset(excluding_mixed, compartment_fungus == "Tt")
