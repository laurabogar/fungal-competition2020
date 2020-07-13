# Trying to analyze FC isotopes in a more cohesive way
# Attempting a linear mixed model a la Argüello et al. 2016
# August 2019

# This tutorial may be helpful: https://ourcodingclub.github.io/2017/03/15/mixed-models.html
# This is also good: http://www.bodowinter.com/tutorial/bw_LME_tutorial.pdf
# Hmm... https://stats.stackexchange.com/questions/35071/what-is-rank-deficiency-and-how-to-deal-with-it

# My model is consistently rank deficient.
# This seems unnecessarily annoying, given
# that the patterns I can see appear to be reasonably clear.

# Consider Bayesian approach? https://stats.stackexchange.com/questions/368272/why-does-logistic-regression-doesnt-give-the-same-result-as-a-chi-square-test-w

# THE PROCESS THAT HAS FINALLY WORKED:
# 1) Load lme4 and then lmerTest
# 2) Build full linear mixed effects model
# 3) Call "anova()" on that model to get a normal looking ANOVA table
# with significance labels.

setwd("~/Documents/2018-2019/Fungal competition/fungal-competition2019/")

# Libraries needed:
require(tidyverse)
require(lme4)

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
# identity?

# Argüello model: 
# log(14C in hyphal compartment A) ~ 
# AMF fungus side A*AMF fungus side B + 
# plant species identity + random pot effect
# with C labeling group "as a covariate." Does this mean as a random effect?

# Plant ID and C-13 labeling batch should  be random effects here.
# Ideally, I'd have a random slope for each,
# but a random intercept is all I can do with my dataset
# and is probably reasonable.

# Let's try an LRT

## trying with transformed data
c13.null = lmer(transmycoC13 ~ compartment_fungus +(1|Plant) + (1|Batch), 
                REML = FALSE,
                data = excluding_mixed)

plot(c13.null) # a little hump shaped, but maybe okay?
qqnorm(resid(c13.null)) # Looks reasonable to me.
qqline(resid(c13.null))

summary(c13.null)
# Note: I think I should be coding plant and batch as nested random effects,
# since they're certainly not crossed... Each plant is uniquely tied
# to a labeling batch.

c13.null2 = lmer(transmycoC13 ~ compartment_fungus +(1|Batch/Plant), 
                 REML = FALSE,
                 data = excluding_mixed)
plot(c13.null2) # looks identical to the plot for the crossed Batch + Plant scheme
qqnorm(resid(c13.null2)) # Also identical
qqline(resid(c13.null2))

summary(c13.null2) # This model gives info on Plant:Batch and Batch as 
# random effects, rather than just showing Plant and Batch separately
# as c13.null did. I think this is slightly more correct,
# although the model results are otherwise indistinguishable.

# Use ML to compare models, but report parameter estimates with REML.

c13.withversus = lmer(transmycoC13 ~ compartment_fungus + versus +(1|Batch/Plant), 
                      REML = FALSE,
                      data = excluding_mixed)
c13.withversusandint = lmer(transmycoC13 ~ compartment_fungus * versus +(1|Batch/Plant), 
                      REML = FALSE,
                      data = excluding_mixed)
c13.vsandintandn = lmer(transmycoC13 ~ compartment_fungus * versus + N_level +(1|Batch/Plant), 
                      REML = FALSE,
                      data = excluding_mixed)

c13.versusandnint = lmer(transmycoC13 ~ compartment_fungus * versus * N_level +(1|Batch/Plant), 
                      REML = FALSE,
                      data = excluding_mixed)
c13.versusandnint = lmer(transmycoC13 ~ compartment_fungus * versus * N_level + (1|Batch/Plant), 
                         REML = FALSE,
                         data = excluding_mixed)
# Okay actually even the above is working now.
# I think my consistent failures before were
# due to some weird problem with data filtering.
# Hooray?

anova(c13.null2, c13.withversus) # withversus is worse
anova(c13.withversus, c13.withversusandint) # withversusandint is worse
anova(c13.withversusandint, c13.vsandintandn) # withversusandint is worse,
# adding N in made the model better by AIC.

# Trying random slopes

c13.null.randslope = lmer(transmycoC13 ~ compartment_fungus +(compartment_fungus|Batch/Plant), 
                                      REML = FALSE,
                                      data = excluding_mixed)

# Error: number of observations (=108) <= number of random effects (=132) for term (compartment_fungus | Plant:Batch); the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
# I think this means I don't have enough data for such a sophisticated model.
# This seems fine to me.

c13.null = lmer(mycoC13ppmexcess ~ compartment_fungus + (mycoC13ppmexcess|Batch), 
                REML = FALSE,
                data = excluding_mixed)
# Error in eval_f(x, ...) : Downdated VtV is not positive definite

#### Zuur method ####

# Start with FULL ("beyond optimal," everything you can think of)
# Then sort out random effects structure using REML AIC/BIC
# Then do the fixed effects structure either using REML -- F or t statistic, or compare nexted MN models -- keep rand effects constant
# Then present final model using REML estimation

# My situation is kind of different, though,
# because I'm not exactly wondering WHAT factors
# might be important, here -- I have clear hypotheses
# that N level, fungal identity, and competitor identity
# should all matter. I think it makes sense, then,
# to present the full model and not bother with this weird
# top-down model selection.

# To be clear:
# Fungal identity should be significant
# Competitor identity should be significant
# N level should be significant
# *Interaction between fungal identity and N level MIGHT be significant
# *Interaction between fungal identity and competitor identity MIGHT be significant
# *Interaction between competitor identity and N level MIGHT be significant
# * Three way interaction between fungus and competitor and N level MIGHT be significant.

c13.full = lmer(transmycoC13 ~ compartment_fungus * versus * N_level + (1|Batch/Plant), 
                data = excluding_mixed) # I don't have any random effects here that I don't think I need

summary(c13.full)

require(stargazer)

stargazer(c13.full, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

c13.full_ml = lmer(transmycoC13 ~ compartment_fungus * versus * N_level + (1|Batch/Plant), 
                   REML = FALSE,
                   data = excluding_mixed) 
# I think reduced1 wins!
c13.reduced1 = lmer(transmycoC13 ~ compartment_fungus + versus + N_level + compartment_fungus:versus + compartment_fungus:N_level + versus:N_level + (1|Batch/Plant), 
                   REML = FALSE,
                   data = excluding_mixed)
c13.reduced2 = lmer(transmycoC13 ~ compartment_fungus + versus + N_level + compartment_fungus:versus + compartment_fungus:N_level + (1|Batch/Plant), 
                    REML = FALSE,
                    data = excluding_mixed) 
c13.reduced3 = lmer(transmycoC13 ~ compartment_fungus + versus + N_level + compartment_fungus:versus + versus:N_level +(1|Batch/Plant), 
                    REML = FALSE,
                    data = excluding_mixed)
c13.reduced4 = lmer(transmycoC13 ~ compartment_fungus + versus + N_level + compartment_fungus:N_level + versus:N_level +(1|Batch/Plant), 
                    REML = FALSE,
                    data = excluding_mixed)

# Likelihood ratio tests to determine factor significance
anova(c13.reduced1, c13.full_ml) # NS, reduced slightly better
# This might mean that the three way interaction is not significant.
anova(c13.reduced2, c13.reduced1) # marginal significance, reduced2 slightly worse
# I think this means that versus:N_level is marginally significant.
anova(c13.reduced3, c13.reduced1) # reduced 1 WAY better
# And this means that compartment_fungus:N_level is definitely significant.
anova(c13.reduced4, c13.reduced1) # NS, reduced4 slightly better.
# While compartment_fungus:versus is not significant.

require(lmerTest)

c13.full = lmer(transmycoC13 ~ compartment_fungus * versus * N_level + (1|Batch/Plant), 
                data = excluding_mixed) # I don't have any random effects here that I don't think I need


anova(c13.full)

# OH THANK GOD This is the tool I've been looking for this whole time.
# Some p values next to my factors.
# Done!!!



#### Bayesian? ####
require(blme)

glmb = blmer(transmycoC13 ~ compartment_fungus + (1|Batch), data=excluding_mixed, 
              fixef.prior = normal(cov = diag(9,3)))
summary(glmb)

glmb = blmer(transmycoC13 ~ compartment_fungus * versus * N_level + (1|Batch), data=excluding_mixed, 
             fixef.prior = normal(cov = diag(9,3)))
