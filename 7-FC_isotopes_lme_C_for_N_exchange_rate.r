# 7 - FC_isotopes_lme_C_for_N_exchange_rate
# June 2020
# Using an LME framework to more rigorously calculate how C for N 
# exchange rates varied in this experiment.

setwd("~/Documents/2018-2019/Fungal competition/fungal-competition2019/")

# Load required packages
require(cowplot)
require(tidyverse)
require(lmerTest)

# Loading required data
forN = read_csv("./FCdata/isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates.csv")
exchangerates = subset(forN, compartment_fungus != "MIXED")
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