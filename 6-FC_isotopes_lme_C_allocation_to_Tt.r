# 6 - FC_isotopes_lme_C_allocation_to_Tt
# June 2020
# Using an LME framework to more rigorously calculate how C allocation ratios
# to Tt shifted with competitive context

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

# Load required packages
require(cowplot)
require(tidyverse)
require(lmerTest)

# Loading required data
together = read_csv("FCdata/isotope_and_plant_metadata_with_competition_coded_clearly.csv")

# Declare constants:
smallconstant = 0.000473 # Need to add this to all C values to get around a negative value below.
artificial_minimum_alloc_ratio = 0.07 # minimum allocation ratio OBSERVED, besides 0, was 0.1428159.
# I'll add 0.07, half this value, to everything that is zero so it can be logged.

#### How does the allocation ratio to Tt change with competition? ####
# Hypothesis: Relative allocation to Tt should
# be greater when Tt is vs. Sp or vs. None.
# When vs. Tt, log alloc ratio should be zero.

allocratios = together[-grep("MIXED", together$Fungi),]
allocratios = allocratios[-grep("Mixed", allocratios$versus),]

allocratios$allocratio = numeric(nrow(allocratios))
allocratios$mycorrhizas.APE13C = allocratios$mycorrhizas.APE13C + smallconstant

for (i in 1:nrow(allocratios)) {
  for (j in 1:nrow(allocratios)) {
    if (allocratios$Plant[i] == allocratios$Plant[j] &
        allocratios$Side[i] != allocratios$Side[j]) {
      if (allocratios$compartment_fungus[i] == "Tt") {
        if (allocratios$compartment_fungus[j] != "None") {
          allocratios$allocratio[i] = allocratios$mycorrhizas.APE13C[i]/allocratios$mycorrhizas.APE13C[j]
          allocratios$allocratio[j] = NA # to ensure one entry per plant
        } else if (allocratios$compartment_fungus[j] == "None") {
          allocratios$allocratio[i] = allocratios$mycorrhizas.APE13C[i]/allocratios$uncolonized_roots.APE13C[j]
        }
      } else if (allocratios$compartment_fungus[i] == "Sp") {
        if (allocratios$compartment_fungus[j] == "Sp") {
          allocratios$allocratio[i] = allocratios$mycorrhizas.APE13C[i]/allocratios$mycorrhizas.APE13C[j]
          allocratios$allocratio[j] = NA # to ensure one entry per plant
        } else if (allocratios$compartment_fungus[j] == "None") {
          allocratios$allocratio[i] = allocratios$mycorrhizas.APE13C[i]/allocratios$uncolonized_roots.APE13C[j]
        }
      }
    }
  }
}

spalloc = allocratios[grep("Sp", allocratios$compartment_fungus),]
spalloc = spalloc[!is.na(spalloc$allocratio),]
spalloc$logallocratio = log(spalloc$allocratio + artificial_minimum_alloc_ratio)

allocratios = allocratios[grep("Tt", allocratios$compartment_fungus),]
allocratios = allocratios[!is.na(allocratios$allocratio),]

# allocratios$logallocratio = log(allocratios$allocratio + artificial_minimum_alloc_ratio) 
# whoops though! also have a negative value.
# It's a Tt/Tt plant. 6073. I guess at least
# one side wasn't really enriched for C?
# just6073 = together[together$Plant == 6073,]
# That value is -0.0004728348.
# Let's just add 0.000473 to everything upstream.


allocratios$logallocratio = log(allocratios$allocratio + artificial_minimum_alloc_ratio) 

labels = c(High = "High N", Low = "Low N")


### Statistical test with lmerTest ###
allocation_ratio.full = lmer(logallocratio ~ versus * N_level + (1|Batch), 
                data = allocratios)

# Boundary(singular) fit... but may not be wrong.

anova(allocation_ratio.full)

### A series of t tests to see if any of these distributions differ from null expectation ###

# I expect that, if the plant had NO preference, each of these distributions 
# would be centered around zero (the dashed line in my plot).
# Because of this expectation, I feel okay with using a t test to see 
# if any of these are actually different from zero. 

# I hypothesize that the log allocation ratio to Tt will be 
# significantly greater than zero when the plant prefers Tt,
# which in this experiment appears to be the High N conditions.

nonehigh = subset(allocratios, versus == "None" & N_level == "High")
t.test(nonehigh$logallocratio) # marginal: p = 0.06
otherhigh = subset(allocratios, versus == "Sp" & N_level == "High")
t.test(otherhigh$logallocratio) # marginal: p = 0.09
selfhigh = subset(allocratios, versus == "Tt" & N_level == "High")
t.test(selfhigh$logallocratio) # NS: p = 0.4899
nonelow = subset(allocratios, versus == "None" & N_level == "Low")
t.test(nonelow$logallocratio) # NS: p = 0.6373
otherlow = subset(allocratios, versus == "Sp" & N_level == "Low")
t.test(otherlow$logallocratio) # significant: p = 0.02295
selflow = subset(allocratios, versus == "Tt" & N_level == "Low")
t.test(selflow$logallocratio) # NS: p = 0.9918

### Plot ###

ggplot(data = allocratios) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = versus, y = logallocratio)) +
  geom_jitter(aes(x = versus, y = logallocratio),
              width = 0.25) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  ylab("Log C allocation ratio to Tt") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Competitor")

# ggplot(data = spalloc) +
#   geom_boxplot(outlier.alpha = 0,
#                aes(x = competition_treatment, y = logallocratio)) +
#   geom_jitter(aes(x = competition_treatment, y = logallocratio),
#               width = 0.25) +
#   geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
#   facet_grid(. ~ N_level, labeller = labeller(N_level = labels))

# VERY interesting.
# Tt gets more C relative to Sp under high N;
# basically same amount vs itself.
# Tt gets LESS N relative to Sp under low N,
# basically same amount vs itself.

justvsmycos = subset(allocratios, competition_treatment != "None")

# STATS for main text
t.test(logallocratio ~ competition_treatment, data = subset(justvsmycos, N_level == "High"))
t.test(logallocratio ~ competition_treatment, data = subset(justvsmycos, N_level == "Low"))

justvsmycos %>% group_by(N_level, competition_treatment) %>% summarize(mu = mean(logallocratio), stdev = sd(logallocratio), med = median(logallocratio))


anovaoverall = aov(logallocratio ~ competition_treatment*N_level, data = justvsmycos)
summary(anovaoverall)
# Significant interaction, not either thing by itself.
anovaflipped = aov(logallocratio ~ N_level*competition_treatment, data = justvsmycos)
summary(anovaflipped) # Basically same result, phew!

TukeyHSD(anovaoverall)
# Other low vs other high
# self high vs other high

figfortalk = ggplot(data = justvsmycos) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competition_treatment, y = logallocratio)) +
  geom_jitter(aes(x = competition_treatment, y = logallocratio),
              width = 0.25) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  xlab("Competition treatment") +
  ylab(bquote(atop("Log exchange rate", "("^13*"C to "^15*"N excess in mycorrhizas)")))

pdf("plots/Allocation_ratio_boxplot_fortalk.pdf", width = 8, height = 4)
figfortalk
dev.off()

median(justvsmycos$allocratio[justvsmycos$competition_treatment == "Other" &
                                justvsmycos$N_level == "High"])
# 1.904
median(justvsmycos$allocratio[justvsmycos$competition_treatment == "Self" &
                                justvsmycos$N_level == "High"])
# 1.0887

median(justvsmycos$allocratio[justvsmycos$competition_treatment == "Other" &
                                justvsmycos$N_level == "Low"])
# 0.4771 -- got 48% as much C as Suillus. 

median(justvsmycos$allocratio[justvsmycos$competition_treatment == "Self" &
                                justvsmycos$N_level == "Low"])
# 1.0814





starinfo = data.frame(x = c((1.5), (1.5)),
                      y = c(2.4, 2.1),
                      N_level = c(rep("High", 1), rep("Low", 1)),
                      labs = c(".", "*"))

seginfo = data.frame(x = c(0.9, 0.9), 
                     xend = c(2.1, 2.1),
                     y = c(2, 2), 
                     yend = c(2, 2), 
                     color = c("black", "black"),
                     N_level = c("High", "Low"))


figforpaper = ggplot(data = justvsmycos) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competition_treatment, y = logallocratio)) +
  geom_jitter(aes(x = competition_treatment, y = logallocratio),
              width = 0.25) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  geom_segment(data=seginfo,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE) +
  geom_text(data = starinfo, aes(x, y, label = labs), size = 8) +
  theme(plot.margin = unit(c(1,1,0.5,0.5), "cm")) +
  xlab("Competition treatment") +
  ylab(expression(atop("Log carbon allocation ratio","("*italic("Thelephora terrestris"*")"))))

save_plot("plots/MAIN_Allocation_ratio_for_paper.pdf",
          figforpaper,
          base_aspect_ratio = 1.8)

pdf("plots/Allocation_ratio_boxplot_withannots_fortalk.pdf", width = 8, height = 4)
figforpaper
dev.off()