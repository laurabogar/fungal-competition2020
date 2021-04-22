# 7 - FC_isotopes_lme_C_allocation_to_Tt
# June 2020
# Using an LME framework to more rigorously calculate how C allocation ratios
# to Tt shifted with competitive context

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

# Load required packages
library(agricolae)
library(cowplot)
library(emmeans)
library(tidyverse)
library(lmerTest)
library(stargazer)

# Loading required data
together = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly_INCLUDING_MIXED_and_pctCN.csv")


#### How does the allocation ratio to Tt change with competition? ####
# Hypothesis: Relative allocation to Tt should
# be greater when Tt is vs. Sp or vs. None.
# When vs. Tt, log alloc ratio should be zero.

allocratios = together[-grep("MIXED", together$Fungi),]
allocratios = allocratios[-grep("Mixed", allocratios$versus),]

allocratios$allocratio = rep(NA, nrow(allocratios))

# allocratios$mycorrhizas.APE13C = allocratios$mycorrhizas.APE13C + smallconstant

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

summary(allocratios$allocratio)

spalloc = allocratios[grep("Sp", allocratios$compartment_fungus),]
spalloc = spalloc[!is.na(spalloc$allocratio),]
spalloc$logallocratio = log(spalloc$allocratio)

allocratios = allocratios[grep("Tt", allocratios$compartment_fungus),]
allocratios = allocratios[!is.na(allocratios$allocratio),]

allocratios$logallocratio = log(allocratios$allocratio)

plantlist = allocratios %>% group_by(competitors, N_level) %>% summarize(list(Plant))
totalplantlist = together %>% group_by(competitors, N_level) %>% summarize(list(Plant))

unlist(plantlist[2,3])
unlist(totalplantlist[14,3])

unlist(plantlist[6,3])
# We lose plant 6030 between the full data and the allocation ratios. Why?
# Based on my notes, it had "lame mycos" (probably old and not numerous) at harvest.
# And I see no evidence of our ever having tin-balled that sample, or sent for analysis.
# It's unclear if the envelope was lost in the shuffle, or what happened, but it can't be included.
length(unlist(totalplantlist[18,3])) # okay, weird, plant 6072 shows up only once here.
# Ah! It's because UCSC lost the samples from side B of that root system,
# so I couldn't use it for pairwise comparison.


test = together %>% filter(Plant == 6072)

nomixed = together %>% filter(!str_detect(versus, "Mixed"), 
                              !str_detect(competitors, "MIXED"), 
                              str_detect(compartment_fungus, "Tt"))
nomixed %>% group_by(competitors, N_level) %>% summarize(n())
### Statistical test with lmerTest ###
allocation_ratio.full = lmer(logallocratio ~ versus * N_level + (1|Batch), 
                data = allocratios)


allocanova = anova(allocation_ratio.full)
allocposthoc = emmeans(allocation_ratio.full, list(pairwise ~ versus*N_level), adjust = "tukey")

sink("stats_tables/Relative_C_allocation_anova.html")

stargazer(allocanova, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()

sink("stats_tables/Relative_C_allocation_anova_posthoc.txt")

allocposthoc

sink()

tx = with(allocratios, interaction(versus, N_level))
forlabels = aov(logallocratio~tx, data = allocratios)
mylabels = HSD.test(forlabels, "tx", group = TRUE)
# although the agricolae labels use a less sophisticated
# anova, the suggested letters match what I would
# need to describe my post-hoc emmeans results. Much
# easier to double check this outcome than try to 
# come up with the significance labels on my own!

### Plot ###
labels = c(High = "High N", Low = "Low N")
annotations = data.frame(x = c((1:3), (1:3)),
                         y = c(2.1, 2.1, 2.1, 1.6, 1.1, 2.1),
                         N_level = c(rep("High", 3), rep("Low", 3)),
                         labs = c(paste(c("a", "a", "a")), paste(c("ab", "b", "a"))))

myplot = ggplot(data = allocratios) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = versus, y = logallocratio)) +
  geom_jitter(aes(x = versus, y = logallocratio),
              width = 0.25) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  ylab("Log C allocation ratio to Tt") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Competitor") +
  scale_x_discrete(labels = c("None" = "None", 
                              "Sp" = "Other",
                              "Tt" = "Self"))+
  geom_text(data = annotations, aes(x, y, label = labs))

save_plot("plots/Relative_C_allocation_to_Tt.pdf", 
          myplot,
          ncol = 1,
          base_aspect_ratio = 1.8)
save_plot("plots/Relative_C_allocation_to_Tt.jpeg", 
          myplot,
          ncol = 1,
          base_aspect_ratio = 1.8)

# VERY interesting.
# Tt gets more C relative to Sp under high N;
# basically same amount vs itself.
# Tt gets LESS N relative to Sp under low N,
# basically same amount vs itself.

median(allocratios$logallocratio[allocratios$versus == "Sp" & allocratios$N_level == "High"]) # 0.658
exp(0.658) #1.93 -- Tt got about 2x the carbon

median(allocratios$logallocratio[allocratios$versus == "Sp" & allocratios$N_level == "Low"]) # - 0.7733
exp(-0.773) # 0.46 -- Tt got about half the carbon.


### T tests ###
#None high p = 0.1
t.test(allocratios$logallocratio[allocratios$versus == "None" & allocratios$N_level == "High"], alternative = "greater")
# Sp high p = 0.1234 Dang it!
t.test(allocratios$logallocratio[allocratios$versus == "Sp" & allocratios$N_level == "High"], alternative = "greater")
# Tt high p = 0.9927
t.test(allocratios$logallocratio[allocratios$versus == "Tt" & allocratios$N_level == "High"])


#### OLD AND IRRELEVANT FROM HERE ON. ####

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