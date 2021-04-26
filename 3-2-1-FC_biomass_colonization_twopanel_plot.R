# 3-2-1 Combining biomass and colonization plots into one figure

# significance labels determined by anova for interaction between N level and fungal treatment.

# setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(agricolae)
library(effects)
library(emmeans)
library(cowplot)
library(lme4)
library(lmerTest)
library(MuMIn)
library(tidyverse)
library(stargazer)

alldata = read_csv("processeddata/percent_col_and_mass_data_by_plant.csv")
compdata = read_csv("processeddata/percent_colonization_and_mass_data_by_compartment.csv")
granular_data_bycompt = read_csv("processeddata/granular_mass_and_colonization_data_by_compartment.csv")
# Note: The below file is produced when you run this script; you could also just start with "colfortest" below.
colfortest = read_csv("processeddata/colonization_and_biomass_data_by_compartment_with_competition.csv")


#### BIOMASS ####

## Stats ##
biomass_lme = lmer(total_biomass ~ N_level * Fungi * percent_col + (1|Batch), data = alldata)
anovaresults = anova(biomass_lme)
anovaresults # Hmm. I should have anticipated: Type III ANOVA table shows 
# only N level having a significant impact on biomass,
# probably because I have such poor replication for Suillus.
summary(biomass_lme)

anotherway = lmer(total_biomass ~ N_level * Fungi + (1|percent_col), data = alldata)
anotheranova = anova(anotherway) # okay, now everything is HIGHLY significant
summary(anotherway)

justNandF = lmer(total_biomass ~ N_level * Fungi + (1|Batch), data = alldata)
anovaagain = anova(justNandF) # super significant again.

# Having compared the models WITH percent colonization
# to those without, I think it looks like Tt plants were
# bigger BECAUSE they had much higher % colonization
# (as far as the model is concerned)

# Would it make more sense to model biomass as a function
# of %Tt, %Sp, and N_level, instead?

# Let's add percent colonization info by fungus to the whole dataset

bio_formodel = granular_data_bycompt[-grep("MIXED", granular_data_bycompt$competitors),]
bio_formodel = subset(bio_formodel, competitors != "FAILED" &
                      competitors != "THETE" &
                      N_level != "None")

bio_and_col = bio_formodel %>%
  group_by(Plant, Side, competitors,
           compartment_fungus, N_level,
           total_root_biomass_compartment) %>%
  summarize(percent_col_overall = 
              (100*(sum(Sp_myco_mass, Tt_myco_mass, na.rm = TRUE))/sum(uncolonized_root_mass, Sp_myco_mass, Tt_myco_mass, na.rm = TRUE)),
            percent_col_Sp = 
              (100*(sum(Sp_myco_mass, na.rm = TRUE))/sum(uncolonized_root_mass, Sp_myco_mass, Tt_myco_mass, na.rm = TRUE)),
            percent_col_Tt = 
              (100*(sum(Tt_myco_mass, na.rm = TRUE))/sum(uncolonized_root_mass, Sp_myco_mass, Tt_myco_mass, na.rm = TRUE)))

bio_and_col_byplant = bio_formodel %>%
  group_by(Plant, competitors, N_level, competitors_attempted) %>%
  summarize(total_root_mass = sum(total_root_biomass_compartment),
            percent_col_overall = 
              (100*(sum(Sp_myco_mass, Tt_myco_mass, na.rm = TRUE))/sum(uncolonized_root_mass, Sp_myco_mass, Tt_myco_mass, na.rm = TRUE)),
            percent_col_Sp = 
              (100*(sum(Sp_myco_mass, na.rm = TRUE))/sum(uncolonized_root_mass, Sp_myco_mass, Tt_myco_mass, na.rm = TRUE)),
            percent_col_Tt = 
              (100*(sum(Tt_myco_mass, na.rm = TRUE))/sum(uncolonized_root_mass, Sp_myco_mass, Tt_myco_mass, na.rm = TRUE)))

justmass = select(alldata, Plant, total_biomass, Batch)

bio_and_col_byplant = left_join(bio_and_col_byplant, justmass)

# Stats with the linear framework:

biomass_by_colonization = lmer(total_biomass ~ N_level * percent_col_Sp * percent_col_Tt + (1|Batch), data = bio_and_col_byplant)
summary(biomass_by_colonization)
biomass_rsquared = round(r.squaredGLMM(biomass_by_colonization)[2], 4)
# Okay, here we have N level as a highly significant predictor
# of biomass, and percent colonization by Tt as marginally significant
# This seems pretty reasonable.

biomass_simple = lm(total_biomass ~ N_level * percent_col_Sp * percent_col_Tt, data = bio_and_col_byplant)
summary(biomass_simple)


# Is it reasonable for me to use the random effect?
# Found this suggestion on Stack Overflow:
biomass_lmer = lmer(total_biomass ~ N_level * percent_col_Sp * percent_col_Tt + (1|Batch), 
                    data = bio_and_col_byplant, REML = FALSE)

AIC(biomass_lmer, biomass_simple)
# Lower AIC for the lmer model, let's do that.



# Now, let's see if being a contaminated plant makes a differences
# for biomass OR colonization.

bio_and_col_byplant$uncontaminated = bio_and_col_byplant$competitors_attempted == bio_and_col_byplant$competitors
bio_and_col_onlyintended = subset(bio_and_col_byplant, uncontaminated == TRUE)

biomass_by_colonization_contamtest = lmer(total_biomass ~ N_level * 
                                 percent_col_Sp * percent_col_Tt *
                                 uncontaminated + 
                                 (1|Batch), data = bio_and_col_byplant)
summary(biomass_by_colonization_contamtest)
anova(biomass_by_colonization_contamtest)

# Well, it looks like the contaminated ones are indeed different.
# Specifically, we've got significant interactions between percent TT and 
# contamination, as well as with N level.

### This is not entirely correct: I only want to know which ones got Tt that
# weren't supposed to... but I think all of these were those actually.

# could I define this contamination a little more specifically?
# Contamination by Tt represents only the following transitions:
# None/None -> anything with Tt
# Tt/None --> Tt/Tt
# Tt/Sp --> Tt/Tt
# Sp/Sp --> anything with Tt
# Sp/None --> anything with Tt

grepl("Tt/Tt", bio_and_col_byplant$competitors_attempted[1])

bio_and_col_byplant$contaminated = numeric(nrow(bio_and_col_byplant))
bio_and_col_byplant$contaminated = FALSE

for (i in 1:nrow(bio_and_col_byplant)) {
  if (bio_and_col_byplant$competitors_attempted[i] == "None/None" |
      bio_and_col_byplant$competitors_attempted[i] == "Sp/None" |
      bio_and_col_byplant$competitors_attempted[i] == "Sp/Sp") {
    bio_and_col_byplant$contaminated[i] = grepl("Tt", bio_and_col_byplant$competitors[i])
  } else if (bio_and_col_byplant$competitors_attempted[i] == "Tt/None" |
             bio_and_col_byplant$competitors_attempted[i] == "Tt/Sp") {
    bio_and_col_byplant$contaminated[i] = grepl("Tt/Tt", bio_and_col_byplant$competitors[i])
    }
  }

bio_and_col_onlyclean = subset(bio_and_col_byplant, contaminated == FALSE)
plantlist = select(bio_and_col_onlyclean, Plant)
write_csv(plantlist, "processeddata/Plants_with_no_Tt_contamination.csv")
write_csv(bio_and_col_onlyclean, "processeddata/Biomass_and_colonization_forplants_with_no_Tt_contamination.csv")


biomass_by_colonization_onlyclean = lmer(total_biomass ~ N_level * 
                                            percent_col_Sp * percent_col_Tt +
                                            (1|Batch), data = bio_and_col_onlyclean)
summary(biomass_by_colonization_onlyclean)
biomass_rsquared = round(r.squaredGLMM(biomass_by_colonization_onlyclean)[2], 4)

plot(allEffects(biomass_by_colonization_onlyclean))
plot(predictorEffect("percent_col_Sp", biomass_by_colonization_onlyclean),
     lines = list(multiline = FALSE))


class(biomass_by_colonization_onlyclean) = "lmerMod"

sink("stats_tables/biomass_lmer_noTtcontam.html")

stargazer(biomass_by_colonization_onlyclean, type = "html",
          dep.var.labels = "Total plant biomass (g)",
          covariate.labels = c("N level (low)",
                               "Percent root mass colonized by Sp",
                               "Percent root mass colonized by Tt",
                               "N level (low):Percent Sp",
                               "N level (low):Percent Tt",
                               "Percent Sp:Percent Tt",
                               "N level (low):Percent Sp:Percent Tt"),
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          star.char = c("*", "**", "***"),
          digit.separator = "",
          summary = TRUE,
          no.space = TRUE,
          add.lines = list(c("Conditional pseudo-$R2$",
                             biomass_rsquared)),
          ci = TRUE,
          notes.append=FALSE)

sink()

### Plot ###

# This function  does automatic letters
tx = with(alldata, interaction(N_level, Fungi))
anovaforplot = aov(total_biomass ~ tx, data = alldata)

# from "agricolae" package
mylabels = HSD.test(anovaforplot, "tx", group = TRUE)

anothertry = data.frame(x = c((1:6), (1:6)),
                        y = c(4, 4, 4, 6, 6, 7.5, 4, 4, 4, 4, 4, 4),
                        N_level = c(rep("High", 6), rep("Low", 6)),
                        labs = c(paste(c("bc", "c", "bc", "b", "b", "a")), paste(c("c", "c", "c", "c", "c", "c"))))

labels = c(High = "High N", Low = "Low N")
massplot = ggplot(data = alldata) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = Fungi, y = total_biomass)) +
  geom_jitter(width = 0.20,
              aes(x = Fungi, y = total_biomass)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Total plant biomass (g)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_text(data = anothertry, aes(x, y, label = labs)) +
  xlab("Fungi on roots at harvest")

massplot_noletters = ggplot(data = bio_and_col_onlyclean) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competitors, y = total_biomass)) +
  geom_jitter(width = 0.20,
              aes(x = competitors, y = total_biomass)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Total plant biomass (g)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungal treatment")

save_plot("plots/Biomass_boxplot_noletters.pdf",
          massplot_noletters,
          base_aspect_ratio = 1.8)

save_plot("plots/Biomass_boxplot_noletters.png",
          massplot_noletters,
          base_aspect_ratio = 1.8)

pdf("plots/Biomass_boxplot.pdf", width = 9, height = 5)
massplot
dev.off()

ggsave("plots/Biomass_boxplot.jpeg", plot = massplot,
       device = "jpeg",
       width = 9, height = 5, units = "in")


#### COLONIZATION ####

### Analyses ###

# colforplot = subset(alldata, Fungi != "None/None")
# remove the "NA" values because these are all failed splits
# except 6105, which is a plant for which I don't have colonization
# data -- I think I harvested it before I had developed that protocol? 
# Either that or I lost the relevant tissue envelopes.
colforplot = granular_data_bycompt[-grep("MIXED", granular_data_bycompt$competitors),]
colforplot = subset(colforplot, competitors != "FAILED" &
                      competitors != "THETE" &
                      N_level != "None")

percent_col = colforplot %>%
  group_by(Plant, Side, competitors,
           compartment_fungus, N_level) %>%
  summarize(percent_col = (100*(sum(Sp_myco_mass, Tt_myco_mass, na.rm = TRUE))/sum(uncolonized_root_mass, Sp_myco_mass, Tt_myco_mass, na.rm = TRUE)))

colforplot = percent_col

# colforplot = full_join(colforplot, percent_col) %>% distinct()

# 
# colforplot = compdata[!is.na(compdata$compartment_fungus),]
# colforplot = subset(colforplot, compartment_fungus != "MIXED" &
#                       compartment_fungus != "OTHER" &
#                       competitors != "THETE" &
#                       competitors != "FAILED")
# colforplot = colforplot[-grep("MIXED", colforplot$competitors),]
# colforplot = subset(colforplot, N_level != "None" &
#                       compartment_fungus != "Failed")
# colforplot$compartment_fungus = recode(colforplot$compartment_fungus,
#                                         "NM" = "None",
#                                         "SUIPU" = "Sp",
#                                         "THETE" = "Tt")
# colforplot$bettercomp = fct_relevel(colforplot$attempted,
#                                    "None/None",
#                                    "Sp/None",
#                                    "Sp/Sp",
#                                    "Tt/Sp",
#                                    "Tt/None",
#                                    "Tt/Tt")
colforplot$competitors_reordered = fct_relevel(colforplot$competitors,
                                    "None/None",
                                    "Sp/None",
                                    "Sp/Sp",
                                    "Tt/Sp",
                                    "Tt/None",
                                    "Tt/Tt")

# colforplot$Fungus_attempted[colforplot$Plant == 6106 & colforplot$Side == "A"] = "Tt" # why not in data? unclear
# colforplot$Fungus_attempted[colforplot$Plant == 6106 & colforplot$Side == "B"] = "Sp"

# dead_mycos = colforplot %>% group_by(N_level, Fungus_attempted, mycofungus) %>% 
#   filter(dead_tissue > 0) %>%
#   summarize(n(), sum(dead_tissue))

# If including mixed:
# dead_mycos = colforplot %>% group_by(N_level, Fungus_attempted, compartment_fungus) %>% 
#   filter(dead_tissue > 0) %>%
#   summarize(n(), sum(dead_tissue))

# what happens if we use only plants that were NOT contaminated by errant Thelephora?

col_onlyclean = subset(colforplot, Plant %in% bio_and_col_onlyclean$Plant)

col_onlyintended = subset(colforplot, Plant %in% bio_and_col_onlyintended$Plant) # I believe this is the same as above.


labels = c(High = "High N", Low = "Low N")

tx = with(colforplot, interaction(N_level, Fungi))
anovaforplot = aov(percent_col ~ tx, data = colforplot)

# from "agricolae" package
mylabels = HSD.test(anovaforplot, "tx", group = TRUE)

# anothertry = data.frame(x = c((1:6), (1:6)),
#                         y = c(4, 4, 4, 6, 6, 7.5, 5, 4, 5, 4, 4, 4),
#                         N_level = c(rep("High", 6), rep("Low", 6)),
#                         labs = c(paste(c("c", "c", "c", "b", "b", "a")), paste(c("bc", "c", "bc", "c", "c", "c"))))

# Okay, actually there are no significant pairwise differences here.
# Maybe no need for the Tukey labels.

# What's up with the Suillus compartments that were supposedly
# 100% colonized?

hundredpercent = percent_col[percent_col$percent_col == max(colforplot$percent_col),]

collabels = data.frame(N_level = c("High", "Low"),
                       x1 = c(1, 1), x2 = c(6, 6), y1 = c(95, 70), y2 = c(95, 71),
                       xstar = c(3.5, 3.5), ystar = c(103, 78),
                       lab = c("a", "b"))


colplot = ggplot(data = colforplot) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competitors, y = percent_col,
                   fill = compartment_fungus)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9,
                                             jitter.width = 0.15),
             aes(x = competitors, y = percent_col,
                  fill = compartment_fungus,
                  shape = compartment_fungus)) +
  # geom_line(aes(group = as.factor(Plant))) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  ylab("Percent fungal colonization of\nroots by compartment") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungi on roots at harvest") +
  scale_fill_manual(values = c("lightgray", "gray46", "white")) +
  scale_shape_manual(values = c(1, 16, 2)) +
  labs(shape = "Fungus", fill = "Fungus") +
  ylim(-2, 105)
  # geom_text(data = collabels, aes(x = xstar,  y = ystar, label = lab)) +
  # geom_segment(data = collabels, aes(x = x1, xend = x2,
  #                                    y = y1, yend = y2),
  #              colour = "black")

colplot_onlyclean = ggplot(data = col_onlyclean) +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(0.9),
               aes(x = competitors, y = percent_col,
                   fill = compartment_fungus)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9,
                                             jitter.width = 0.15),
             aes(x = competitors, y = percent_col,
                 fill = compartment_fungus,
                 shape = compartment_fungus)) +
  # geom_line(aes(group = as.factor(Plant))) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  ylab("Percent root mass colonized\nby compartment") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungal treatment") +
  scale_fill_manual(values = c("lightgray", "gray46", "white")) +
  scale_shape_manual(values = c(1, 16, 2)) +
  labs(shape = "Fungus", fill = "Fungus") +
  ylim(-2, 105)

colplot_onlyintended = ggplot(data = col_onlyintended) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competitors, y = percent_col,
                   fill = compartment_fungus)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9,
                                             jitter.width = 0.15),
             aes(x = competitors, y = percent_col,
                 fill = compartment_fungus,
                 shape = compartment_fungus)) +
  # geom_line(aes(group = as.factor(Plant))) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  ylab("Percent root mass colonized\nby compartment") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungal treatment") +
  scale_fill_manual(values = c("lightgray", "gray46", "white")) +
  scale_shape_manual(values = c(1, 16, 2)) +
  labs(shape = "Fungus", fill = "Fungus") +
  ylim(-2, 105)

save_plot("plots/Colonization_boxplot_by_compartment_updated.pdf",
          colplot,
          base_aspect_ratio = 1.8)

Figure = plot_grid(massplot, colplot, 
                   nrow = 2,
                   ncol = 1,
                   align = "v",
                   axis = "l",
                   rel_heights = c(2,2),
                   labels = c("a", "b"))
save_plot("plots/Mass_and_colonization_two_panel_boxplot_vertical_updated.pdf",
          Figure)

Figure = plot_grid(massplot_noletters, colplot_onlyclean, ncol = 2, align = "h",
                    labels = c("a", "b"),
                   rel_widths = c(1, 1.6))
save_plot("plots/Mass_and_colonization_two_panel_boxplot_onlyclean.pdf",
          Figure, ncol = 2)

### Trying out col plot with color indicating intended fungi ####

colplot_color = ggplot(data = colforplot) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competitors, y = percent_col,
                   fill = compartment_fungus)) +
  geom_jitter(width = 0.2,
              aes(x = competitors, y = percent_col,
                  color = bettercomp,
                  shape = compartment_fungus)) +
  # geom_line(aes(group = as.factor(Plant))) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Percent fungal colonization of\nroots by compartment") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungal treatment") +
  scale_fill_manual(values = c("lightgray", "gray46", "white")) +
  scale_shape_manual(values = c(1, 16, 2)) +
  scale_color_manual(values = c("white", "lightgoldenrod1", "darkgoldenrod1", "tan4", "plum2", "darkorchid4"),
                    name = "Fungi applied") +
  labs(shape = "Fungus", fill = "Fungus") +
  ylim(-2, 105)

colplot_color_bycompt = ggplot(data = colforplot) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competitors_reordered, 
                   y = percent_col,
                   fill = compartment_fungus)) +
  geom_point(aes(x = competitors_reordered, y = percent_col,
                  color = Fungus_attempted,
                  shape = compartment_fungus),
             position = position_jitterdodge()) +
  # geom_line(aes(group = as.factor(Plant))) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Percent fungal colonization of\nroots by compartment") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungi on roots at harvest") +
  scale_shape_manual(values = c(21, 23, 24)) +
  scale_color_manual(values = c("black", "darkgoldenrod1", "darkorchid4"),
                    name = "Fungus applied") +
  labs(shape = "Fungus at harvest", fill = "Fungus at harvest") +
  scale_fill_manual(values = c("lightgray", "gray46", "white")) +
  ylim(-2, 105)



#### COLONIZATION STATS ####
# colfortest = colforplot
colfortest = col_onlyclean
colfortest$versus = numeric(nrow(colfortest))

for (i in 1:nrow(colfortest)) {
  if (colfortest$compartment_fungus[i] == "Sp") {
    if (colfortest$competitors[i] == "Sp/None") {
      colfortest$versus[i] = "None"
    } else if (colfortest$competitors[i] == "Sp/Sp") {
      colfortest$versus[i] = "Sp"
    } else if (colfortest$competitors[i] == "Tt/Sp") {
      colfortest$versus[i] = "Tt"
    } else if (grepl("MIXED", colfortest$competitors[i])){
      colfortest$versus[i] = "Mixed"
    }
  } else if (colfortest$compartment_fungus[i] == "Tt") {
    if (colfortest$competitors[i] == "Tt/None") {
      colfortest$versus[i] = "None"
    } else if (colfortest$competitors[i] == "Tt/Tt") {
      colfortest$versus[i] = "Tt"
    } else if (colfortest$competitors[i] == "Tt/Sp") {
      colfortest$versus[i] = "Sp"
    } else if (grepl("MIXED", colfortest$competitors[i])) {
      colfortest$versus[i] = "Mixed"
    }
  } else if (colfortest$compartment_fungus[i] == "None") {
    if (colfortest$competitors[i] == "Sp/None") {
      colfortest$versus[i] = "Sp"
    } else if (colfortest$competitors[i] == "None/None") {
      colfortest$versus[i] = "None"
    } else if (colfortest$competitors[i] == "Tt/None") {
      colfortest$versus[i] = "Tt"
    }
  } else if (colfortest$compartment_fungus[i] == "MIXED") {
    if (colfortest$competitors[i] == "MIXED/Sp") {
      colfortest$versus[i] = "Sp"
    } else if (colfortest$competitors[i] == "MIXED/Tt") {
      colfortest$versus[i] = "Tt"
    }
  }
}
# 
# write_csv(colfortest, "processeddata/colonization_and_biomass_data_by_compartment_with_competition.csv")

colfortest = subset(colfortest, compartment_fungus != "None") # exclude compartments with no fungi,
# since "significant differences" in colonization based on the no spores
# vs yes spores compartments will not be informative.

colonization_test = lmer(percent_col ~ compartment_fungus * N_level * versus + (1|Plant),
                         data = colfortest)

summary(colonization_test) # nothing significant.
anovaresults = anova(colonization_test) #compartment fungus highly significant,
#marginal interaction with competitor identity
anovaresults

sink("stats_tables/colonization_by_compartment_lme_anova_updated.html")

stargazer(anovaresults, type = "html",
          # covariate.labels = c("compartment fungus",
          #                      "N level",
          #                      "competitor",
          #                      "compartment fungus:N level",
          #                      "compartment fungus:competitor",
          #                      "N level:competitor",
          #                      "compartment fungus:N level:competitor"),
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()

colonizationposthoc = emmeans(colonization_test, list(pairwise ~ compartment_fungus*N_level*versus), adjust = "tukey")

