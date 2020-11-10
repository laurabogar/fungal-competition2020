# 3-2-1 Combining biomass and colonization plots into one figure

#### BIOMASS ####

# significance labels determined by anova for interaction between N level and fungal treatment.

# setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(agricolae)
library(emmeans)
library(cowplot)
library(lme4)
library(lmerTest)
library(tidyverse)
library(stargazer)

alldata = read_csv("processeddata/percent_col_and_mass_data_by_plant.csv")
compdata = read_csv("processeddata/percent_colonization_and_mass_data_by_compartment.csv")


#### BIOMASS ####
# This function maybe does automatic letters?
tx = with(alldata, interaction(N_level, Fungi))
anovaforplot = aov(total_biomass ~ tx, data = alldata)

# from "agricolae" package
mylabels = HSD.test(anovaforplot, "tx", group = TRUE)
# Oh thank goodness this matches my prior
# results and simplifies my life a lot.

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
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Total plant biomass (g)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_text(data = anothertry, aes(x, y, label = labs)) +
  xlab("Fungi on roots at harvest")

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
colforplot = compdata[!is.na(compdata$compartment_fungus),]
colforplot = subset(colforplot, compartment_fungus != "MIXED" &
                      compartment_fungus != "OTHER" &
                      competitors != "THETE" &
                      competitors != "FAILED")
colforplot = colforplot[-grep("MIXED", colforplot$competitors),]
colforplot = subset(colforplot, N_level != "None" &
                      compartment_fungus != "Failed")
colforplot$compartment_fungus = recode(colforplot$compartment_fungus,
                                        "NM" = "None",
                                        "SUIPU" = "Sp",
                                        "THETE" = "Tt")
colforplot$bettercomp = fct_relevel(colforplot$attempted,
                                   "None/None",
                                   "Sp/None",
                                   "Sp/Sp",
                                   "Tt/Sp",
                                   "Tt/None",
                                   "Tt/Tt")
colforplot$competitors_reordered = fct_relevel(colforplot$competitors,
                                    "None/None",
                                    "Sp/None",
                                    "Sp/Sp",
                                    "Tt/Sp",
                                    "Tt/None",
                                    "Tt/Tt")

colforplot$Fungus_attempted[colforplot$Plant == 6106 & colforplot$Side == "A"] = "Tt" # why not in data? unclear
colforplot$Fungus_attempted[colforplot$Plant == 6106 & colforplot$Side == "B"] = "Sp"

dead_mycos = colforplot %>% group_by(N_level, Fungus_attempted, mycofungus) %>% 
  filter(dead_tissue > 0) %>%
  summarize(n(), sum(dead_tissue))

# If including mixed:
# dead_mycos = colforplot %>% group_by(N_level, Fungus_attempted, compartment_fungus) %>% 
#   filter(dead_tissue > 0) %>%
#   summarize(n(), sum(dead_tissue))


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


collabels = data.frame(N_level = c("High", "Low"),
                       x1 = c(1, 1), x2 = c(6, 6), y1 = c(95, 70), y2 = c(95, 71),
                       xstar = c(3.5, 3.5), ystar = c(103, 78),
                       lab = c("a", "b"))


colplot = ggplot(data = colforplot) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competitors, y = percent_col,
                   fill = compartment_fungus)) +
  geom_point(aes(x = competitors, y = percent_col,
                  fill = compartment_fungus,
                  shape = compartment_fungus),
             position = position_jitterdodge()) +
  # geom_line(aes(group = as.factor(Plant))) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
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

save_plot("plots/Colonization_boxplot_by_compartment.pdf",
          colplot,
          base_aspect_ratio = 1.8)

Figure = plot_grid(massplot, colplot, 
                   nrow = 2,
                   ncol = 1,
                   align = "v",
                   axis = "l",
                   rel_heights = c(2,2),
                   labels = c("a", "b"))
save_plot("plots/Mass_and_colonization_two_panel_boxplot_vertical.pdf",
          Figure)

Figure = plot_grid(massplot, colplot, ncol = 2, align = "h",
                    labels = c("a", "b"),
                   rel_widths = c(1, 1.6))
save_plot("plots/Mass_and_colonization_two_panel_boxplot.pdf",
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
colfortest = colforplot
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
anovaresults = anova(colonization_test)
summary(anovaresults)

sink("stats_tables/colonization_by_compartment_lme_anova.html")

stargazer(anovaresults, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          summary = FALSE,
          no.space = TRUE)

sink()

colonizationposthoc = emmeans(colonization_test, list(pairwise ~ compartment_fungus*N_level*versus), adjust = "tukey")

