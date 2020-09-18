# 3-2-1 Combining biomass and colonization plots into one figure

test = alldata[alldata$Fungus == "Tt/Sp",]

#### BIOMASS ####

# significance labels determined by anova for interaction between N level and fungal treatment.

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(agricolae)
library(cowplot)
library(tidyverse)

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
  xlab("Fungal treatment")

pdf("plots/Biomass_boxplot.pdf", width = 9, height = 5)
massplot
dev.off()

ggsave("plots/Biomass_boxplot.jpeg", plot = massplot,
       device = "jpeg",
       width = 9, height = 5, units = "in")

#### COLONIZATION ####

### Analyses ###
test = colforplot[is.na(colforplot$compartment_fungus),]

# colforplot = subset(alldata, Fungi != "None/None")
# remove the "NA" values because these are all failed splits
# except 6105, which is a plant for which I don't have colonization
# data -- I think I harvested it before I had developed that protocol? 
# Either that or I lost the relevant tissue envelopes.
colforplot = compdata[!is.na(compdata$compartment_fungus),]
colforplot = subset(colforplot, compartment_fungus != "MIXED" &
                      compartment_fungus != "OTHER" &
                      competitors != "THETE")
colforplot = colforplot[-grep("MIXED", colforplot$competitors),]
colforplot = subset(colforplot, N_level != "None")

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
                       x1 = c(1, 1), x2 = c(5, 5), y1 = c(80, 50), y2 = c(81, 51),
                       xstar = c(3, 3), ystar = c(88, 58),
                       lab = c("a", "b"))


colplot = ggplot(data = colforplot) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competitors, y = percent_col,
                   fill = compartment_fungus)) +
  geom_jitter(width = 0.20,
              aes(x = competitors, y = percent_col,
                  fill = compartment_fungus,
                  shape = compartment_fungus)) +
  # geom_line(aes(group = as.factor(Plant))) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Percent fungal colonization of\nroot system (by mass)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungal treatment") +
  scale_fill_grey() +
  scale_shape_manual(values = c(16, 1, 2))
  # geom_text(data = collabels, aes(x = xstar,  y = ystar, label = lab)) +
  # geom_segment(data = collabels, aes(x = x1, xend = x2,
  #                                    y = y2, yend = y2),
  #              colour = "black")

gray.colors(3)
pdf("plots/Colonization_boxplot.pdf", width = 7, height = 5)
colplot
dev.off()

Figure = plot_grid(massplot, colplot, ncol = 2, align = "h",
                    labels = c("A", "B"))
save_plot("plots/MAIN_Mass_and_colonization_two_panel_boxplot.pdf",
          Figure, ncol = 2)