# 3-2 FC biomass boxplot

# significance labels determined by anova for interaction between N level and fungal treatment.

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(agricolae)
library(cowplot)
library(tidyverse)

alldata = read_csv("processeddata/percent_col_and_mass_data_by_plant.csv")


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


