#3-4-FC_plant_response_boxplot

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(agricolae)
library(tidyverse)
library(cowplot)

alldata = read_csv("processeddata/percent_col_and_mass_data_by_plant.csv")


#### PLANT RESPONSE TO COLONIZATION ####

noNMplants = subset(alldata, Fungi != "None/None")

tx = with(noNMplants, interaction(N_level, Fungi))
anovanoNMplants = aov(plant_response ~ tx, data = noNMplants)

# from "agricolae" package
mylabels = HSD.test(anovanoNMplants, "tx", group = TRUE)
# Manually fed this result into the below manual annotation:
annotations = data.frame(x = c((1:5), (1:5)),
                         y = c(1, 1, 1.3, 1.3, 1.7, 1, 1, 1, 1, 1),
                         N_level = c(rep("High", 5), rep("Low", 5)),
                         labs = c(paste(c("bc", "abc", "ab", "bc", "a")), paste(c("c", "bc", "bc", "bc", "bc"))))


labels = c(High = "High N", Low = "Low N")

responseplot = ggplot(data = noNMplants) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = Fungi, y = plant_response)) +
  geom_jitter(width = 0.20,
              aes(x = Fungi, y = plant_response)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Plant response to fungal\ncolonization (log response ratio)") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_text(data = annotations, aes(x, y, label = labs)) +
  xlab("Fungal treatment")


pdf("plots/Plant_response_boxplot.pdf", width = 7, height = 5)
responseplot
dev.off()

ggsave("plots/Plant_response_boxplot.jpeg", plot = responseplot,
       device = "jpeg",
       width = 7, height = 5, units = "in")


