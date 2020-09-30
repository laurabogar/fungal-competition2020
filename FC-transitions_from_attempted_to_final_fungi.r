
# setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(tidyverse)
library(stargazer)

plantleveldata = read_csv("processeddata/plant_level_biomass_colonization_data.csv")
compartmentleveldata = read_csv("processeddata/percent_colonization_and_mass_data_by_compartment.csv")
compartmentleveldata = subset(compartmentleveldata, 
                              is.na(compartment_fungus) == F &
                                is.na(Fungus_attempted) == F &
                                is.na(N_level) == F)

plantleveldata$Fungal_treatment = recode(plantleveldata$Fungal_treatment,
                                         "THETE/THETE" = "Tt/Tt",
                                         "THETE/NM" = "Tt/None",
                                         "THETE/SUIPU" = "Tt/Sp",
                                         "SUIPU/SUIPU" = "Sp/Sp",
                                         "SUIPU/NM" = "Sp/None",
                                         "NM/NM" = "None/None")

#### TRANSITIONS BY COMPARTMENT ####
compt_transitions = select(compartmentleveldata, plant = Plant,
                           N_level,
                           # competitors_attempted = attempted,
                           # competitors_realized = competitors,
                           compartment_attempted = Fungus_attempted,
                           realized = compartment_fungus)

compt_summary = compt_transitions %>% group_by(N_level, 
                                               # competitors_attempted, 
                                               # competitors_realized, 
                                               compartment_attempted, 
                                               realized) %>% summarize(count = n())

sink("stats_tables/fungal_transition_summary_table.html")

stargazer(compt_summary, type = "html",
          summary = FALSE,
          no.space = TRUE)

sink()

# It looks to me like Sp transition to non-Sp compartments
# way more often in high N than low. Let's check with a Fisher's exact test.

# Fungus at harvest:	
#   N level	Sp	Not Sp
# High	11	19
# Low	25	6

mytable = rbind(c(11, 19), c(25, 6))
mytable = matrix(dimnames = list(N_level = c("High", "Low"),
                                 Fungus_at_harvest = c("Sp", "Not Sp")))
fisher.test(mytable)
# p = 0.0006726

### TRANSITIONS TO WRONG FUNGI BY PLANT ####
# plant level data appears to exclude mixed compartments. I don't want this.
transitions_tbl = select(compartmentleveldata, 
                         plant = Plant,
                         N_level,
                         competitors_attempted = attempted,
                         competitors_realized = competitors)

transitions_tbl = transitions_tbl[!duplicated(transitions_tbl$plant),] # get one row per plant


transsummary = transitions_tbl %>% group_by(N_level, competitors_attempted, competitors_realized) %>% summarize(count = n())
transsummary$realized = as.character(transsummary$realized)

sink("stats_tables/fungal_transition_summary_table_byplant.html")

stargazer(transsummary, type = "html",
          summary = FALSE,
          no.space = TRUE)

sink()

#### Fisher's exact test with independent samples ####
# Is my significant result robust to
# using only truly independent samples?
compt_transitions_oneperplant = subset(compt_transitions, compartment_attempted == "Sp")
compt_transitions_oneperplant = compt_transitions_oneperplant[!duplicated(compt_transitions_oneperplant$plant),] # arbitrarily take one compartment per plant
compt_summary_oneperplant = compt_transitions_oneperplant %>% group_by(N_level, 
                                               # competitors_attempted, 
                                               # competitors_realized, 
                                               compartment_attempted, 
                                               realized) %>% summarize(count = n())



# Fungus at harvest:	
# N level	  Sp	  Not Sp
# High	    9	    12
# Low	      18	    5

table2 = rbind(c(9,12), c(18,5))
table2 = matrix(dimnames = list(N_level = c("High", "Low"),
                                 Fungus_at_harvest = c("Sp", "Not Sp")))
fisher.test(table2)
# p = 0.02907

#### Making a plot ####

compartmentleveldata$percent_Tt = numeric(nrow(compartmentleveldata))
compartmentleveldata$percent_Sp = numeric(nrow(compartmentleveldata))

for (i in 1:nrow(compartmentleveldata)) {
  if (compartmentleveldata$compartment_fungus[i] == "Tt") {
    compartmentleveldata$percent_Tt[i] = compartmentleveldata$percent_col[i]
  }
}

transition_plot = ggplot(data = compartmentleveldata) +
  geom_jitter(aes(x =))

colplot = ggplot(data = colforplot) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competitors, y = percent_col,
                   fill = compartment_fungus)) +
  geom_jitter(width = 0.2,
              aes(x = competitors, y = percent_col,
                  fill = compartment_fungus,
                  shape = compartment_fungus)) +
  # geom_line(aes(group = as.factor(Plant))) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Percent fungal colonization of\nroots by compartment") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungal treatment") +
  scale_fill_manual(values = c("lightgray", "gray46", "white")) +
  scale_shape_manual(values = c(1, 16, 2)) +
  labs(shape = "Fungus", fill = "Fungus") +
  ylim(-2, 105)