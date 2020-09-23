
setwd("~/Documents/Fungal competition project/fungal-competition2020/")

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
# How did I lose 18 plants?

summary(metadata$Failed_split == "Y") # 8/2 = 4 failed splits
hmm = metadata[metadata$Failed_split == "Y",]


sum(grepl("MIXED|OTHER", metadata$Actual_fungi_at_harvest)) #16/2 = 8 microcosms where we had a mixed result
more = metadata[grepl("MIXED|OTHER|FAILED", metadata$Actual_fungi_at_harvest),]

# I'm still missing six plants that need explaining

letslook = metadata[!metadata$Plant %in% plantleveldata$Plant,]
letslook = letslook[!letslook$Plant %in% more$Plant,]
letslook = letslook[!letslook$Plant %in% hmm$Plant,]

# > length(unique(more$Plant))
# [1] 9
# > length(unique(hmm$Plant))
# [1] 6
# > length(unique(letslook$Plant))
# [1] 4 # but one of these doesn't count because it was the "no N" control plant

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