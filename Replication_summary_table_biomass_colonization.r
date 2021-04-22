# Creating replication summary table for publication
# 2 March 2021 LB

library(cowplot)
library(tidyverse)

mydata = read_csv("processeddata/Biomass_and_colonization_forplants_with_no_Tt_contamination.csv") # from script 3-2-1
# These are data from which I've filtered:
# - any plant that had Tt on it at the end of the experiment, and shouldn't have
# - any plant with mixed colonization anywhere on its roots.
# so these are only the clean plants, the ones included in my final versions of the 
# colonization and biomass boxplots.

mysummarydata = select(mydata, Plant, N_level, competitors)

checking = mysummarydata %>% group_by(Plant) %>% summarize(count = n()) %>%
  filter(count != 1) # only one plant per ID, that's good.

thetable = mysummarydata %>% group_by(N_level, competitors) %>%
  summarize(count = n())

write_csv(thetable, "processeddata/bioandcol_replication_summary_by_compartment.csv")

