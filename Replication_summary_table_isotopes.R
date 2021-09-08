# Creating replication summary table for publication
# 2 March 2021 LB

library(cowplot)
library(tidyverse)

mydata = read_csv("processeddata/isotope_and_plant_metadata_with_competition_coded_clearly_INCLUDING_MIXED_and_pctCN.csv") # from script 4-alternative
leafinfo = read_csv("processeddata/Cleaned_processed_FC_leaf_isotopes.csv")

mydata = mydata %>% mutate(competitors = Fungi)

mysummarydata = select(mydata, Plant, Side, N_level, competitors, compartment_fungus, received15N)

mysum_noNAs = mysummarydata %>% filter(is.na(competitors) == FALSE)

mysum_filtered = mysum_noNAs %>% distinct()

# mysummarydata_nomixed = mysum_filtered[!grepl("MIXED", mysum_filtered$competitors),]

mysummarydata_nomixed = mysum_filtered[!grepl("Mixed", mysum_filtered$compartment_fungus),]

checking = mysummarydata_nomixed %>% group_by(Plant) %>% summarize(count = n()) %>%
  filter(count != 2) # This used to show just two plants; now there are 3, because of the mixed compt in 6043b.

test = mydata[mydata$Plant %in% checking$Plant,]
# 6043 was revealed to be a mixed compartment with DNA data.
# But why do I only have one entry each for plants 6013 and 6072?

# Based on my notes: Both plants were from the low N treatment.
# 6013 was Sp/NM, and is missing the NM side. 
# This sample did not give usable data at the UCSC facility.
# 6072 (Tt/Tt) is missing one of its Tt sides, side b. Those samples
# also did not give usable data at the facility.

# I have also just double checked -- these plants did NOT make it into the
# "isotope_and_plant_metadata_FOR_N_ANALYSES_and_exchange_rates_withpctN.csv"
# file that I've been using for nitrogen data analyses. This makes sense,
# since neither of these root compartments received 15N label.


finaldata = mysummarydata_nomixed %>% filter(!Plant %in% checking$Plant)

# Looking into weird asymmetries in the table:

# # I have a low-N Tt/None plant showing up as enriched on the none side,
# # but nonexistent on the Tt side.
# myTtnone = finaldata %>% filter(competitors == "Tt/None", N_level == "Low")
# # Examining the data, this problem plant is number 6061 (side a).
# # What happened to side b?
# fullTtnone = mysummarydata %>% filter(competitors == "Tt/None", N_level == "Low")
# # Well, it's absent even in this data pre-filtering.
# # FOUND IT!! And fixed the underlying scripts. Shouldn't be a problem now.

# # I also have a low N Tt/Tt plant showing up that either shouldn't be there, or
# # is missing its buddy....
# 
# myTtTt = finaldata %>% filter(competitors == "Tt/Tt", N_level == "Low")
# 
# # But it seems to have disappeared? I think the fixes I put in place while
# # tracking down the other missing plant took care of the issue -- the
# # lopsided entry was a Tt/None compartment miscoded as Tt/Tt.

thetable = finaldata %>% group_by(competitors, compartment_fungus, N_level, received15N) %>%
  summarize(count = n())

write_csv(thetable, "processeddata/replication_summary_by_compartment.csv")


### Summarizing plants used for leaf isotopes ####

plants_for_leaf_isotopes = finaldata %>% 
  group_by(Plant) %>% 
  select(Plant, N_level, competitors) %>%
  unique()

write_csv(plants_for_leaf_isotopes, "processeddata/plants_for_leaf_isotope_measurements.csv")

# Other notes:

# I am going to have to manually remove plant 6076 (Tt/None, low N) from 
# this table, because ultimately this plant wasn't usable for mycorrhizal calculations. It should have been Sp/None,
# but was contaminated with Tt. The only thing it might have been good for was isotope work;
# however, I also couldn't find any Tt mycos when I harvested it, despite seeing Tt crust on the surface
# of the soil. (I probably missed them, since young ones look a lot like roots
# -- otherwise Tt is a better saprobe than we'd realized?)

# In Tt/Tt low N, I will NOT remove plant 6030 from the count, but we should be aware
# that it is a plant for which we have isotopes for only half the root system.
# Based on my notes, it had "lame mycos" (probably old and not numerous) at harvest.
# And I see no evidence of our ever having tin-balled that sample, or sent for analysis.
# It's unclear if the envelope was lost in the shuffle, or what happened, 
# but we don't have data for mycos from half the root system, so it's not useful
# specifically for paired analyses of carbon flow. This plant still might be useful
# for figuring out how much C Tt mycos got, though. (The non-15N-labeled side is 
# the one for which we have data.)

# Plant 6072 (Tt/Tt low N) is also available as only half a root system,
# in this case because the sample yielded unusable data at UCSC. Still potentially
# useful for other non-pairwise isotope calculations, though, so I should keep it around.
# Adding it back into the table (was excluded due to missing data earlier)
