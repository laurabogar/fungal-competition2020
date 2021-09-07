# Initial cleaning of isotope data for resource trading/"fungal competition" project

# setwd("~/Documents/2018-2019/Fungal competition/")
# setwd("~/Documents/Fungal competition project/fungal-competition2020/")

require(tidyverse)

R13C = 0.0112372 # Slater et al. 2011
R15N = 0.0036765 # Slater et al. 2011

# fcdata = read_csv("FC_isotope_data_from_UCSC.csv") # to revert to old analysis, uncomment
leafdata = read_csv("rawdata/leaf_isotopes_UCSB.csv")
metadata = read_csv("rawdata/Fungal_competition_plant_tracking.csv")

leafdata = mutate(leafdata, Plant = plant)

metadata_oneperplant = metadata %>%
  select(Plant, Fungal_treatment, N_level, `Harvest date`, Batch, Harvest_time, Actual_fungi_at_harvest, Failed_split)

# metadata[is.na(metadata$Batch),]

isowithmeta = left_join(leafdata, metadata_oneperplant)
isowithmeta$enriched = numeric(nrow(isowithmeta))
isowithmeta$enriched[grep("FC", isowithmeta$Batch)] = 1

isowithmeta = isowithmeta[!duplicated(isowithmeta$Plant),]

### Calculating other isotope measures ####
atompct = function(delta, R) {
  output = 100/(1+(1/(((delta/1000)+1)*R)))
  return(output)
}

isowithmeta$atmpct13C = numeric(nrow(isowithmeta))
isowithmeta$atmpct15N = numeric(nrow(isowithmeta))

for (i in 1:nrow(isowithmeta)) {
  isowithmeta$atmpct13C[i] = atompct(isowithmeta$d13C[i], R13C)
  isowithmeta$atmpct15N[i] = atompct(isowithmeta$d15N[i], R15N)
}

baseline = isowithmeta[isowithmeta$enriched == 0,]


meanbaseline = baseline %>% summarize(meand13C = mean(d13C), 
                                      meand15N = mean(d15N),
                                      meanatmpct13C = mean(atmpct13C),
                                      meanatmpct15N = mean(atmpct15N))

# isowithmeta$mmol13Cexcess = numeric(nrow(isowithmeta))
# isowithmeta$mmol15Nexcess = numeric(nrow(isowithmeta))
# 
# for (i in 1:nrow(isowithmeta)) {
#   deltadiff = (isowithmeta$d13C[i] - meanbaseline$meand13C)
#   isowithmeta$mmol13Cexcess[i] = R13C*deltadiff
#   ndeltadiff = (isowithmeta$d15N[i] - meanbaseline$meand15N)
#   isowithmeta$mmol15Nexcess[i] = R15N*ndeltadiff
# }

# Well, at least I've got a couple of representatives of each host
# in this data set.

isowithmeta$mmol13Cexcess = numeric(nrow(isowithmeta))
isowithmeta$mmol15Nexcess = numeric(nrow(isowithmeta))

isowithmeta$APE13C = numeric(nrow(isowithmeta))
isowithmeta$APE15N = numeric(nrow(isowithmeta))

for (i in 1:nrow(isowithmeta)) {
  deltadiff = (isowithmeta$d13C[i] - meanbaseline$meand13C)
  isowithmeta$mmol13Cexcess[i] = R13C*deltadiff
  
  ndeltadiff = (isowithmeta$d15N[i] - meanbaseline$meand15N)
  isowithmeta$mmol15Nexcess[i] = R15N*ndeltadiff
  
  isowithmeta$APE13C[i] = isowithmeta$atmpct13C[i] - meanbaseline$meanatmpct13C
  isowithmeta$APE15N[i] = isowithmeta$atmpct15N[i] - meanbaseline$meanatmpct15N
  
}

isowithmeta_leafnames = select(isowithmeta, everything(), -plant) %>%
  rename(leaf.d13C = d13C,
         leaf.d15N = d15N,
         leaf.pctC = pctC,
         leaf.pctN = pctN,
         leaf.atmpct13C = atmpct13C,
         leaf.atmpct15N = atmpct15N,
         leaf.mmol13Cexcess = mmol13Cexcess,
         leaf.mmol15Nexcess = mmol15Nexcess,
         leaf.APE13C = APE13C,
         leaf.APE15N = APE15N)

write_csv(isowithmeta_leafnames, "processeddata/Cleaned_processed_FC_leaf_isotopes.csv")
