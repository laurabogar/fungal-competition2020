# Initial cleaning of isotope data for resource trading/"fungal competition" project

# setwd("~/Documents/2018-2019/Fungal competition/")
setwd("~/Documents/Fungal competition project/fungal-competition2020/")

require(tidyverse)

R13C = 0.0112372 # Slater et al. 2011
R15N = 0.0036765 # Slater et al. 2011

# fcdata = read_csv("FC_isotope_data_from_UCSC.csv") # to revert to old analysis, uncomment
leafdata = read_csv("rawdata/leaf_isotopes_UCSB.csv")
metadata = read_csv("rawdata/Fungal_competition_plant_tracking.csv")

fcdata = select(fcdata, Plant = plant, Side = side, everything())
metadata$Side = tolower(metadata$Side)
fcdata$Side = tolower(fcdata$Side)

# metadata[is.na(metadata$Batch),]

isowithmeta = left_join(fcdata, metadata)
isowithmeta$enriched = numeric(nrow(isowithmeta))
isowithmeta$enriched[grep("FC", isowithmeta$Batch)] = 1

summary(as.factor(isowithmeta$tissue))
isowithmeta$tissue = tolower(isowithmeta$tissue)
summary(as.factor(isowithmeta$tissue))
isowithmeta$tissue[isowithmeta$tissue == "hyphae/rhizomorphs"] = "hyphae"
isowithmeta$tissue[grep("mycos", isowithmeta$tissue)] = "mycorrhizas"
isowithmeta$tissue[grep("nm|tfr", isowithmeta$tissue)] = "uncolonized_roots"
isowithmeta$mycofungus = isowithmeta$Actual_fungus_by_compartment

mixed_to_check = subset(isowithmeta, Actual_fungus_by_compartment == "MIXED")
mixed_to_check = select(mixed_to_check, Plant, Side, "Identifier 2")
tocheck = select(mixed_to_check, Plant, Side)
tocheck = arrange(tocheck, Plant)

labels = paste(tocheck$Plant, tocheck$Side)

mixed_ones_to_check = unique(labels)

mixed_ones_to_check

closerlook = arrange(mixed_to_check, Plant, Side)
closerlook = select(closerlook, Plant, Side, Harvest_notes)
write_csv(closerlook, "harvest_notes.csv")

# [1] "6002 b" "6028 a" "6034 b" "6052 b"
# [5] "6056 a" "6057 a" "6067 b" "6098 a"

# We seem to have no isotope data for 6034b (uncolonized) or 6098a (mycos). I think these
# are the samples that failed at the facility.

#6002b: Identifier 2 = 3_F10, it's SUIPU, based on mass and notes
#6002b: Identifier 2 = 4_C3, it's THETE, based on mass and notes
# 6028a: Identifier 2 = 2_B12, it's SUIPU
# 6028a: Identifier 2 = 4_H9, it's THETE
# 6034b: Identifier 2 = 3_B4, it's THETE
# 6034b: Identifier 2 = 4_D8, it's SUIPU 
# 6052b: Identifier 2 = 1_D4, it's SUIPU
# 6052b: Identifier 2 = 2_C12, it's THETE
# 6056a: Identifier 2 = 2_D2, it's SUIPU
# 6056a: Identifier 2 = ... looks like I don't have THETE tips for some reason.
# 6057a: Identifier 2 = 3_B12, it's SUIPU
# 6057a: Identifier 2 = 2_C6, it's THETE
# 6067b: Identifier 2 = 4_E9, it's SUIPU
# 6067b: Identifier 2 = 4_H4, it's THETE
# 6098a: Identifier 2 = 1_C9, it's SUIPU... but no myco isotope data

mixed_compartments_to_update = c("6002b", "6028a", "6052b", "6057a", "6067b")

#6002b: Identifier 2 = 3_F10, it's SUIPU, based on mass and notes
#6002b: Identifier 2 = 4_C3, it's THETE, based on mass and notes
#6057a: ID 2 = 3_B12, it's SUIPU
#6057a: ID2 = 2_C6, it's THETE
# 6067b: ID2 = 4_E9, it's SUIPU
# 6067b: ID2 = 4_H4, it's THETE
# 6098a: Apparently only ran one mycorrhiza sample, and it failed. Boo.

SUIPU_ids = c("3_F10", "2_B12", "4_D8", "1_D4", "2_D2", "3_B12", "4_E9")
THETE_ids = c("4_C3", "4_H9", "3_B4", "2_C12", "2_C6", "4_H4")


# Let's straighten out those mixed compartments

isowithmeta$mycofungus[isowithmeta$`Identifier 2` %in% THETE_ids] = "THETE"

isowithmeta$mycofungus[isowithmeta$`Identifier 2` %in% SUIPU_ids] = "SUIPU"

isowithmeta$uniqueIds = paste(isowithmeta$Plant, isowithmeta$Side, isowithmeta$tissue, isowithmeta$mycofungus)

isowithmeta$replicated = duplicated(isowithmeta$uniqueIds)

isowithmeta[isowithmeta$replicated == TRUE,] # 15 of these.
# Some were TFRs (terminal fine roots) that I recoded as 
# uncolonized roots; others appear to just have been samples
# that I genuinely replicated for some reason.

isowithmeta = subset(isowithmeta, is.na(isowithmeta$d13C) == FALSE)

# Let's average out the replicates.

for (i in 1:nrow(isowithmeta)) {
  for (j in 1:nrow(isowithmeta)) {
    if (isowithmeta$uniqueIds[i] == isowithmeta$uniqueIds[j]) {
      isowithmeta$d13C[i] = mean(isowithmeta$d13C[i], isowithmeta$d13C[j])
      isowithmeta$pctC[i] = mean(isowithmeta$pctC[i], isowithmeta$pctC[j])
      isowithmeta$mcgC[i] = mean(isowithmeta$mcgC[i], isowithmeta$mcgC[j])
      isowithmeta$d15N[i] = mean(isowithmeta$d15N[i], isowithmeta$d15N[j])
      isowithmeta$pctN[i] = mean(isowithmeta$pctN[i], isowithmeta$pctN[j])
      isowithmeta$mcgN[i] = mean(isowithmeta$mcgN[i], isowithmeta$mcgN[j])
    }
  }
}

isowithmeta = subset(isowithmeta, replicated == FALSE)

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

baseline = isowithmeta[isowithmeta$enriched == 0 & isowithmeta$tissue != "hyphae",]


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

# write_csv(isowithmeta, "Cleaned_processed_FC_isotope_data.csv")
write_csv(isowithmeta, "processeddata/Cleaned_processed_FC_isotope_data_July.csv")
