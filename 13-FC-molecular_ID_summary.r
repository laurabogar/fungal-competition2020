# Summarizing molecular ID results

library(tidyverse)

mol_IDs = read_csv("rawdata/Fungal competition -- molecular IDs summary.csv")
compartmentleveldata = read_csv("processeddata/percent_colonization_and_mass_data_by_compartment.csv")

molID_plantside = separate(mol_IDs, `Envelope ID`, c("Plant", "Side"), sep = 4) %>%
  mutate(Plant = as.numeric(Plant))

compartment_contents = select(compartmentleveldata, 
                              Plant,
                              Side,
                              N_level,
                              Batch,
                              Fungus_attempted,
                              compartment_fungus,
                              competitors,
                              attempted,
                              enriched)

alltogether = left_join(molID_plantside, compartment_contents)

alltogether$tissue = alltogether$Contents

alltogether$tissue[grepl("myco", alltogether$Contents)] = "mycorrhiz"
alltogether$tissue[grepl("hypha", alltogether$Contents)] = "hyphae"

sum(grepl("hypha", alltogether$Contents))

nodupes = alltogether %>% unique()
nodupes$unique_number = c(1:nrow(nodupes))

toexamine = nodupes %>% select(Plant, Side, tissue, Morphospecies, `Actual identity`)
Tt = toexamine[grepl("Thelephora", toexamine$`Actual identity`),]
Sp = toexamine[grepl("Suillus", toexamine$`Actual identity`),]
justTtSp = rbind(Tt, Sp)

finalset = nodupes %>% filter(str_detect(`Actual identity`, "Thelephora|Suillus|Hebeloma"))
finalset$plantside = paste(finalset$Plant, finalset$Side)

scrutiny = finalset$plantside[duplicated(finalset$plantside)]
sameplantside = finalset[finalset$plantside %in% scrutiny,]
# I'm checking here for any root compartments where I ended up
# with multiple sequence IDs for a single compartment.
# (Tried amplifying both hyphae and mycos wherever available, and
# in most cases only one of the two approaches worked, so worth keeping.)
# Duplicate entries: 6027A (Tt mycos and hyphae, mycos matched to genus and hyphae to species)
# 6052B: Collected at least 4 envelopes of mycos (Sp %, Tt %, Sp extra, Tt extra)
# Extracted Sp envelopes 3 times. One extraction came up as Tt the first time, Sp the second time it was sequenced (same extraction, diff dilutions).
# Another Sp extraction came out as Tt, and a third as Sp. The Tt extraction came out as Hebeloma.
# I think it's worth keeping one instance of the Sp envelope in the table
# (mark as "wrong ID" but acknowledge in the table that a 1:10 dilution gave Sp),
# and keep the record of Tt coming up as Hebeloma to be upfront about potential contamination.
# 6094B: Hyphae and mycos both came out as Tt.
toremove= c(55, 18, 36, 37)

nodupes_filtered = nodupes %>%
  filter(unique_number %in% toremove == FALSE)

mysummary = nodupes_filtered %>%
  group_by(compartment_fungus, `Matched prediction?`, `At level`) %>%
  summarize(count = n())

MAINSUMMARY = nodupes_filtered %>%
  group_by(compartment_fungus, Morphospecies, `Actual identity`) %>%
  summarize(count = n())

bytissue = nodupes_filtered %>%
  group_by(compartment_fungus, tissue, Morphospecies, `Actual identity`) %>%
  summarize(count = n())
