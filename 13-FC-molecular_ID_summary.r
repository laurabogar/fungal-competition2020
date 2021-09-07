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

nodupes = alltogether %>% unique()

mysummary = alltogether %>%
  group_by(compartment_fungus, `Matched prediction?`, `At level`) %>%
  summarize(count = n())

anotherone = nodupes %>%
  group_by(compartment_fungus, Morphospecies, `Actual identity`) %>%
  summarize(count = n())
