# 2-1 Leak check t test

# I'd like to perform a simple t-test to see if my mesh panels
# allowed for passive leakage of 15N without fungal transport.

# setwd("~/Documents/Fungal competition project/fungal-competition2020/")

library(tidyverse)

data = read_csv("processeddata/minimally_processed_isotope_data_July.csv")

unenriched = subset(data, enriched == 0 &
                      (tissue == "uncolonized_roots"|tissue == "mycorrhizas"))


nm15N = subset(data,
               Actual_fungus_by_compartment == "NM" &
                 received15N == "Y")

##test
library(plotly)

sp15N = subset(data,
              Actual_fungus_by_compartment == "SUIPU" &
                received15N == "Y")

forplot = subset(data,
                 (Actual_fungus_by_compartment == "NM" |
                    Actual_fungus_by_compartment == "SUIPU") &
                   received15N == "Y")

APE15N_roots_plot = ggplot(subset(data, received15N == "Y" &
                tissue == "uncolonized_roots" &
                Actual_fungus_by_compartment != "MIXED" &
                Actual_fungus_by_compartment != "OTHER")) + 
  geom_boxplot(outlier.alpha = 0,
               aes(x = Actual_fungus_by_compartment,
                 y = APE15N,
                 color = N_level)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15),
             aes(x = Actual_fungus_by_compartment,
                 y = APE15N,
                 color = N_level))

APE15N_mycorrhizas_plot = ggplot(subset(data, received15N == "Y" &
                                    tissue == "mycorrhizas" &
                                    Actual_fungus_by_compartment != "MIXED" &
                                    Actual_fungus_by_compartment != "OTHER" &
                                      Actual_fungus_by_compartment != "NM")) + 
  geom_boxplot(outlier.alpha = 0,
               aes(x = Actual_fungus_by_compartment,
                   y = APE15N,
                   color = N_level)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15),
             aes(x = Actual_fungus_by_compartment,
                 y = APE15N,
                 color = N_level))

percent_N_mycos_plot = ggplot(subset(data, tissue == "mycorrhizas" &
                                       Actual_fungus_by_compartment != "MIXED" &
                                       Actual_fungus_by_compartment != "OTHER" &
                                       Actual_fungus_by_compartment != "NM" &
                                       N_level != "NA")) + 
  geom_boxplot(outlier.alpha = 0,
               aes(x = Actual_fungus_by_compartment,
                   y = pctN,
                   color = N_level)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15),
             aes(x = Actual_fungus_by_compartment,
                 y = pctN,
                 color = N_level))

granular_data_bycompt = read_csv("processeddata/granular_mass_and_colonization_data_by_compartment.csv")
granular_nomix = granular_data_bycompt[-grep("MIXED", granular_data_bycompt$competitors),]
granular_nomix$Side = tolower(granular_nomix$Side)
granular_nomix = subset(granular_nomix, compartment_fungus != "Failed")
granular_tojoin = select(granular_nomix, Plant, Side, 
                         total_root_biomass_compartment,
                         compartment_fungus)
scalingN = right_join(data, granular_tojoin) %>% distinct()


percent_N_roots_plot = ggplot(subset(data, tissue == "uncolonized_roots" &
                Actual_fungus_by_compartment != "MIXED" &
                Actual_fungus_by_compartment != "OTHER" &
                N_level != "NA")) + 
  geom_boxplot(outlier.alpha = 0,
               aes(x = Actual_fungus_by_compartment,
                   y = pctN,
                   color = N_level)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15),
             aes(x = Actual_fungus_by_compartment,
                  y = pctN,
                  color = N_level)) +
  ylab("Percent N in uncolonized roots")

ggplotly(percent_N_roots_plot) # helpful for exploring extreme values, but kinda
# weird looking as a plot.

# What's going on with the NM outlier?

myoutlier = data[data$pctN == 3.4,] # It seems very possible to me
# that this was actually a Thelephora mycorrhiza that I failed to recognize.
# Or perhaps a root that had some other fungus -- mold? -- contamination.

C13_roots_plot = ggplot(subset(data, tissue == "uncolonized_roots" &
                                       Actual_fungus_by_compartment != "MIXED" &
                                       Actual_fungus_by_compartment != "OTHER" &
                                       N_level != "NA")) + 
  geom_boxplot(outlier.alpha = 0,
               aes(x = Actual_fungus_by_compartment,
                   y = APE13C,
                   color = N_level)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15),
             aes(x = Actual_fungus_by_compartment,
                 y = APE13C,
                 color = N_level))

pctC_roots_plot = ggplot(subset(data, tissue == "uncolonized_roots" &
                                 Actual_fungus_by_compartment != "MIXED" &
                                 Actual_fungus_by_compartment != "OTHER" &
                                 N_level != "NA")) + 
  geom_boxplot(outlier.alpha = 0,
               aes(x = Actual_fungus_by_compartment,
                   y = pctC,
                   color = N_level)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15),
             aes(x = Actual_fungus_by_compartment,
                 y = pctC,
                 color = N_level))

C13_mycorrhizas_plot = ggplot(subset(data, tissue == "mycorrhizas" &
                                          Actual_fungus_by_compartment != "MIXED" &
                                          Actual_fungus_by_compartment != "OTHER" &
                                          Actual_fungus_by_compartment != "NM" &
                                       N_level != "NA")) + 
  geom_boxplot(outlier.alpha = 0,
               aes(x = Actual_fungus_by_compartment,
                   y = APE13C,
                   color = N_level)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15),
             aes(x = Actual_fungus_by_compartment,
                 y = APE13C,
                 color = N_level))

ggplot(subset(data, tissue == "mycorrhizas" &
                Actual_fungus_by_compartment != "MIXED" &
                Actual_fungus_by_compartment != "OTHER" &
                N_level != "NA")) + 
  geom_boxplot(outlier.alpha = 0,
               aes(x = Actual_fungus_by_compartment,
                   y = pctN,
                   color = N_level)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15),
             aes(x = Actual_fungus_by_compartment,
                 y = pctN,
                 color = N_level))






##endtest


my_t_test = t.test(nm15N$APE15N, unenriched$APE15N)

# Save output to file
sink("stats_tables/leak_check_t_test.txt")

my_t_test

sink()

# We get a qualitatively similar result if we limit ourselves
# to only processing uncolonized roots (not mycorrhizas) in unenriched
# plants:
# nmunenriched = subset(unenriched, tissue == "uncolonized_roots")
# t.test(nm15N$APE15N, nmunenriched$APE15N)
