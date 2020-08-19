# 2-1 Leak check t test

# I'd like to perform a simple t-test to see if my mesh panels
# allowed for passive leakage of 15N without fungal transport.

setwd("~/Documents/Fungal competition project/fungal-competition2020/")

data = read_csv("processeddata/minimally_processed_isotope_data_July.csv")

unenriched = subset(data, enriched == 0 &
                      (tissue == "uncolonized_roots"|tissue == "mycorrhizas"))


nm15N = subset(data,
               Actual_fungus_by_compartment == "NM" &
                 received15N == "Y")


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
