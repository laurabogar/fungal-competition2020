# FC QC: How much 15N label leaked across mesh panels?

setwd("~/Documents/Fungal competition project/fungal-competition2020/")


require(tidyverse)

minimally_processed_isotopes = read_csv("FCdata/Cleaned_processed_FC_isotope_data_July.csv")
metadata_byplant = read_csv("FCdata/percent_col_and_mass_data_by_plant.csv")

#### How much did the N-15 label leak across mesh panels? ####

unenriched = subset(minimally_processed_isotopes, enriched == 0 & 
                      (tissue == "uncolonized_roots"|tissue == "mycorrhizas"))

nmunenriched = subset(unenriched, tissue == "uncolonized_roots")

mean(unenriched$APE15N)
sd(unenriched$APE15N)

# How were unenriched plants distributed across treatments?

formerge = select(metadata_byplant, Plant, N_level)
unwithn = left_join(unenriched, formerge)
unwithn_nodup = unwithn[!duplicated(unwithn$Plant),]

unwithn_nodup %>% group_by(Actual_fungi_at_harvest, N_level) %>% summarize(total = n())
# Oh gosh, I've only got three of these. Shouldn't matter a 
# lot, though, since they're just a baseline -- as long
# as all plants have the SAME baseline for this study, it's
# low stakes.

unwithn %>% group_by(tissue, N_level) %>% summarize(total = n())


# How many unenriched plants did I have overall?
summary(as.factor(unwithn_nodup$Batch))
# Okay, this is still just three. I have more in my plant spreadsheet.

nm15N = subset(minimally_processed_isotopes, 
               Actual_fungus_by_compartment == "NM" &
                 `Receives 15N label?` == "Y")

mean(nm15N$APE15N)
sd(nm15N$APE15N)


t.test(nm15N$APE15N, unenriched$APE15N)


# This is encouraging:
# > t.test(nm15N$APE15N, unenriched$APE15N)
# 
# Welch Two Sample t-test
# 
# data:  nm15N$APE15N and unenriched$APE15N
# t = -1.1281, df = 9.3342,
# p-value = 0.2874
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.0017046817  0.0005660773
# sample estimates:
#   mean of x     mean of y 
# -5.693022e-04  1.850373e-17 