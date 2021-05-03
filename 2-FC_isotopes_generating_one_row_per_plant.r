# Continued isotopic data wrangling, fungal competition project. 
# Producing one-row-per-plant and one-row-per-compartment dataframes for analyses.

# setwd("~/Documents/2018-2019/Fungal competition/")
setwd("~/Documents/Fungal competition project/fungal-competition2020/")

require(tidyverse)
require(cowplot)
library(naniar)

# This function was important for data wrangling 10/28/2020:
# Got it here: https://stackoverflow.com/questions/18142117/how-to-replace-nan-value-with-zero-in-a-huge-data-frame/18143097

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

# isotopes = read_csv("Cleaned_processed_FC_isotope_data.csv")
isotopes = read_csv("processeddata/Cleaned_processed_FC_isotope_data_July.csv")
percent_col = read_csv("rawdata/percent_colonization_and_mass_data_by_compartment.csv")

# evencleaner = isotopes[!isotopes$Failed_split == "Y",] # This is a mistake? because I need one of these as a baseline plant.
evencleaner = subset(isotopes, is.na(Plant) == FALSE)

sum(evencleaner$Roots_dead_at_harvest == "Y")
evencleaner[evencleaner$Limited_colonization == "Y",]

evencleaner$Side[evencleaner$Plant == 6004] = "b" # found in notes -- failed split, just one side.
evencleaner$mycofungus[evencleaner$Plant == 6004] = "SUIPU"
evencleaner$Actual_fungus_by_compartment[evencleaner$Plant == 6004] = "SUIPU"
evencleaner$Actual_fungi_at_harvest[evencleaner$Plant == 6004] = "SUIPU"


evencleaner$mycofungus[evencleaner$Plant == 6059 & evencleaner$Side == "b"] = "THETE" # found in notes; roots were dead but collected THETE rhizomorphs
evencleaner$Actual_fungus_by_compartment[evencleaner$Plant == 6059 & evencleaner$Side == "b"] = "THETE" # found in notes; roots were dead but collected THETE rhizomorphs
evencleaner$Actual_fungi_at_harvest[evencleaner$Plant == 6059 & evencleaner$Side == "b"] = "THETE" # found in notes; roots were dead but collected THETE rhizomorphs

evencleaner$Actual_fungi_at_harvest[evencleaner$Plant == 6061 & evencleaner$Side == "b"] = "THETE/NM" # coded wrong, but notes clearly define side a as NM and b as THETE.


# evencleaner$mycofungus[evencleaner$tissue != "mycorrhizas"] = NA
# I'm not super sure what I meant by this... Keep it in for
# now, I guess.


mydata = select(evencleaner, 
                Plant, Side, tissue, Amount,
                Fungal_treatment:Fungus_attempted,
                received15N = `Receives 15N label?`,
                Batch, Actual_fungi_at_harvest,
                Actual_fungus_by_compartment, mycofungus,
                enriched, pctC, pctN, CNratio, Failed_split, atmpct13C:APE15N)



# This isn't making a lot of sense because
# I have a bunch of negative APE15N values,
# which screw up the signs of these ratios.

min(mydata$APE15N)
# Let's just add 0.0025 to each APE15N value to get them
# all to positive/more meaningful numbers...
# Actually, the above was only valuable when I was trying
# to visualize APE with log-transformation. I think,
# in reality, this is less helpful than it is misleading.

# Question 2: Does the ratio of 13C in mycorrhizas to 
# 15N in coarse roots on the 15N-receiving (focal)
# sides of the microcosms track with N level or 
# competition treatment?

# 1) Calculate this per-plant metric

onerowperplant = select(mydata, Plant, Side, tissue, enriched,
                        received15N, N_level, Failed_split, Actual_fungi_at_harvest,
                        Actual_fungus_by_compartment, mycofungus, APE13C,
                        APE15N, pctC, pctN, CNratio)

onerowperplant$mycofungus[onerowperplant$mycofungus == "MIXED"] = "NM"


# write_csv(onerowperplant, "minimally_processed_isotope_data.csv")
write_csv(onerowperplant, "processeddata/minimally_processed_isotope_data_July.csv")

carbon = select(onerowperplant, Plant, Side, tissue, mycofungus, APE13C, pctC)
carbonspread_13C = carbon %>% spread(tissue, APE13C)
carbonspread_pctC = carbon %>% spread(tissue, pctC)
carbonfor13C = select(carbonspread_13C, Plant, Side, 
                      mycofungus, 
                      hyphae.APE13C = hyphae,
                      mycorrhizas.APE13C = mycorrhizas,
                      uncolonized_roots.APE13C = uncolonized_roots)

carbontomerge = select(carbonspread_pctC,
                       Plant, Side, mycofungus,
                       hyphae.pctC = hyphae,
                       mycorrhizas.pctC = mycorrhizas,
                       uncolonized_roots.pctC = uncolonized_roots)

carbonspread = left_join((carbonfor13C %>% group_by(Plant, Side, mycofungus) %>%
                    mutate(id = row_number())), carbontomerge %>% 
                    group_by(Plant, Side, mycofungus) %>%
                    mutate(id = row_number()))


# carbonspread$compid = paste(carbonspread$Plant, carbonspread$Side)
# for (i in 1:nrow(carbonspread)) {
#   if (carbonspread$mycofungus[i] == "MIXED") {
#     carbonspread$uncolonized_roots[carbonspread$compid %in% carbonspread$compid[i]] = carbonspread$uncolonized_roots[i]
#   }
# }

carbonincludingmixed = carbonspread

carbon_onerow = carbonspread %>% group_by(Plant, Side, mycofungus) %>% 
  summarize(hyphae.APE13C = mean(hyphae.APE13C, na.rm = TRUE),
            mycorrhizas.APE13C = mean(mycorrhizas.APE13C, na.rm = TRUE),
            uncolonized_roots.APE13C = mean(uncolonized_roots.APE13C, na.rm = TRUE),
            hyphae.pctC = mean(hyphae.pctC, na.rm = TRUE),
            mycorrhizas.pctC = mean(mycorrhizas.pctC, na.rm = TRUE),
            uncolonized_roots.pctC = mean(uncolonized_roots.pctC, na.rm = TRUE))

carbon_onerow[is.nan.data.frame(carbon_onerow)] <- NA

# carbonspread = subset(carbonspread, mycofungus != "MIXED")

# carbonfinal = select(carbonspread, Plant, Side,
#                      hyphae.APE13C = hyphae,
#                      mycorrhizas.APE13C = mycorrhizas,
#                      uncolonized_roots.APE13C = uncolonized_roots,
#                      mycofungus)
# 
# carbonfinal_includingmixed = select(carbonincludingmixed, Plant, Side,
#                      hyphae.APE13C = hyphae,
#                      mycorrhizas.APE13C = mycorrhizas,
#                      uncolonized_roots.APE13C = uncolonized_roots,
#                      mycofungus)

nitrogen = select(onerowperplant, Plant, Side, tissue, mycofungus, APE15N, pctN)

nitrospread_15N = nitrogen %>% spread(tissue, APE15N)
nitrospread_pctN = nitrogen %>% spread(tissue, pctN)

nitro_15N_tomerge = select(nitrospread_15N, Plant, Side, 
                      mycofungus, 
                      hyphae.APE15N = hyphae,
                      mycorrhizas.APE15N = mycorrhizas,
                      uncolonized_roots.APE15N = uncolonized_roots)

nitro_pctN_tomerge = select(nitrospread_pctN,
                       Plant, Side, mycofungus,
                       hyphae.pctN = hyphae,
                       mycorrhizas.pctN = mycorrhizas,
                       uncolonized_roots.pctN = uncolonized_roots)

nitrospread = left_join((nitro_15N_tomerge %>% group_by(Plant, Side, mycofungus) %>%
                            mutate(id = row_number())), nitro_pctN_tomerge %>% 
                           group_by(Plant, Side, mycofungus) %>%
                           mutate(id = row_number()))

nitrogen_onerow = nitrospread %>% group_by(Plant, Side, mycofungus) %>% 
  summarize(hyphae.APE15N = mean(hyphae.APE15N, na.rm = TRUE),
            mycorrhizas.APE15N = mean(mycorrhizas.APE15N, na.rm = TRUE),
            uncolonized_roots.APE15N = mean(uncolonized_roots.APE15N, na.rm = TRUE),
            hyphae.pctN = mean(hyphae.pctN, na.rm = TRUE),
            mycorrhizas.pctN = mean(mycorrhizas.pctN, na.rm = TRUE),
            uncolonized_roots.pctN = mean(uncolonized_roots.pctN, na.rm = TRUE))

nitrogen_onerow[is.nan.data.frame(nitrogen_onerow)] <- NA
# nitrospread = nitrogen %>% spread(tissue, APE15N)
# nitrospread$compid = paste(nitrospread$Plant, nitrospread$Side)

# for (i in 1:nrow(nitrospread)) {
#   if (nitrospread$mycofungus[i] == "MIXED") {
#     nitrospread$uncolonized_roots[nitrospread$compid %in% nitrospread$compid[i]] = carbonspread$uncolonized_roots[i]
#   }
# }

nitrospreadincludingmixed = nitrospread
# nitrospread = subset(nitrospread, mycofungus != "MIXED")

# nitrofinal = select(nitrospread, Plant, Side, hyphae.APE15N = hyphae,
#                     mycorrhizas.APE15N = mycorrhizas,
#                     uncolonized_roots.APE15N = uncolonized_roots,
#                     mycofungus)
# 
# nitrofinal_includingmixed = select(nitrospreadincludingmixed, Plant, Side, hyphae.APE15N = hyphae,
#                     mycorrhizas.APE15N = mycorrhizas,
#                     uncolonized_roots.APE15N = uncolonized_roots,
#                     mycofungus)

# nandc = full_join(carbonspread, nitrospread)
nandc = full_join(carbon_onerow, nitrogen_onerow)


onerow = select(onerowperplant, everything(), -APE13C, -APE15N, -pctC, -pctN, -tissue, -CNratio)
everything = left_join(nandc %>% group_by(Plant, Side, mycofungus) %>% 
                   mutate(id = row_number()), 
                 onerow %>% group_by(Plant, Side, mycofungus) %>% 
                   mutate(id = row_number()))
everything = select(everything, everything(), -id)

write_csv(everything, "processeddata/isotopes_two_rows_per_plant_updated_with_pctCN.csv")


# Old version 10/28/2020:
   everything = left_join(nandc, onerow)
   shouldbeunique = paste(everything$Plant, everything$Side, everything$mycofungus)
   sum(duplicated(shouldbeunique))
   everything = everything[!duplicated(shouldbeunique),]
nandc = full_join(carbonfinal, nitrofinal)
nandc_includingmixed = full_join(carbonfinal_includingmixed, nitrofinal_includingmixed)

onerow = select(onerowperplant, everything(), -APE13C, -APE15N)

everything = left_join(nandc, onerow)
everything_includingmixed = left_join(nandc_includingmixed, onerow)

shouldbeunique = paste(everything$Plant, everything$Side, everything$mycofungus)

everything = everything[!duplicated(shouldbeunique),]
everything_includingmixed = everything_includingmixed[!duplicated(shouldbeunique),]

# write_csv(everything, "isotopes_one_row_per_plant_including_unenriched.csv")
write_csv(everything, "processeddata/isotopes_one_row_per_plant_including_unenriched_July.csv")
write_csv(everything_includingmixed, "processeddata/isotopes_one_row_per_plant_including_unenriched_and_mixed.csv")

# UPDATE the below is great for N analyses and C for N,
# but actually I don't need that side to have received 15N
# for carbon calculations.
# actualonerow = subset(everything, received15N == "Y")
actualonerow = everything
actualonerow$CforN = actualonerow$mycorrhizas.APE13C/actualonerow$uncolonized_roots.APE15N
# But actually it's more subtle than this, right? 
# If the compartment is uncolonized, C for N should
# just reflect C in uncolonized roots vs N in uncolonized roots.

for (i in 1:nrow(actualonerow)) {
  if (actualonerow$Actual_fungus_by_compartment[i] == "NM") {
    actualonerow$CforN[i] = actualonerow$uncolonized_roots.APE13C[i]/actualonerow$uncolonized_roots.APE15N[i]
  } else if (actualonerow$Actual_fungus_by_compartment[i] != "NM") {
    actualonerow$CforN[i] = actualonerow$mycorrhizas.APE13C[i]/actualonerow$uncolonized_roots.APE15N[i]
  }
}

mynas = actualonerow[is.na(actualonerow$CforN),]
# I guess these are the incomplete cases, for now.
# What a bummer!
allCforNs = actualonerow
# allCforNs = subset(actualonerow, is.na(CforN) == FALSE)
# I included the below originally, but actually I don't think
# I need to a priori elminiate ALL mixed compartments.
# allCforNs = allCforNs[-grep("MIXED", allCforNs$Actual_fungi_at_harvest),]

allCforNs = subset(allCforNs, !is.na(allCforNs$Actual_fungi_at_harvest))

allCforNs$competition_treatment = numeric(nrow(allCforNs))
for (i in 1:nrow(allCforNs)) {
  if (allCforNs$Actual_fungus_by_compartment[i] == "THETE") {
    if (allCforNs$Actual_fungi_at_harvest[i] == "THETE/THETE") {
      allCforNs$competition_treatment[i] = "SELF"
    } else if (allCforNs$Actual_fungi_at_harvest[i] == "THETE/SUIPU") {
      allCforNs$competition_treatment[i] = "OTHER"
    } else if (allCforNs$Actual_fungi_at_harvest[i] == "THETE/NM") {
      allCforNs$competition_treatment[i] = "NM"
    }
  } else if (allCforNs$Actual_fungus_by_compartment[i] == "SUIPU") {
    if (allCforNs$Actual_fungi_at_harvest[i] == "SUIPU/SUIPU") {
      allCforNs$competition_treatment[i] = "SELF"
    } else if (allCforNs$Actual_fungi_at_harvest[i] == "THETE/SUIPU") {
      allCforNs$competition_treatment[i] = "OTHER"
    } else if (allCforNs$Actual_fungi_at_harvest[i] == "SUIPU/NM") {
      allCforNs$competition_treatment[i] = "NM"
    }
  }
}

allCforNs$mycoC13forN15 = allCforNs$mycorrhizas.APE13C/allCforNs$mycorrhizas.APE15N

# write_csv(allCforNs, "isotopwrite_csv(allCforNs, "isotope_data_two_rows_per_plant.csv")
write_csv(allCforNs, "processeddata/isotope_data_two_rows_per_plant_July.csv")

forNanalysis = subset(allCforNs, received15N == "Y")

# write_csv(forNanalysis, "isotope_data_one_row_per_plant.csv")
write_csv(forNanalysis, "processeddata/isotope_data_one_row_per_plant_July.csv")

#### JUST EXPLORATORY PLOTS FROM HERE ON OUT ####



allCforNs_noNM = subset(allCforNs, Actual_fungus_by_compartment != "NM")

ggplot(data = allCforNs) +
  geom_boxplot(aes(x = Actual_fungi_at_harvest, y = (CforN))) +
  geom_jitter(aes(x = Actual_fungi_at_harvest, 
                  y = CforN, 
                  color = Actual_fungus_by_compartment))

# CforNbyN_level = ggplot(data = allCforNs) +
#   geom_boxplot(aes(x = N_level, 
#                    y = CforN,
#                    color = Actual_fungus_by_compartment)) +
#   geom_jitter(aes(x = N_level, 
#                   y = CforN, 
#                   color = Actual_fungus_by_compartment))
# 
# pdf("plots/boxplot_C_for_N_by_N_level.pdf", width = 7, height = 5)
# CforNbyN_level
# dev.off()

# Let's look at this more closely. I'm not interested
# in the NM roots.


justmycos = ggplot(data = allCforNs) +
  geom_boxplot(aes(x = N_level, y = C13forN15,
                   color = Actual_fungus_by_compartment))
  # geom_point(aes(x = N_level, 
  #                 y = C13forN15, 
  #                 color = Actual_fungus_by_compartment))

pdf("plots/boxplot_C13_for_N15_mycos_by_N_level.pdf", width = 7, height = 5)
justmycos
dev.off()

t.test(C13forN15 ~ N_level, data = allCforNs)

test = allCforNs[is.na(allCforNs$C13forN15) == FALSE,]

ggplot(data = allCforNs) +
  geom_point(aes(x =uncolonized_roots.APE15N,
                  y = mycorrhizas.APE13C, 
                  color = Actual_fungus_by_compartment,
                  shape = N_level)) +
  geom_smooth(method = "lm", aes(x =uncolonized_roots.APE15N,
                                 y = mycorrhizas.APE13C,
                  color = Actual_fungus_by_compartment))

max(allCforNs$mycorrhizas.APE13C[!is.na(allCforNs$mycorrhizas.APE13C)])
allCforNs[allCforNs$mycorrhizas.APE13C > 0.12,] #6041b

# What if I remove the outlier?
allCforNs = read_csv("isotope_data_one_row_per_plant.csv")
nooutlier = subset(allCforNs, Plant != 6041)
nooutlier_noNMs = subset(allCforNs_noNM, Plant != 6041)
ggplot(data = nooutlier) +
  geom_point(aes(x =mycorrhizas.APE15N,
                 y = mycorrhizas.APE13C, 
                 color = Actual_fungus_by_compartment,
                 shape = N_level)) +
  geom_smooth(method = "lm", aes(x =uncolonized_roots.APE15N,
                                 y = mycorrhizas.APE13C,
                                 color = Actual_fungus_by_compartment))

ggplot(data = nooutlier) +
  geom_point(aes(x =uncolonized_roots.APE15N,
                 y = mycorrhizas.APE13C, 
                 color = Actual_fungus_by_compartment,
                 shape = N_level)) +
  geom_smooth(method = "lm", aes(x =uncolonized_roots.APE15N,
                                 y = mycorrhizas.APE13C,
                                 color = Actual_fungus_by_compartment))
ggplot(data = nooutlier) +
  geom_point(aes(x = mycorrhizas.APE13C,
                 y = hyphae.APE13C, 
                 color = Actual_fungus_by_compartment,
                 shape = N_level)) +
  geom_smooth(method = "lm", aes(x = mycorrhizas.APE13C,
                                 y = hyphae.APE13C, 
                                 color = Actual_fungus_by_compartment))

myhyphaelm = lm(hyphae.APE13C ~ mycorrhizas.APE13C, data = allCforNs)
summary(myhyphaelm)

myhyphaelm_withN = lm(hyphae.APE13C ~ mycorrhizas.APE13C*N_level, data = allCforNs)
summary(myhyphaelm_withN) # Oh my, that's much worse.

# What if I try including ALL the samples?
# Maybe I have some hyphae-myco pairs on the non-focal sides.

forplot = subset(everything,
                 Actual_fungus_by_compartment == "SUIPU" |
                   Actual_fungus_by_compartment == "THETE")
forplot = subset(forplot, is.na(forplot$hyphae.APE13C) == FALSE)
ggplot(data = forplot) +
  geom_point(aes(x = mycorrhizas.APE13C,
                 y = hyphae.APE13C, 
                 color = Actual_fungus_by_compartment,
                 shape = N_level)) +
  geom_smooth(method = "lm", aes(x = mycorrhizas.APE13C,
                                 y = hyphae.APE13C))

fulldatahyphaelm = lm(hyphae.APE13C ~ mycorrhizas.APE13C, data = forplot)
summary(fulldatahyphaelm)

fulldatahyphaelm_withfungus = lm(hyphae.APE13C ~ mycorrhizas.APE13C*Actual_fungus_by_compartment, data = forplot)
summary(fulldatahyphaelm_withfungus)
# Oh my, this model is way worse if we add
# fungal species as a predictor.

ggplot(data = nooutlier) +
  geom_point(aes(x = mycorrhizas.APE15N,
                 y = hyphae.APE15N, 
                 color = Actual_fungus_by_compartment,
                 shape = N_level)) +
  geom_smooth(method = "lm", aes(x = mycorrhizas.APE15N,
                                 y = hyphae.APE15N))

myhyphaeNlm = lm(hyphae.APE15N ~ mycorrhizas.APE15N, data = allCforNs)
summary(myhyphaeNlm)

ggplot(data = nooutlier) +
  geom_point(aes(x = hyphae.APE15N,
                 y = uncolonized_roots.APE15N, 
                 color = Actual_fungus_by_compartment,
                 shape = N_level)) +
  geom_smooth(method = "lm", aes(x = hyphae.APE15N,
                                 y = uncolonized_roots.APE15N))

Ntoplantlm = lm(uncolonized_roots.APE15N ~ hyphae.APE15N, data = nooutlier)
summary(Ntoplantlm)

ggplot(data = nooutlier) +
  geom_point(aes(x = mycorrhizas.APE15N,
                 y = uncolonized_roots.APE15N, 
                 color = Actual_fungus_by_compartment,
                 shape = N_level)) +
  geom_smooth(method = "lm", aes(x = mycorrhizas.APE15N,
                                 y = uncolonized_roots.APE15N))

Nmycotoplantlm = lm(uncolonized_roots.APE15N ~ mycorrhizas.APE15N, data = nooutlier)
summary(Nmycotoplantlm)

# would I see clearer signal of C and N corresponding
# in mycorrhizas if I controlled for competition treatment?

subtoplot = subset(nooutlier, Actual_fungi_at_harvest == "THETE/THETE")
ggplot(data = subtoplot) +
  geom_point(aes(x = mycorrhizas.APE15N,
                 y = mycorrhizas.APE13C, 
                 color = Actual_fungus_by_compartment,
                 shape = N_level)) +
  geom_smooth(method = "lm", aes(x =mycorrhizas.APE15N,
                                 y = mycorrhizas.APE13C))
# Wow, that looks WAAAY better

CforNlm_THETETHETE = lm(mycorrhizas.APE13C ~ mycorrhizas.APE15N, data = subtoplot)
summary(CforNlm_THETETHETE)

# How about THETE/NM?
THETENM = subset(nooutlier, Actual_fungi_at_harvest == "THETE/NM")
ggplot(data = THETENM) +
  geom_point(aes(x = mycorrhizas.APE15N,
                 y = mycorrhizas.APE13C,
                 shape = N_level)) +
  geom_smooth(method = "lm", aes(x =mycorrhizas.APE15N,
                                 y = mycorrhizas.APE13C))
CforNlm_THETENM = lm(mycorrhizas.APE13C ~ mycorrhizas.APE15N, data = THETENM)
summary(CforNlm_THETENM)

# How about THETE/SUIPU?
THETESUIPU = subset(nooutlier, Actual_fungi_at_harvest == "THETE/SUIPU")
ggplot(data = THETESUIPU) +
  geom_point(aes(x = mycorrhizas.APE15N,
                 y = mycorrhizas.APE13C,
                 shape = N_level,
                 color = Actual_fungus_by_compartment)) +
  geom_smooth(method = "lm", aes(x =mycorrhizas.APE15N,
                                 y = mycorrhizas.APE13C))

CforNlm_THETESUIPU = lm(mycorrhizas.APE13C ~ mycorrhizas.APE15N, data = THETESUIPU)
summary(CforNlm_THETESUIPU)

CforNeverything = (lm(mycorrhizas.APE13C ~ mycorrhizas.APE15N*Actual_fungi_at_harvest, data = nooutlier))
summary(CforNeverything)

CforNnocompinfo = (lm(mycorrhizas.APE13C ~ mycorrhizas.APE15N, data = nooutlier))
summary(CforNnocompinfo)

#### Mycorrhiza 13C vs mycorrhiza 15N ####

# This curve no longer looks good when it's possible 
# to have negative APE15N values, because
# the log calculation breaks.

forcurveplot = allCforNs_noNM
mymin = min(forcurveplot$mycorrhizas.APE15N)
# Coerce all values to 0.0000001 or greater -- smaller than
# min(abs(allCforNs_noNM$mycorrhizas.APE15N))
forcurveplot$mycorrhizas.APE15N = allCforNs_noNM$mycorrhizas.APE15N + abs(mymin) + 0.0000001

curveplot = ggplot(data = forcurveplot) +
  geom_point(aes(x = mycorrhizas.APE15N,
                 y = mycorrhizas.APE13C, 
                 color = Actual_fungus_by_compartment,
                 shape = N_level)) +
  geom_smooth(method = "glm", 
              formula = y ~ log(x), 
              aes(x =mycorrhizas.APE15N,
                                 y = mycorrhizas.APE13C))

pdf("plots/myco_c_for_n_withcurve.pdf", width = 7, height = 5)
curveplot
dev.off()

forcurveplot_nooutlier = subset(forcurveplot, Plant != 6041)

curveplot_nooutlier = ggplot(data = forcurveplot_nooutlier) +
  geom_point(aes(x = mycorrhizas.APE15N,
                 y = mycorrhizas.APE13C, 
                 color = Actual_fungus_by_compartment,
                 shape = N_level)) +
  geom_smooth(method = "glm", 
              formula = y ~ log(x), 
              aes(x =mycorrhizas.APE15N,
                  y = mycorrhizas.APE13C))

pdf("plots/myco_c_for_n_withcurve_nooutlier.pdf", width = 7, height = 5)
curveplot_nooutlier
dev.off()


mymodel = glm(mycorrhizas.APE13C ~ log(mycorrhizas.APE15N), data = forcurveplot_nooutlier)
summary(mymodel)

linearmodel = glm(mycorrhizas.APE13C ~ (mycorrhizas.APE15N), data = forcurveplot_nooutlier)
summary(linearmodel)

# forcurveplot_nooutlier$propcol = forcurveplot_nooutlier$percent_col/100
# bionomialmodel = glm(mycorrhizas.APE13C ~ )

ggplot(data = nooutlier) +
  geom_point(aes(x = uncolonized_roots.APE15N,
                 y = mycorrhizas.APE13C, 
                 color = Actual_fungus_by_compartment,
                 shape = N_level)) +
  geom_smooth(method = "lm", aes(x = uncolonized_roots.APE15N,
                                 y = mycorrhizas.APE13C))
# 
# CforNagain = (lm(mycorrhizas.APE13C ~ mycorrhizas.APE15N, data = nooutlier))
# summary(CforNagain)
# 
# anothertest = subset(nooutlier, Actual_fungus_by_compartment != "NM")
# testwithKP = lm(mycorrhizas.APE13C ~ mycorrhizas.APE15N*Actual_fungus_by_compartment, data = anothertest)
# summary(testwithKP)
# 
# THETEfocus = subset(nooutlier, Actual_fungus_by_compartment == "THETE")
# summary(as.factor(THETEfocus$Actual_fungi_at_harvest))
# 
# CforNTHETEonly = lm(mycorrhizas.APE13C ~ mycorrhizas.APE15N*Actual_fungi_at_harvest, data = THETEfocus)
# summary(CforNTHETEonly)
# 
# CforNTHETEnocomp = lm(mycorrhizas.APE13C ~ mycorrhizas.APE15N, data = THETEfocus)
# summary(CforNTHETEnocomp)
# 
# CforNTHETEvsroots = lm(mycorrhizas.APE13C ~ uncolonized_roots.APE15N, data = THETEfocus)
# summary(CforNTHETEvsroots)
# 
# CforNhyphaevsroots = lm(hyphae.APE13C ~ uncolonized_roots.APE15N, data = THETEfocus)
# summary(CforNhyphaevsroots) # Wow, REALLY bad.



toplot = subset(allCforNs, competition_treatment != 0)

comptx = ggplot(data = toplot) +
  geom_boxplot(aes(x = competition_treatment, y = C13forN15)) +
  geom_jitter(aes(x = competition_treatment, 
                  y = C13forN15, 
                  color = Actual_fungus_by_compartment,
                  shape = N_level)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf("plots/boxplot_C13_for_N15_mycos_by_competition_treatment.pdf", width = 7, height = 5)
comptx
dev.off()

atest = aov(C13forN15 ~ competition_treatment, data = toplot)
summary(atest)

ggplot(data = allCforNs) +
  geom_boxplot(aes(x = N_level, y = C13forN15,
                   color = Actual_fungus_by_compartment))

amodel = glm(C13forN15 ~ N_level*Actual_fungus_by_compartment, data = allCforNs)
summary(amodel)
myanova = aov(C13forN15 ~ N_level*Actual_fungus_by_compartment, data = allCforNs)
summary(myanova)
TukeyHSD(myanova)
# Hmm, EVERYTHING is significant!
# N level, fungal identity, and interaction.
#

fungalidtradingrate = glm(data = allCforNs)
