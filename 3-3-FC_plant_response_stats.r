# 3-3-FC_plant_response stats

#### PLANT RESPONSE TO COLONIZATION ####

forplot = subset(plantleveldata, Fungi != "None/None")


tx = with(forplot, interaction(N_level, Fungi))
anovaforplot = aov(plant_response ~ tx, data = forplot)

# from "agricolae" package
mylabels = HSD.test(anovaforplot, "tx", group = TRUE)

annotations = data.frame(x = c((1:5), (1:5)),
                         y = c(1, 1, 1.3, 1.3, 1.7, 1, 1, 1, 1, 1),
                         N_level = c(rep("High", 5), rep("Low", 5)),
                         labs = c(paste(c("bc", "abc", "ab", "bc", "a")), paste(c("c", "bc", "bc", "bc", "bc"))))


labels = c(High = "High N", Low = "Low N")

responseplot = ggplot(data = forplot) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = Fungi, y = plant_response)) +
  geom_jitter(width = 0.20,
              aes(x = Fungi, y = plant_response)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Plant response to fungal\ncolonization (log response ratio)") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_text(data = annotations, aes(x, y, label = labs)) +
  xlab("Fungal treatment")


pdf("plots/Plant_response_boxplot.pdf", width = 7, height = 5)
responseplot
dev.off()

responseanova = aov(plant_response ~ N_level * Fungi, data = forplot)
summary(responseanova)
responseTukey = TukeyHSD(responseanova)

write.csv(responseTukey$N_level, "Statistical_tables/Plant_response_Tukey_Nlevel.csv")
write.csv(responseTukey$Fungi, "Statistical_tables/Plant_response_Tukey_Fungi.csv")
write.csv(responseTukey$`N_level:Fungi`, "Statistical_tables/Plant_response_Tukey_Nlevel-Fungi.csv")

# Do these differ from zero?
summary(forplot$Fungi)

# High N:
hightx = subset(forplot, N_level == "High")
t.test(hightx$plant_response[hightx$Fungi == "Sp/None"]) # not enough observations
t.test(hightx$plant_response[hightx$Fungi == "Sp/Sp"]) # p = 0.1185
t.test(hightx$plant_response[hightx$Fungi == "Tt/Sp"]) # p = 0.02562
t.test(hightx$plant_response[hightx$Fungi == "Tt/None"]) # p = 0.1547
t.test(hightx$plant_response[hightx$Fungi == "Tt/Tt"]) # p = 7.478e-11

# Low N:
lowtx = subset(forplot, N_level == "Low")
t.test(lowtx$plant_response[lowtx$Fungi == "Sp/None"]) # p = 0.1654
t.test(lowtx$plant_response[lowtx$Fungi == "Sp/Sp"]) # p = 0.6937
t.test(lowtx$plant_response[lowtx$Fungi == "Tt/Sp"]) # p = 0.2881
t.test(lowtx$plant_response[lowtx$Fungi == "Tt/None"]) # p = 0.399
t.test(lowtx$plant_response[lowtx$Fungi == "Tt/Tt"]) # p = 0.1192

# Alpha should be 0.05/10 tests = 0.005, per Bonferroni correction. So the only significant 
# non-zero growth effect difference is Tt/Tt in high N




