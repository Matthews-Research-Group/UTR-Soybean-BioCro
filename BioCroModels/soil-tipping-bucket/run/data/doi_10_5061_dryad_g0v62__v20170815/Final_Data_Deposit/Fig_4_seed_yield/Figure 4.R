#Infile and format seed yield data
setwd("/Users/sharongray/Documents/Grad School Manuscripts/Soil Moisture and driFACE manuscript/Final Submission/Final/Fig_4_seed_yield")
library(ggplot2)
seed_yield <- read.csv("seed_yield_data_2004thru_2011_2010CPrm.csv", header = TRUE, sep = ',', stringsAsFactors=FALSE)
seed_yield$P.PET..m. <- as.numeric(seed_yield$P.PET..mm.)
seed_yield$seed.yield.CO2.response <- as.numeric(seed_yield$seed.yield.CO2.response)
seed_yield$Cultivar <- as.factor(seed_yield$Cultivar)
##############################################################FIGURE 4########################################
#plot ratio of seed yield in elevated CO2: seed yield in ambient CO2 against precipitation - potential evapotranspiration 
Fig4 <- ggplot(seed_yield, aes(x=P.PET..m., y=seed.yield.CO2.response), shape=Cultivar) +
  geom_point(aes(shape= Cultivar), size=5, colour="black", stroke=2) +
  geom_smooth(method = 'lm',size = 1, colour = 'black') + theme_bw() +
  theme_classic(base_size=25) +
  xlab(expression('Mean Daily P - PET (mm)')) +
  xlim(-3.5, -1.1) +
  ylab(expression('CO'[2]*' '*'Effect on Seed Yield')) +
  ylim(0.7, 1.5) +
  scale_shape_manual(values=c(0, 1, 2, 3, 4, 5, 6)) +
  coord_fixed(ratio=3) +
  theme(legend.margin = unit(0.5, "cm")) +
  theme(legend.key.size = unit(1, "cm")) +
  theme(axis.line = element_line(colour = "black", size = 2, linetype = "solid")) +
  theme(axis.text = element_text(colour = "black", size = 18)) +
  theme(axis.ticks = element_line(colour="black", size=rel(2)), axis.ticks.length = unit(0.5, "cm")) +
  theme(axis.title.y=element_text(vjust=0.5)) 
 Fig4 
 ggsave("Fig4_CI.pdf", width = 10, height = 8)

#Linear Regression Analysis: Regress the ratio of seed yield in elevated CO2: seed yield in ambient CO2 against precipitation - potential evapotranspiration
Fig4_reg = lm(seed.yield.CO2.response ~ P.PET..m., seed_yield)               
summary(Fig4_reg) 







