# Clear workspace
rm(list=ls())
# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Load Data
load('SoyFACE_results_and_measurements.RData')
TNC.data <- read.csv('../Data/2022_Carb_data/2022_LD11_TNC.csv')

# years
years <- c('2002','2004','2005','2006')

result <- results[[1]]

# Plot simulated TNC per m2 
sim_subC_per_m2 <- data.frame(
  time = result$time,
  DOY = result$doy,
  Leaf = result$Leaf_substrate_carbon,
  Stem = result$Stem_substrate_carbon
)

sim_subC_per_m2_long <- tidyr::gather(sim_subC_per_m2[1:2600,c('time',
                                                               'Leaf', 'Stem')], 
                                      key="Type", value="Value", -time)
sim_subC_per_m2_long$DOY <- as.integer(sim_subC_per_m2$time[1:2600])
View(sim_subC_per_m2_long)

ggplot(sim_subC_per_m2_long, aes(time, Value, group = Type)) + 
  geom_point(aes(color = Type, size = Type))+
  scale_shape_manual(values=c(18, 16)) +
  scale_size_manual(values=c(1, 0.5)) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  # scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  scale_x_continuous(breaks = seq(180,280,30))+
  labs(title=element_blank(), 
       x='Day of Year (2022)',
       y='Leaf Substrate C (mol C/m^2)')


# Plot measured vs simulated TNC per mass 
# Simulated alone
sim_substrate_C_by_mass <- data.frame(
  time = result$time,
  Leaf = 10^4 * result$Leaf_substrate_carbon / result$Leaf,
  Stem = 10^4 * result$Stem_substrate_carbon / result$Stem
)

data_long_per_mass <- tidyr::gather(sim_substrate_C_by_mass[1:2600,c('time', 'Leaf', 
                                                                     'Stem')], 
                                    key="Type", value="Value", -time)

# Create line plot
ggplot(data_long_per_mass, aes(x=time, y=Value/1000, color=Type)) +
  theme_classic() +
  scale_color_discrete(labels = c("Leaf", "Stem"))+
  theme(legend.position = "bottom", 
        panel.background = element_rect(fill = "transparent",colour = NA))+
  geom_line()+
  labs(title = 'Substrate C per mass', 
       x='Day of Year (2022)',
       y='Substrate C (mol C/ kg)')

# Simulated + Measured
# Leaf
sim_substrate_C_by_mass$Source <- 'Simulated'

Leaf.carb.data <- rbind(sim_substrate_C_by_mass[,c('time','Leaf','Source')], 
                        TNC.data[, c('time','Leaf','Source')])

ggplot(Leaf.carb.data, aes(time, Leaf, group = Source)) + 
  geom_point(aes(shape=Source, color=Source, size=Source))+
  scale_shape_manual(values=c(18, 16)) +
  scale_size_manual(values=c(2, 0.5)) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  #scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 1)) +
  scale_x_continuous(breaks = seq(180,280,30))+
  labs(title=element_blank(), 
       x='Day of Year (2022)',
       y='Leaf Substrate C (mol C/Mg)')

# Stem
Stem.carb.data <- rbind(sim_substrate_C_by_mass[,c('time','Stem','Source')], 
                        TNC.data[, c('time','Stem','Source')])

ggplot(Stem.carb.data, aes(time, Stem, group = Source)) + 
  geom_point(aes(shape=Source, color=Source, size=Source))+
  scale_shape_manual(values=c(18, 16)) +
  scale_size_manual(values=c(2, 0.5)) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  # scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  scale_x_continuous(breaks = seq(180,280,30))+
  labs(title=element_blank(), 
       x='Day of Year (2022)',
       y='Stem Substrate C (mol C/Mg)')
