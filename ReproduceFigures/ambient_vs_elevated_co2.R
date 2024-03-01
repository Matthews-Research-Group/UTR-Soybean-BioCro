library(ggplot2)
year <- c('2002', '2005', '2004', '2006')
size.title <- 24
size.axislabel <-18
size.axis <- 20
size.legend <- 20

for (i in 1:length(year)) { 
  yr = year[i]
  
  data_ambient <- read.csv(file=paste0('../Data/SoyFACE_data/',yr, '_ambient_', 'biomass.csv'))
  colnames(data_ambient)<-c("DOY","Leaf","Stem","Pod")
  data_ambient <- reshape2::melt(data_ambient, id.vars = "DOY")
  data_ambient$co2_opt = 'Ambient CO2 (372 ppm)'
  
  data_co2 <- read.csv(file=paste0('../Data/SoyFACE_data/',yr, '_co2_', 'biomass.csv'))
  colnames(data_co2)<-c("DOY","Leaf","Stem","Pod")
  data_co2 <- reshape2::melt(data_co2, id.vars = "DOY")
  data_co2$co2_opt = 'Elevated CO2 (550 ppm)'
  
  # Combine data set and reshape
  data_combined <- rbind(data_ambient, data_co2)
  
  # Plot
  f <- ggplot() + theme_classic()
  f <- f + theme(plot.title=element_text(size=size.title, hjust=0.5),
                 axis.text=element_text(size=size.axis),
                 axis.title=element_text(size=size.axislabel),
                 legend.position = c(.25,.85), legend.title = element_blank(),
                 legend.text=element_text(size=size.legend),
                 legend.background = element_rect(fill = "transparent",colour = NA),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
                 plot.background = element_rect(fill = "transparent", colour = NA))
  f <- f + geom_line(data=data_combined, aes(x=DOY, y=value,
                                             color=variable,
                                             linetype = co2_opt))
  f <- f + geom_point(data=data_combined, aes(x=DOY, y=value,
                                              color=variable,
                                              shape = co2_opt))
  f <- f + labs(title=element_blank(), x=paste0('Day of Year (',yr,')'),y='Biomass (Mg / ha)')
  ggsave(paste0(yr , 'biomass_ambient_vs_elevated_co2.png'))
  print(f)
}