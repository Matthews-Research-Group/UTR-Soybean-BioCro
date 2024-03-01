result <- soybean_optsolver[[i]](optim_params_conversion(optim_params8))

# organize simulated data
r.lsrp.doy <- reshape2::melt(result[,c("time","Root","Leaf","Stem","Pod")],id.vars="time")
r.lsrp.doy$value<-r.lsrp.doy$value

# Leaf
# organize the expirimental data (mean and std)
s.exp.leaf <- cbind(ExpBiomass[[i]][,c("DOY","Leaf")])
colnames(s.exp.leaf) <- c("time","Leaf") # DOY renamed as time
r.exp.leaf <- reshape2::melt(s.exp.leaf, id.vars = "time") # DOY renamed as time

s.exp.std.leaf <- cbind(ExpBiomass.std[[i]][,c("DOY","Leaf")])
colnames(s.exp.std.leaf) <- c("time","Leaf") # DOY renamed as time
r.exp.std.leaf <- reshape2::melt(s.exp.std.leaf, id.vars = "time") # DOY renamed as time
r.exp.std.leaf$ymin <- r.exp.leaf$value - r.exp.std.leaf$value
r.exp.std.leaf$ymax <- r.exp.leaf$value + r.exp.std.leaf$value

# Stem
s.exp.stem <- cbind(ExpBiomass[[i]][,c("DOY","Stem")]) # DOY renamed as time
colnames(s.exp.stem) <- c("time","Stem") # DOY renamed as time
r.exp.stem <- reshape2::melt(s.exp.stem, id.vars = "time") # DOY renamed as time

s.exp.std.stem <- cbind(ExpBiomass.std[[i]][,c("DOY","Stem")]) # DOY renamed as time
colnames(s.exp.std.stem) <- c("time","Stem") # DOY renamed as time
r.exp.std.stem <- reshape2::melt(s.exp.std.stem, id.vars = "time") # DOY renamed as time
r.exp.std.stem$ymin <- r.exp.stem$value - r.exp.std.stem$value
r.exp.std.stem$ymax <- r.exp.stem$value + r.exp.std.stem$value

# Pod
s.exp.pod <- cbind(ExpBiomass[[i]][,c("DOY","Pod")]) # DOY renamed as time
colnames(s.exp.pod) <- c("time","Pod") # DOY renamed as time
r.exp.pod <- reshape2::melt(s.exp.pod, id.vars = "time") # DOY renamed as time

s.exp.std.pod <- cbind(ExpBiomass.std[[i]][,c("DOY","Pod")]) # DOY renamed as time
colnames(s.exp.std.pod) <- c("time","Pod") # DOY renamed as time
r.exp.std.pod <- reshape2::melt(s.exp.std.pod, id.vars = "time") # DOY renamed as time
r.exp.std.pod$ymin <- r.exp.pod$value - r.exp.std.pod$value
r.exp.std.pod$ymax <- r.exp.pod$value + r.exp.std.pod$value

# Combine
r.exp.ls <- rbind(r.exp.leaf, r.exp.stem, r.exp.pod)
r.exp.ls$source = "Observed"
r.lsrp.doy$source = "Simulated"
r.all <- rbind(r.lsrp.doy, r.exp.ls)
# Reverse the order as follow
r.all$variable <- factor(r.all$variable, levels = rev(levels(r.all$variable)))

# combine the simulated and experimental data

# Colorblind friendly color palette (https://personal.sron.nl/~pault/)
col.palette.muted <- c("#332288", "#117733", "#999933", "#882255")

size.title <- 16
size.axislabel <-14
size.axis <- 16
size.legend <- 14

f <- ggplot() + theme_classic()

f <- f + geom_point(data=r.all, aes(x=time, y=value,
                                    color=variable,
                                    size = source, shape = source),
                    show.legend = TRUE, stroke=0.5) +
  scale_shape_manual(values = c(15, 16)) +
  scale_size_manual(values = c(4, 0.5)) +
  scale_color_manual(values = col.palette.muted)

# for leaf
f <- f + geom_errorbar(data=r.exp.std.leaf, aes(x=time, ymin=ymin, ymax=ymax),  # DOY renamed as time
                       width=3.5, size=0.25, show.legend = FALSE)
# for stem
f <- f + geom_errorbar(data=r.exp.std.stem, aes(x=time, ymin=ymin, ymax=ymax),   # DOY renamed as time
                       width=3.5, size=0.25, show.legend = FALSE)

# for pod
f <- f + geom_errorbar(data=r.exp.std.pod, aes(x=time, ymin=ymin, ymax=ymax),   # DOY renamed as time
                       width=3.5, size=0.25, show.legend = FALSE)

# change the plot labels and theme
f <- f + labs(title=element_blank(), x=paste0('Day of Year (',year[i],')'),y='Biomass (Mg / ha)')
f <- f + theme(plot.title=element_text(size=size.title, hjust=0.5),
               axis.text=element_text(size=size.axis),
               axis.title=element_text(size=size.axislabel),
               legend.position = c(.1,.85), legend.title = element_blank(),
               legend.text=element_text(size=size.legend),
               legend.background = element_rect(fill = "transparent",colour = NA),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA))

f <- f + scale_x_continuous(breaks = seq(150,280,30))
f <- f + guides(shape = guide_legend(override.aes = list(size = 2)))
f <- f + guides(color = guide_legend(override.aes = list(size=2)))
f

# Compare mass fractions = substrate carbon / structural carbon
key <- list(text = c('Leaf', 
                     'Stem',
                     'Root',
                     'Pod'),
            x = .8, y = .9)
total_hours = dim(result)[1]
xyplot(result$Leaf_mass_fraction[1:(total_hours-24*20)]
       +result$Stem_mass_fraction[1:(total_hours-24*20)]
       +result$Root_mass_fraction[1:(total_hours-24*20)]
       +result$Pod_mass_fraction[1:(total_hours-24*20)]
       ~result$time[1:(total_hours-24*20)],
       main = 'Mass Fractions of Organs', 
       xlab = 'Time (Day of Year)',
       ylab = 'Mass Fraction = Substrate C:Structural C',
       auto.key = key)

# Compare mass fractions = substrate carbon / Dry Mass
key <- list(text = c('Leaf', 
                     'Stem',
                     'Root',
                     'Pod'),
            x = .8, y = .9)
xyplot(result$Leaf_substrate_carbon/result$Leaf
       +result$Stem_substrate_carbon/result$Stem
       +result$Root_substrate_carbon/result$Root
       +result$Pod_substrate_carbon/result$Pod
       ~result$time,
       main = 'Mass Fractions of Organs', 
       xlab = 'Time (Day of Year)',
       ylab = 'Mass Fraction = Substrate C/Dry Mass (mol/gDM-1)',
       auto.key = key)

# Compare structural carbon
c_key <- list(text = c('Leaf', 
                       'Stem',
                       'Root',
                       'Pod'),
              x = .05, y = .9)
xyplot(result$Leaf_structural_carbon
       +result$Stem_structural_carbon
       +result$Root_structural_carbon
       +result$Pod_structural_carbon
       ~result$time, 
       xlab = 'DOY',
       ylab = 'Structural Carbon (mol/m^2)',
       main = 'Structural Carbon in Organs', auto.key=c_key)
# Compare substrate carbon
xyplot(result$Leaf_substrate_carbon
       +result$Stem_substrate_carbon
       +result$Root_substrate_carbon
       +result$Pod_substrate_carbon
       ~result$time, 
       xlab = 'DOY',
       ylab = 'Substrate Carbon (mol/m^2)',
       main = 'Substrate Carbon in Organs', auto.key=c_key)
# Compare substrate carbon
xyplot(result$Leaf_substrate_carbon / result$lai
       ~result$time, 
       xlab = 'DOY',
       ylab = 'Substrate Carbon (mol/m^2)',
       main = 'Substrate Carbon in Organs', auto=TRUE)

# Compare Transport Rate in All Links
t_key <- list(text = c('Leaf to Stem', 
                       'Stem to Root',
                       'Stem to Pod'),
              x = .8, y = .9)
xyplot(result$substrate_transport_Leaf_to_Stem
       +result$substrate_transport_Stem_to_Root
       +result$substrate_transport_Stem_to_Pod
       ~result$time,
       xlab = 'DOY',
       ylab = 'Transport rate (mol/m^2/hr)',
       main = 'Transport Rates', auto.key = t_key)

# Compare Leaf Assimilation Rate, Export, Utilization rates, and percentages
xyplot(result$canopy_assimilation_rate*0.6/0.18
       +result$Leaf_utilization_rate
       +result$substrate_transport_Leaf_to_Stem
       +result$Leaf_senescence_rate
       ~result$time,
       main = 'Leaf Input and Outputs Snapshot', 
       xlab = 'DOY',
       ylab = 'Rate (mol/m^2/hr)',
       auto.key = list(text = c('Canopy Assimilation', 
                                'Utilization',
                                'Transport Leaf to Stem',
                                'Senescence'),
                       x = .8, y = .98))

start.time = 24*30+1
end.time = start.time+24*5
xyplot((result$canopy_assimilation_rate*0.6/0.18)[start.time:end.time]
       +result$Leaf_utilization_rate[start.time:end.time]
       +result$substrate_transport_Leaf_to_Stem[start.time:end.time]
       +result$Leaf_senescence_rate[start.time:end.time]
       ~result$time[start.time:end.time],
       main = 'Leaf Input and Outputs Snapshot', 
       xlab = 'DOY',
       ylab = 'Rate (mol/m^2/hr)',
       par.settings = list(superpose.symbol = list(pch = 19, cex = 1,
                                                   col = c("blue", "magenta", "green", "red"))),
       auto.key = list(text = c('Canopy Assimilation', 
                                'Utilization',
                                'Transport Leaf to Stem',
                                'Senescence'),
                       x = .01, y = 0.99))

print("Daily Net Assimilation")
tot_net_assim = sum((result$canopy_assimilation_rate*0.6/0.18)[start.time:end.time])
print(tot_net_assim/5)
total_util = sum(result$Leaf_utilization_rate[start.time:end.time])
print("Daily Utilization")
print(total_util/5)
print("Toal Utilization percentage")
print(total_util/tot_net_assim)

total_trans = sum(result$substrate_transport_Leaf_to_Stem[start.time:end.time])
print("Toal Export percentage")
print(total_trans/tot_net_assim)

total_sene = sum(result$Leaf_senescence_rate[start.time:end.time])
print("Toal Senescence percentage")
print(total_sene/tot_net_assim)

xyplot(result$Leaf_mass_fraction[start.time:end.time]
       +result$Stem_mass_fraction[start.time:end.time]
       ~result$time[start.time:end.time],
       main = 'Leaf Mass Fraction Snapshot', 
       xlab = 'DOY',
       ylab = 'Rate (mol/m^2/hr)',
       auto.key = list(text = c('Leaf Mass Fraction',
                                'Stem Mass Fraction'),
                       x = .01, y = .98))

transport_percent = result$substrate_transport_Leaf_to_Stem/result$canopy_assimilation_rate*(60/18)
utilization_percent = result$Leaf_utilization_rate/result$canopy_assimilation_rate*(60/18)
xyplot(transport_percent + utilization_percent ~ result$time, ylim=c(-1,5), auto=TRUE)
# Need to look at daily sum

# Compare Stem Import, Export, Utilization rates, and percentages
xyplot(result$substrate_transport_Leaf_to_Stem
       +result$Stem_utilization_rate
       +result$substrate_transport_Stem_to_Root
       +result$substrate_transport_Stem_to_Pod
       +result$Stem_senescence_rate~result$time, 
       main = 'Stem Input and Outputs', 
       xlab = 'DOY',
       ylab = 'Rate (mol/m^2/hr)',
       auto.key = list(text = c('Transport Leaf to Stem', 
                                'Utilization',
                                'Transport Stem to Root',
                                'Transport Stem to Pod',
                                'Senescence'),
                       x = .7, y = .98))

# Compare Root Import, Utilization rates, and percentages
xyplot(result$substrate_transport_Stem_to_Root
       +result$Root_utilization_rate
       +result$Root_senescence_rate~result$time, 
       main = 'Root Input and Utilization', 
       xlab = 'DOY',
       ylab = 'Rate (mol/m^2/hr)',
       auto.key = list(text = c('Transport Stem to Root', 
                                'Utilization',
                                'Senescence'),
                       x = .05, y = .95))

# Compare Pod Import, Utilization rates, and percentages
xyplot(result$assimilation_rate
       +result$substrate_transport_Stem_to_Pod
       +result$Pod_utilization_rate
       +result$Pod_senescence_rate~result$time, 
       main = 'Pod Input and Utilization', 
       xlab = 'DOY',
       ylab = 'Rate (mol/m^2/hr)',
       auto.key = list(text = c('Leaf Assimilation rate',
                                'Transport Stem to Pod', 
                                'Utilization',
                                'Senescence'),
                       x = .05, y = .95))

# Plot senescence
s_key <- list(text = c('Leaf', 
                       'Stem',
                       'Sum',
                       'Observation'),
              x = .05, y = .95)
total_senescence = result$Leaf_senescence_loss*1/3+result$Stem_senescence_loss*1/3
biocro_senescence_plot <- xyplot(result$Leaf_senescence_loss*1/3
                                 +result$Stem_senescence_loss*1/3
                                 +total_senescence
                                 ~result$time, 
                                 main = 'Senescence Loss/Litter', xlab = 'DOY', ylab='Senescence Loss (Mg/ha)', 
                                 ylim = c(0,5), auto.key=s_key)
experiment_plot <- xyplot(c(0,   0,   0,   0,   0.38,0.66,0.98,2.58)~
                            c(179, 189, 203, 217, 231, 246, 259, 288),cex = 2,col = "red")
library(latticeExtra)
biocro_senescence_plot+as.layer(experiment_plot)
biocro_senescence_plot+as.layer(experiment_plot)

