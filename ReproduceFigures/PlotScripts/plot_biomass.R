size.title <- 12
size.axislabel <- 12
size.axis <- 12
size.legend <- 12

plot_SoyFACE_biomass <- function(result, mea, mea.std, co2_opt, yr){
  # organize simulated data
  r.lsrp.doy <- reshape2::melt(result[,c("fractional_doy","Root","Leaf","Stem","Pod")],id.vars="fractional_doy")
  
  # Leaf
  # organize the experimental data (mean and std)
  s.exp.leaf <- cbind(mea[,c("DOY","Leaf")])
  colnames(s.exp.leaf) <- c("fractional_doy","Leaf") 
  r.exp.leaf <- reshape2::melt(s.exp.leaf, id.vars = "fractional_doy") 
  
  s.exp.std.leaf <- cbind(mea.std[,c("DOY","Leaf")])
  colnames(s.exp.std.leaf) <- c("fractional_doy","Leaf") 
  r.exp.std.leaf <- reshape2::melt(s.exp.std.leaf, id.vars = "fractional_doy") 
  r.exp.std.leaf$ymin <- r.exp.leaf$value - r.exp.std.leaf$value
  r.exp.std.leaf$ymax <- r.exp.leaf$value + r.exp.std.leaf$value
  
  # Stem
  s.exp.stem <- cbind(mea[,c("DOY","Stem")]) 
  colnames(s.exp.stem) <- c("fractional_doy","Stem") 
  r.exp.stem <- reshape2::melt(s.exp.stem, id.vars = "fractional_doy") 
  
  s.exp.std.stem <- cbind(mea.std[,c("DOY","Stem")]) 
  colnames(s.exp.std.stem) <- c("fractional_doy","Stem") 
  r.exp.std.stem <- reshape2::melt(s.exp.std.stem, id.vars = "fractional_doy") 
  r.exp.std.stem$ymin <- r.exp.stem$value - r.exp.std.stem$value
  r.exp.std.stem$ymax <- r.exp.stem$value + r.exp.std.stem$value
  
  # Pod
  s.exp.pod <- cbind(mea[,c("DOY","Pod")]) 
  colnames(s.exp.pod) <- c("fractional_doy","Pod") 
  r.exp.pod <- reshape2::melt(s.exp.pod, id.vars = "fractional_doy") 
  
  s.exp.std.pod <- cbind(mea.std[,c("DOY","Pod")]) 
  colnames(s.exp.std.pod) <- c("fractional_doy","Pod") 
  r.exp.std.pod <- reshape2::melt(s.exp.std.pod, id.vars = "fractional_doy") 
  r.exp.std.pod$ymin <- r.exp.pod$value - r.exp.std.pod$value
  r.exp.std.pod$ymax <- r.exp.pod$value + r.exp.std.pod$value
  
  # Combine
  r.exp.ls <- rbind(r.exp.leaf, r.exp.stem, r.exp.pod)
  r.exp.ls$Source = "Observed"
  r.lsrp.doy$Source = "Simulated"
  r.all <- rbind(r.lsrp.doy, r.exp.ls)
  # Reverse the order as follow
  r.all$Organ <- factor(r.all$variable, levels = rev(levels(r.all$variable)))
  # Colorblind friendly color palette (https://personal.sron.nl/~pault/)
  col.palette.muted.organs <- c("Leaf"="#117733", 
                                "Stem"= "#999933", 
                                "Root"="#332288", 
                                "Pod"= "#882255")
  
  f <- ggplot() + theme_classic()
  f <- f +
    geom_line(data = subset(r.all, Source == "Simulated"), na.rm = TRUE, 
              aes(x=fractional_doy,y=value, color=Organ), linewidth=1, alpha = 0.8) +
    geom_point(data = subset(r.all, Source == "Observed"), na.rm = TRUE,
               aes(x=fractional_doy, y=value, color=Organ), shape=15, size=2, stroke=.5) +
    
    scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, 2)) +
    scale_color_manual(values = col.palette.muted.organs)
  
  
  # for leaf
  f <- f + geom_errorbar(data=r.exp.std.leaf, aes(x=fractional_doy, ymin=ymin, ymax=ymax),  
                         width=3.5, size=0.25, show.legend = FALSE)
  # for stem
  f <- f + geom_errorbar(data=r.exp.std.stem, aes(x=fractional_doy, ymin=ymin, ymax=ymax),   
                         width=3.5, size=0.25, show.legend = FALSE)
  
  # for pod
  f <- f + geom_errorbar(data=r.exp.std.pod, aes(x=fractional_doy, ymin=ymin, ymax=ymax),   
                         width=3.5, size=0.25, show.legend = FALSE)
  
  f <- f + theme(plot.title=element_text(size=size.title, hjust=0.5),
                 axis.text=element_text(size=size.axis),
                 axis.title=element_text(size=size.axislabel),
                 axis.title.y = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 panel.background = element_rect(fill = "transparent",colour = NA))
  
  # change the plot labels and theme
  if(yr != '2003'){
    f <- f + labs(title= yr , x='Day of Year', y=NULL)
  }
  
  if(co2_opt == '_ambient_'){
    if (yr == '2002' || yr == '2005'){
      f <- f + theme(plot.background = element_rect(fill = "grey90", colour = NA))
    }else{
      f <- f + theme(plot.background = element_rect(fill = "transparent", colour = NA))
    }
  }
  f <- f + scale_x_continuous(breaks = seq(150,280,30))
  
  # calculate mse
  names(r.lsrp.doy) <- c('time', 'Organ', 'biomass', 'Source')
  names(r.exp.ls) <- c('time', 'Organ', 'biomass', 'Source')
  calculate_mse_soyFACE(yr,r.lsrp.doy, r.exp.ls)
  
  return(f)
}

plot_ld11_biomass <- function(result, mea){
  # Save plot into the list
  biocro_organ_biomass <- result[c('time', 'Leaf', 'Stem', 'Root', 'Pod')]
  biocro_organ_biomass_tall <- melt(biocro_organ_biomass, id.vars = 'time')
  names(biocro_organ_biomass_tall) <- c('time','Organ', 'biomass')
  
  field_organ_biomass <- mea[c('time', 'leaf', 'stem', 'root', 'pod')]
  names(field_organ_biomass)[names(field_organ_biomass) %in% c('leaf', 'stem', 'root', 'pod')] <- c('Leaf', 'Stem', 'Root', 'Pod')
  field_organ_biomass_tall <- melt(field_organ_biomass, id.vars = 'time')
  names(field_organ_biomass_tall) <- c('time','Organ', 'biomass')
  
  col.palette.muted <- c( "#117733", "#999933", "#332288", "#882255")
  
  f <- ggplot() + theme_classic() +
    geom_line(data = biocro_organ_biomass_tall, na.rm = TRUE, aes(x = time/24 + 1, y = biomass, color = Organ), linewidth = 1, alpha=0.8) +
    geom_point(data = field_organ_biomass_tall, na.rm = TRUE, aes(x = time/24 + 1, y = biomass, color = Organ), shape = 15, size = 2, stroke=.5)+
    theme(plot.title=element_text(size=size.title, hjust=0.5),
          axis.text=element_text(size=size.axis),
          axis.title.x =element_text(size=size.axislabel),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA))+
    scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, 2)) +
    scale_x_continuous(breaks = seq(150,280,30))+
    labs(title=years[i], 
         x=paste0('Day of Year'), 
         y='Biomass (Mg/ha)')+
    scale_color_manual(values = col.palette.muted)
  
  calculate_mse(years[[i]],biocro_organ_biomass_tall, field_organ_biomass_tall)
  
  return(f)
}