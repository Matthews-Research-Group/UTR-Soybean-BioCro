# Plot partitioning
col.palette.muted <- c( "#117733", "#999933",  "#882255", "#332288")
plot_partitioning <- function(result, year){
  canopy_assim_daily <- aggregate(result$canopy_assimilation_rate,list(result$doy), 
                                  FUN=sum) * 0.6 / 180.156e-3
  leaf_reuse_daily <- aggregate(result$Leaf_senescence_rate * 
                                  parameters$Leaf_senescence_reuse_factor,
                                list(result$doy), FUN=sum)
  stem_reuse_daily <- aggregate(result$Stem_senescence_rate * 
                                  parameters$Stem_senescence_reuse_factor,
                                list(result$doy), FUN=sum)
  root_reuse_daily <- aggregate(result$Root_senescence_rate * 
                                  parameters$Root_senescence_reuse_factor,
                                list(result$doy), FUN=sum)
  avg_dvi_daily <- aggregate(result$DVI, list(result$doy), FUN=mean)
  
  net_subC_input <- data.frame(Group.1 = avg_dvi_daily$x, # canopy_assim_daily$Group.1
                               x = canopy_assim_daily$x + 
                                 leaf_reuse_daily$x + 
                                 stem_reuse_daily$x + 
                                 root_reuse_daily$x)
  leaf_export_daily <- aggregate(result$substrate_transport_Leaf_to_Stem,
                                 list(result$doy), FUN=sum)
  pod_allocation_daily <- aggregate(result$substrate_transport_Stem_to_Pod,
                                    list(result$doy), FUN=sum) 
  root_allocation_daily <- aggregate(result$substrate_transport_Stem_to_Root,
                                     list(result$doy), FUN=sum) + root_reuse_daily
  stem_allocation_daily <- leaf_export_daily + stem_reuse_daily - 
    pod_allocation_daily - root_allocation_daily
  leaf_allocation_daily <- canopy_assim_daily + leaf_reuse_daily - leaf_export_daily
  
  allocation_percentage <- data.frame(DVI = avg_dvi_daily$x, # DOY = leaf_export_daily$Group.1, 
                                      Leaf = 100 * leaf_allocation_daily$x / net_subC_input$x,
                                      Stem = 100 * stem_allocation_daily$x / net_subC_input$x,
                                      Root = 100 * root_allocation_daily$x / net_subC_input$x,
                                      Pod = 100 * pod_allocation_daily$x / net_subC_input$x)
  # To check if the sum is 100%
  # allocation_percentage$tot_percentage <- allocation_percentage$Leaf+
  #                                         allocation_percentage$Stem+
  #                                         allocation_percentage$Root+
  #                                         allocation_percentage$Pod
  allocation_percentage <- allocation_percentage[1:which.min(abs(allocation_percentage$DVI-parameters$stop_growth_dvi)),]
  
  allocation_percentage_tall <- melt(allocation_percentage, id.vars = 'DVI') # 'DOY')
  names(allocation_percentage_tall) <- c('DVI','Organ', 'Percentage') # DOY
  p <- ggplot() + theme_classic() +
    geom_point(data = allocation_percentage_tall, aes(x=DVI, y=Percentage, color=Organ), # x=DOY
               size = 0.8)+
    scale_color_manual(values = col.palette.muted)+
    theme(plot.title=element_text(size=16, hjust=0.5),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA))+
    scale_y_continuous(limits = c(-20, 120), breaks = seq(-20, 120, 20)) +
    scale_x_continuous(breaks = seq(0,2,0.5))+
    labs(title=element_blank(), 
         x= paste0('DVI (', year, ')'), # paste0('Day of Year (', year, ')'), 
         y='Allocation %')
  ggsave(paste0("allocation__percentage_", year, '.png'), width = 4, height = 3, units = "in")
  root_reuse_contribution_percentage <- root_reuse_daily$x/root_allocation_daily$x
  root_reuse_contribution_percentage <- root_reuse_contribution_percentage[
    -which(root_reuse_contribution_percentage>400)]
  
  # print(xyplot(root_reuse_daily$x/root_allocation_daily$x*100~
  #                root_reuse_daily$Group.1[
  #                  -which(root_reuse_contribution_percentage>400)],
  #              xlab='DOY', ylab='Percentage contribution from reuse'))
  
  # reuse percentage of all C source
  reuse_percentage <- data.frame(DVI = avg_dvi_daily$x, # DOY = leaf_reuse_daily$Group.1, 
                                 Remobilization_rate = 100*(leaf_reuse_daily$x+
                                                              stem_reuse_daily$x+
                                                              root_reuse_daily$x)/net_subC_input$x)
  reuse_percentage <- reuse_percentage[1:which.min(abs(reuse_percentage$DVI-parameters$stop_growth_dvi)),]
  reuse.p <- ggplot() + theme_classic() +
    geom_point(data = reuse_percentage, aes(x=DVI, y=Remobilization_rate), # DOY
               size = 0.8)+
    scale_color_manual(values = col.palette.muted)+
    theme(plot.title=element_text(size=16, hjust=0.5),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA))+
    scale_y_continuous(limits = c(-20, 120), breaks = seq(-20, 120, 20)) +
    scale_x_continuous(breaks = seq(0,2,0.5))+
    labs(title=element_blank(), 
         x=paste0('DVI (', year, ')'),# paste0('Day of Year (', year, ')'), 
         y='Remobolization %')
  ggsave(paste0("allocation__percentage_", year, '.png'), width = 4, height = 3, units = "in") # allocation__percentage_ or reuse__percentage_
  return(p) # possible returns: allocation_percentage_tall, p, reuse.p, depending on different purposes
}

