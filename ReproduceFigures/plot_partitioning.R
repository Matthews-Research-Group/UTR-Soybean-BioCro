# Plot partitioning
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
  net_subC_input <- data.frame(Group.1 = canopy_assim_daily$Group.1,
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
  
  allocation_percentage <- data.frame(DOY = leaf_export_daily$Group.1, 
                                      Leaf = 100 * leaf_allocation_daily$x / net_subC_input$x,
                                      Stem = 100 * stem_allocation_daily$x / net_subC_input$x,
                                      Root = 100 * root_allocation_daily$x / net_subC_input$x,
                                      Pod = 100 * pod_allocation_daily$x / net_subC_input$x)
  # To check if the sum is 100%
  # allocation_percentage$tot_percentage <- allocation_percentage$Leaf+
  #                                         allocation_percentage$Stem+
  #                                         allocation_percentage$Root+
  #                                         allocation_percentage$Pod
  allocation_percentage <- allocation_percentage[1:which(allocation_percentage$DOY==result$doy
                                                         [which.min(abs(result$DVI-parameters$stop_growth_dvi))])-1,]
  
  allocation_percentage_tall <- melt(allocation_percentage, id.vars = 'DOY')
  names(allocation_percentage_tall) <- c('DOY','Organ', 'Percentage')
  p <- ggplot() + theme_classic() +
    geom_point(data = allocation_percentage_tall, aes(x=DOY, y=Percentage, color=Organ),
               size = 0.8)+
    theme(plot.title=element_text(size=16, hjust=0.5),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA))+
    scale_y_continuous(limits = c(-20, 120), breaks = seq(-20, 120, 20)) +
    scale_x_continuous(breaks = seq(180,280,30))+
    labs(title=element_blank(), 
         x=paste0('Day of Year (', year, ')'), 
         y='Allocation %')
  print(p)
  root_reuse_contribution_percentage <- root_reuse_daily$x/root_allocation_daily$x
  root_reuse_contribution_percentage <- root_reuse_contribution_percentage[
    -which(root_reuse_contribution_percentage>400)]
  print(xyplot(root_reuse_daily$x/root_allocation_daily$x*100~
                 root_reuse_daily$Group.1[
                   -which(root_reuse_contribution_percentage>400)]),
        xlab='DOY',
        ylab='Percentage contribution from reuse')
  
  return(allocation_percentage_tall)
}