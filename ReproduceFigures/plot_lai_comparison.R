# Colorblind friendly color palette (https://personal.sron.nl/~pault/)
col.palette.muted <- c("#332288", "#117733", "#999933", "#882255")
size.title <- 12
size.axislabel <-10
size.axis <- 10
size.legend <- 12
plot_amb_elev_lai <- function(res, elev_res, year, lai, elev_lai) {
  
  
  res$lai_tot <- (res$Leaf_substrate_carbon+res$Leaf_structural_carbon)* 0.3 * parameters$iSp
  elev_res$lai_tot <- (elev_res$Leaf_substrate_carbon+elev_res$Leaf_structural_carbon)* 0.3 * parameters$iSp
  
  s.lai <- cbind(res[,c("fractional_doy","lai_tot")],elev_res[,"lai_tot"])
  colnames(s.lai) <- c("fractional_doy","Amb","Elev")
  r.lai <- reshape2::melt(s.lai, id.vars = "fractional_doy")
  
  s.exp.lai <- cbind(lai[,c("DOY","LAI_mean")],elev_lai[,"LAI_mean"])
  colnames(s.exp.lai) <- c("DOY","Amb","Elev")
  r.exp.lai <- reshape2::melt(s.exp.lai, id.vars = "DOY")
  
  s.exp.std.lai <- cbind(lai[,c("DOY","LAI_std")],elev_lai[,"LAI_std"])
  colnames(s.exp.std.lai) <- c("DOY","Amb","Elev")
  r.exp.std.lai <- reshape2::melt(s.exp.std.lai, id.vars = "DOY")
  r.exp.std.lai$ymin <- r.exp.lai$value - r.exp.std.lai$value
  r.exp.std.lai$ymax <- r.exp.lai$value + r.exp.std.lai$value
  
  f <- ggplot() + theme_classic()
  f <- f + geom_point(data=r.lai, aes(x=fractional_doy, y=value, colour=variable),show.legend = FALSE,size=0.25)
  f <- f + geom_errorbar(data=r.exp.std.lai, aes(x=DOY, ymin=ymin, ymax=ymax), width=3.5, size=0.25, show.legend = FALSE)
  f <- f + geom_point(data=r.exp.lai, aes(x=DOY, y=value, fill=variable), shape=22, size=2, show.legend = FALSE, stroke=.5)
  f <- f + coord_cartesian(ylim = c(0,10)) + scale_x_continuous(breaks = seq(150,275,30))
  f <- f + labs(x='Day of Year' ,y=bquote("LAI"~(m^2~"/"~m^2)))
  f <- f + theme(plot.title=element_text(size=size.title, hjust=0.5),
                 axis.text=element_text(size=size.axis),
                 axis.title=element_text(size=size.axislabel),
                 axis.title.y = element_blank(),
                 legend.position = c(.25,.85), legend.title = element_blank(),
                 legend.text=element_text(size=size.legend),
                 legend.background = element_rect(fill = "transparent",colour = NA),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
                 plot.background = element_rect(fill = "transparent", colour = NA))
  f <- f + guides(colour = guide_legend(override.aes = list(size=2)))
  f <- f + scale_fill_manual(values = col.palette.muted[2:3], guide = "none")
  f <- f + scale_colour_manual(values = col.palette.muted[2:3],labels=c('Ambient',bquote(Elevated~CO[2])))
  
  return(f)
}