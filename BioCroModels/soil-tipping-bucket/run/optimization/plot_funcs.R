# Define functions to create plots
plot_all_tissues <- function(res, year, biomass, biomass.std) {
  
  r <- reshape2::melt(res[, c("fractional_doy","Root","Leaf","Stem","Shell","Grain")], id.vars="fractional_doy")
  r.exp <- reshape2::melt(biomass[, c("DOY", "Leaf", "Stem", "Shell","Seed")], id.vars = "DOY")
  r.exp.std <- reshape2::melt(biomass.std[, c("DOY", "Leaf", "Stem", "Shell","Seed")], id.vars = "DOY")
  r.exp.std$ymin<-r.exp$value-r.exp.std$value
  r.exp.std$ymax<-r.exp$value+r.exp.std$value
  
  # Colorblind friendly color palette (https://personal.sron.nl/~pault/)
  col.palette.muted <- c("#332288", "#117733", "#999933", "#882255","#EE3377")
  
  size.title <- 12
  size.axislabel <-16
  size.axis <- 12
  size.legend <- 12
  
  f <- ggplot() + theme_classic()
  f <- f + geom_point(data=r, aes(x=fractional_doy,y=value, colour=variable), show.legend = TRUE, size=0.25)
  f <- f + geom_errorbar(data=r.exp.std, aes(x=DOY, ymin=ymin, ymax=ymax), width=3.5, size=.25, show.legend = FALSE)
  f <- f + geom_point(data=r.exp, aes(x=DOY, y=value, fill=variable), shape=22, size=2, show.legend = FALSE, stroke=.5)
  f <- f + labs(title=element_blank(), x=paste0('Day of Year (',year,')'),y='Biomass (Mg / ha)')
  f <- f + coord_cartesian(ylim = c(0,10)) + scale_y_continuous(breaks = seq(0,10,2)) + scale_x_continuous(breaks = seq(150,275,30))
  f <- f + theme(plot.title=element_text(size=size.title, hjust=0.5),
                 axis.text=element_text(size=size.axis),
                 axis.title=element_text(size=size.axislabel),
                 legend.position = c(.15,.85), legend.title = element_blank(),
                 legend.text=element_text(size=size.legend),
                 legend.background = element_rect(fill = "transparent",colour = NA),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
                 plot.background = element_rect(fill = "transparent", colour = NA))
  f <- f + guides(colour = guide_legend(override.aes = list(size=2)))
  f <- f + scale_fill_manual(values = col.palette.muted[2:5], guide = FALSE)
  f <- f + scale_colour_manual(values = col.palette.muted, labels=c('Root','Leaf','Stem','Shell','Seed'))
  
  return(f)
}

# Define functions to create plots
plot_all_tissues_TWO <- function(res1,res2, year, biomass, biomass.std) {
  
  r_control <- reshape2::melt(res1[, c("time","Root","Leaf","Stem","Shell","Grain")], id.vars="time")
  r_new     <- reshape2::melt(res2[, c("time","Root","Leaf","Stem","Shell","Grain")], id.vars="time")
  r.exp <- reshape2::melt(biomass[, c("DOY", "Leaf", "Stem", "Shell","Seed")], id.vars = "DOY")
  r.exp.std <- reshape2::melt(biomass.std[, c("DOY", "Leaf", "Stem", "Shell","Seed")], id.vars = "DOY")
  r.exp.std$ymin<-r.exp$value-r.exp.std$value
  r.exp.std$ymax<-r.exp$value+r.exp.std$value
  
  # Colorblind friendly color palette (https://personal.sron.nl/~pault/)
  col.palette.muted <- c("#332288", "#117733", "#999933", "#882255","#EE3377")
  
  size.title <- 12
  size.axislabel <-16
  size.axis <- 12
  size.legend <- 12
  
  f <- ggplot() + theme_classic()
  f <- f + geom_line(data=r_control, aes(x=time,y=value, colour=variable), show.legend = TRUE, size=1.5)
  f <- f + geom_line(data=r_new, aes(x=time,y=value, colour=variable), show.legend = FALSE, size=1.5,
                     linetype = "dotted")
  
  f <- f + geom_errorbar(data=r.exp.std, aes(x=DOY, ymin=ymin, ymax=ymax), width=3.5, size=.25, show.legend = FALSE)
  f <- f + geom_point(data=r.exp, aes(x=DOY, y=value, fill=variable), shape=22, size=2, show.legend = FALSE, stroke=.5)
  f <- f + labs(title=element_blank(), x=paste0('Day of Year (',year,')'),y='Biomass (Mg / ha)')
  f <- f + coord_cartesian(ylim = c(0,10)) + scale_y_continuous(breaks = seq(0,10,2)) + scale_x_continuous(breaks = seq(150,275,30))
  f <- f + theme(plot.title=element_text(size=size.title, hjust=0.5),
                 axis.text=element_text(size=size.axis),
                 axis.title=element_text(size=size.axislabel),
                 legend.position = c(.15,.85), legend.title = element_blank(),
                 legend.text=element_text(size=size.legend),
                 legend.background = element_rect(fill = "transparent",colour = NA),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
                 plot.background = element_rect(fill = "transparent", colour = NA))
  f <- f + guides(colour = guide_legend(override.aes = list(size=2)))
  f <- f + scale_fill_manual(values = col.palette.muted[2:5], guide = "none")
  f <- f + scale_colour_manual(values = col.palette.muted, labels=c('Root','Leaf','Stem','Shell','Seed'))
  
  return(f)
}

plot_litters <- function(res, year, biomass,biomass.std) {

  r <- reshape2::melt(res[, c("time","LeafLitter","StemLitter","TotalLitter")], id.vars="time")
  r.exp <- reshape2::melt(biomass[, c("DOY", "CumLitter")], id.vars = "DOY")
  r.exp.std <- reshape2::melt(biomass.std[, c("DOY", "CumLitter")], id.vars = "DOY")
  r.exp.std$ymin<-r.exp$value-r.exp.std$value
  r.exp.std$ymax<-r.exp$value+r.exp.std$value

  # Colorblind friendly color palette (https://personal.sron.nl/~pault/)
  col.palette.muted <- c("#332288", "#117733", "#999933", "#882255")

  size.title <- 12
  size.axislabel <-10
  size.axis <- 10
  size.legend <- 12

  f <- ggplot() + theme_classic()
  f <- f + geom_point(data=r, aes(x=time,y=value, colour=variable), show.legend = TRUE, size=0.25)
  f <- f + geom_point(data=r.exp, aes(x=DOY, y=value, fill=variable), shape=22, size=2, show.legend = FALSE, stroke=.5)
  f <- f + geom_errorbar(data=r.exp.std, aes(x=DOY, ymin=ymin, ymax=ymax), width=3.5, size=.25, show.legend = FALSE)
  f <- f + labs(title=element_blank(), x=paste0('Day of Year (',year,')'),y='Biomass (Mg / ha)')
  f <- f + coord_cartesian(ylim = c(0,5)) + scale_y_continuous(breaks = seq(0,5,1)) + scale_x_continuous(breaks = seq(150,275,30))
  f <- f + theme(plot.title=element_text(size=size.title, hjust=0.5),
                 axis.text=element_text(size=size.axis),
                 axis.title=element_text(size=size.axislabel),
                 legend.position = c(.15,.85), legend.title = element_blank(),
                 legend.text=element_text(size=size.legend),
                 legend.background = element_rect(fill = "transparent",colour = NA),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
                 plot.background = element_rect(fill = "transparent", colour = NA))
  f <- f + guides(colour = guide_legend(override.aes = list(size=2)))
  f <- f + scale_fill_manual(values = c("LeafLitter" = col.palette.muted[2],"StemLitter" = col.palette.muted[3],"TotalLitter"="red","CumLitter"="red"), guide = FALSE)
  f <- f + scale_colour_manual(values = c("LeafLitter" = col.palette.muted[2],"StemLitter" = col.palette.muted[3],"TotalLitter"="red"), labels=c('LeafLitter','StemLitter','TotalLitter'))

  return(f)
}

multiyear_biomass_plot <- function(predefine_eCO2,my_soybean_parameter,opt_result_name){
  using_Gray_precip  = TRUE
  obs_weather_Gray<-
      read.csv('../data/doi_10_5061_dryad_g0v62__v20170815/Final_Data_Deposit/ExtendedDataFig1_Weather_Data_2004thru2011/soyFACE_weather_data_2004thru2011.csv')
  
  new_par   = readRDS(paste0('opt_results/opt_result_',opt_result_name,'.rds'))
  arg_names = names(new_par)
  soybean_parameters0          = soybean2$parameters
  df=data.frame(arg_names,NEW=format(as.numeric(new_par),scientific=F),CTL=as.numeric(soybean_parameters0[arg_names]))
  print(df)
  
  my_soybean_parameter$parameters[arg_names] = new_par
  my_soybean_parameter$parameters['mrc_stem'] = my_soybean_parameter$parameters['mrc_leaf']
  
  years <- c('2002','2004','2005','2006')
  
  # sowing and harvest DOYs for each growing season
  dates <- data.frame("year" = c(2002, 2004:2006),"sow" = c(152,149,148,148), "harvest" = c(288, 289, 270, 270))
  
  # initialize variables
  results <- list()
  results_CTL <- list()
  weather_growing_season <- list()
  ExpBiomass <- list()
  ExpBiomass.std <- list()
  RootVals <- list()
  weights <- list()
  numrows <- vector()
  
  soybean_steadystate_modules0 = my_soybean_parameter$direct_modules 
  soybean_derivative_modules0  = my_soybean_parameter$differential_modules 
  soybean_initial_state0       = my_soybean_parameter$initial_values 
  soybean_parameters0          = my_soybean_parameter$parameters 
  print(soybean_parameters0[c('iSp','mrc_leaf','mrc_stem')])
  
  if(predefine_eCO2){
    #pre-define
    soybean_parameters0$RL_at_25    = soybean_parameters0$RL_at_25 * (1+0.37*0.42)
    soybean_parameters0$Vcmax_at_25 = soybean_parameters0$Vcmax_at_25*0.85
    soybean_parameters0$Jmax_at_25  = soybean_parameters0$Jmax_at_25*0.95
    soybean_parameters0$Catm  = 570 
  }
  
  for (i in 1:length(years)) {
    yr <- years[i]
    weather <- read.csv(file = paste0('Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))
  
    sowdate <- dates$sow[which(dates$year == yr)]
    harvestdate <- dates$harvest[which(dates$year == yr)]
    sd.ind <- which(weather$doy == sowdate)[1]
    hd.ind <- which(weather$doy == harvestdate)[24]
    
    weather_growing_season <- weather[sd.ind:hd.ind,]
    weather_growing_season$time_zone_offset = -6
  
    if(using_Gray_precip  & yr>=2004){ #Gray's does not have year 2002
      obs_weather_Gray_yeari = obs_weather_Gray[obs_weather_Gray$Year==yr,]
      obs_weather_Gray_yeari = obs_weather_Gray_yeari[obs_weather_Gray_yeari$DOY>=sowdate &obs_weather_Gray_yeari$DOY<=harvestdate,]
      obs_weather_Gray_yeari_hourly <- obs_weather_Gray_yeari %>%
        # Add an hourly sequence per day
        uncount(weights = 24, .id = "hour") %>%
        # Adjust the hour (0 to 23)
        mutate(hour = hour - 1,
               precip = precip.mm. / 24)
      
      # replace precip with Gray's
      # If Gray's data is shorter than the growing season, we only replace the shorter part
      weather_growing_season$precip[weather_growing_season$doy>=obs_weather_Gray_yeari$DOY[1] &weather_growing_season$doy<=tail(obs_weather_Gray_yeari$DOY,1)] = obs_weather_Gray_yeari_hourly$precip
    }
  
  #calculate CTL for diagnostic purposes
    results_CTL[[i]] <- run_biocro(
      soybean2$initial_values,
      soybean2$parameters,
      weather_growing_season,
      soybean2$direct_modules,
      soybean2$differential_modules,
      soybean2$ode_solver
    )
    
    results[[i]] <- run_biocro(
                              my_soybean_parameter$initial_values,
                              soybean_parameters0,
                              weather_growing_season,
                              my_soybean_parameter$direct_modules,
                              my_soybean_parameter$differential_modules,
                              my_soybean_parameter$ode_solver)
    
    print(paste("year",yr,".peak Leaf is,",max(results[[i]]$Leaf)))
    results[[i]]$time = results[[i]]$doy + results[[i]]$hour/24
    
    if(predefine_eCO2){
      ExpBiomass[[i]] <- read.csv(file=paste0('Data/biomasses_with_seed/',yr,'_co2_biomass.csv'))
    }else{
      ExpBiomass[[i]] <- read.csv(file=paste0('Data/biomasses_with_seed/',yr,'_ambient_biomass.csv'))
    }
    colnames(ExpBiomass[[i]])<-c("DOY","Leaf","Stem","Shell0","Seed","Litter","CumLitter")
    Shell = ExpBiomass[[i]]$Shell0 - ExpBiomass[[i]]$Seed
  #  Shell[which.max(Shell):length(Shell)] = max(Shell) #make Shell not decline
    ExpBiomass[[i]]$Shell = Shell 
  
    if(predefine_eCO2){
      ExpBiomass.std[[i]] <- read.csv(file=paste0('Data/biomasses_with_seed/',yr,'_co2_biomass_std.csv'))
    }else{
      ExpBiomass.std[[i]] <- read.csv(file=paste0('Data/biomasses_with_seed/',yr,'_ambient_biomass_std.csv'))
    }
    colnames(ExpBiomass.std[[i]])<-c("DOY","Leaf","Stem","Shell0","Seed","Litter","CumLitter")
    ExpBiomass.std[[i]]$Shell = sqrt((ExpBiomass.std[[i]]$Shell0)^2 + (ExpBiomass.std[[i]]$Seed)^2)
    
    RootVals[[i]] <- data.frame("DOY"=ExpBiomass[[i]]$DOY[5], "Root"=0.17*sum(ExpBiomass[[i]][5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
    
    numrows[i] <- nrow(weather_growing_season)
    invwts <- ExpBiomass.std[[i]]
    weights[[i]] <- log(1/(invwts[,2:ncol(invwts)]+1e-5))
  }
  
  #saveRDS(results,paste0("rds_results/result_",opt_result_name,".rds"))
  #stop()
  plot_list = list()
  plot_list_ctl = list()
  plot_list_litter = list()
  for (i in 1:length(years)){
  	yr = years[i]
  #	FigA <- plot_all_tissues_TWO(results_CTL[[i]],results[[i]], yr, ExpBiomass[[i]], ExpBiomass.std[[i]])
  	FigA <- plot_all_tissues(results[[i]], yr, ExpBiomass[[i]], ExpBiomass.std[[i]])
  	plot_list[[i]] = FigA
  
  	FigAA <- plot_all_tissues(results_CTL[[i]], yr, ExpBiomass[[i]], ExpBiomass.std[[i]])
  	plot_list_ctl[[i]] = FigAA
          
          results[[i]]$TotalLitter = results[[i]]$LeafLitter + results[[i]]$StemLitter
  
  	FigB <- plot_litters(results[[i]],yr, ExpBiomass[[i]], ExpBiomass.std[[i]])
  	plot_list_litter[[i]] = FigB
  }
  
  if(predefine_eCO2){
    pdf_name = paste0("figs/fig_",opt_result_name,"_eCO2.pdf")
  }else{
    pdf_name = paste0("figs/fig_",opt_result_name,"_aCO2.pdf")
   # pdf_name2 = paste0("figs/fig_",opt_result_name,"_aCO2_ctl.pdf")
   # pdf(pdf_name2,height = 8, width=8,bg='transparent')
   # grid.arrange(grobs = plot_list_ctl,nrow=2,ncol=2)
   # dev.off()
  }
  pdf(pdf_name,height = 8, width=8,bg='transparent')
  grid.arrange(grobs = plot_list,nrow=2,ncol=2)
  dev.off()
  
  #litter biomass plot
  if(predefine_eCO2){
    pdf_name = paste0("figs/fig_",opt_result_name,"_litter_eCO2.pdf")
  }else{
    pdf_name = paste0("figs/fig_",opt_result_name,"_litter_aCO2.pdf")
  }
  pdf(pdf_name,height = 8, width=8,bg='transparent')
  grid.arrange(grobs = plot_list_litter,nrow=2,ncol=2)
  dev.off()
}
