CheckLeapYear<-function(year)
  {
    if( (year%%400==0) || ((year%%100 !=0) && (year%%4==0))) temp<-1 else temp<-0
    return (temp)
  }
years = c(2002,2004:2006)
td = -5
for (yr in years){
#	fn = paste0(yr,"_Bondville_IL_daylength.csv")
	fn = paste0(yr,"_nasa_climate.csv")
	x = read.csv(fn)
	hour = x$hour
	hr <- hour-td
	hr <- hr%%24
	x$hour = hr
	x1 = x[c(1,1,1,1,1,1:(dim(x)[1]+td)),]	
	x1$hour = c(0:4,x1$hour[-(1:5)])
	doy = rep(1:365,each=24)
	x1$doy = doy
#	if(CheckLeapYear(yr)==1){
#		tmp = x1[x1$doy==59,]
#		x1 = rbind(x1[x1$doy<=59,],tmp,x1[x1$doy>59,])	
#		doy = rep(1:366,each=24)
#		x1$doy = doy
#	}
	write.csv(x1,paste0(yr,"_nasa_GMT.csv"),row.names = FALSE)
}
