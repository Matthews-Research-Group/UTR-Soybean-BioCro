Data cant_day_2009_dailymeans;
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\Table_S3_canopy_temp\canopy_temp_2009_daily_ring_means.csv' 
dlm=',' firstobs=2 lrecl=20000; input block ring CO2$ DOY canopy_t;
run;

proc sort data =cant_day_2009_dailymeans;
by block co2 ring;
run ;

proc mixed data =cant_day_2009_dailymeans order =data covtest;
class block co2 doy ring ;
model canopy_t = co2|doy/ ddfm =kr;
random block ;
repeated doy / type =ar(1) subject =block*co2;
lsmeans co2|doy/ pdiff slice=doy ; 
ods output diffs = daycanopyT_pdiffs2009;
ods output lsmeans = daycanopyT_lsmeans2009;
*ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Canopy Temperature\2009daylightcanopyTpdiffs.csv';
proc print data = daycanopyT_pdiffs2009; where DOY=_DOY; run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Canopy Temperature\2009daylightcanopyTlsmeans.csv';
proc print data = daycanopyT_lsmeans2009;run;ODS CSV close;
