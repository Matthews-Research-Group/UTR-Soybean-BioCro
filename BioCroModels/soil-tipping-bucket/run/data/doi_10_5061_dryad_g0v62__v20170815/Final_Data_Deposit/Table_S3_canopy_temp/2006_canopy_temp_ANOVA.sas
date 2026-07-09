Data cant_day_2006_dailymeans;
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\Table_S3_canopy_temp\canopy_temp_2006_daily_ring_means.csv'
dlm=',' firstobs=2 lrecl=20000; input block ring CO2$ DOY canopy_t;
run;

proc sort data =cant_day_2006_dailymeans;
by block co2 ring;
run ;

proc mixed data =cant_day_2006_dailymeans order =data covtest;
class block co2 doy ring ;
model canopy_t = co2|doy/ ddfm =kr;
random block ;
repeated doy / type =ar(1) subject =block*co2;
lsmeans co2|doy/ pdiff ; 
ods output diffs = daycanopyT_pdiffs2006;
ods output lsmeans = daycanopyT_lsmeans2006;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Canopy Temperature\2006daylightcanopyTpdiffs.csv';
proc print data = daycanopyT_pdiffs2006;run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Canopy Temperature\2006daylightcanopyTlsmeans.csv';
proc print data = daycanopyT_lsmeans2006;run;ODS CSV close;
