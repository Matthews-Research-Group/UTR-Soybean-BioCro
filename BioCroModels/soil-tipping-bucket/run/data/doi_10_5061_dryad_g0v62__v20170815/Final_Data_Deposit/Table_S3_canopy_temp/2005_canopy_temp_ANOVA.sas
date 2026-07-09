Data cant_2005_dailymeans;
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\Table_S3_canopy_temp\canopy_temp_2005_daily_ring_means.csv' 
dlm=',' firstobs=2 lrecl=20000; input block ring CO2$ DOY canopy_t;
run;

proc sort data =cant_2005_dailymeans;
by block co2 ring;
run ;
proc mixed data =cant_2005_dailymeans order =data covtest;
class block co2 doy ring ;
model canopy_t = co2|doy/ ddfm =kr;
random block ;
repeated doy / type =ar(1) subject =block*co2;
lsmeans co2|doy/ pdiff ; 
ods output diffs = daycanopyT_pdiffs2005;
ods output lsmeans = daycanopyT_lsmeans2005;
*ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Canopy Temperature\2005daylightcanopyTpdiffs.csv';
proc print data = daycanopyT_pdiffs2005; where DOY=_DOY; run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Canopy Temperature\2005daylightcanopyTlsmeans.csv';
proc print data = daycanopyT_lsmeans2005;run;ODS CSV close;
