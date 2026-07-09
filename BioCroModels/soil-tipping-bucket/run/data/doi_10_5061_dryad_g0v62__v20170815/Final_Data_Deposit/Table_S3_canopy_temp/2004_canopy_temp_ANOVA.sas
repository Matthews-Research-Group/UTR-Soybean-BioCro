Data soyFACE2004canopyT;
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\Table_S3_canopy_temp\canopy_temp_2004_daily_ring_means.csv' 
dlm=',' firstobs=2 lrecl=20000; input block ring CO2$ DOY canopy_t;
run;

proc mixed data =soyFACE2004canopyT order =data covtest;
class block co2 doy ring ;
model canopy_t = co2|doy/ ddfm =kr;
random block ;
repeated doy / type =ar(1) subject =block*co2;
lsmeans co2|doy/ pdiff ; 
ods output diffs = middaycanopyT_pdiffs2004;
ods output lsmeans = middaycanopyT_lsmeans2004;
*ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Canopy Temperature\2004daylightcanopyTpdiffs.csv';
proc print data = middaycanopyT_pdiffs2004;run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Canopy Temperature\2004daylightcanopyTlsmeans.csv';
proc print data = middaycanopyT_lsmeans2004;run;ODS CSV close;
