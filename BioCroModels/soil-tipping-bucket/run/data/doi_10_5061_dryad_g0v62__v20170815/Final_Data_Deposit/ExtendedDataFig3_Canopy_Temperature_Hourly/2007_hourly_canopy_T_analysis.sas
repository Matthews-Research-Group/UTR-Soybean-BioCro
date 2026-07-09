Data soyFACE2007canopyT;
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\Fig_S3_Canopy_Temperature_Hourly\hourly_canopy_T_2007.csv'
dlm=',' firstobs=2 lrecl=20000; input block ring CO2$ hour canopy_t;
run;

proc mixed data =soyFACE2007canopyT order =data covtest;
class block co2 hour ring ;
model canopy_t = co2|hour/ ddfm =kr;
random block ;
repeated hour / type =ar(1) subject =block*co2;
lsmeans co2|hour/ pdiff slice=hour ; 
ods output diffs = hourlycanT_pdiffs2007;
ods output lsmeans = hourlycanT_lsmeans2007;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Canopy Temperature\canopyTlsmeans_hourly2007.csv';
proc print data = hourlycanT_lsmeans2007;run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Canopy Temperature\canopyTpdiffs_hourly2007.csv';
proc print data = hourlycanT_pdiffs2007; where hour=_hour; run;ODS CSV close;
