Data dailyavg_gasex11;
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\ExtendedDataFig9_Table_S11_driFACE_gas_exchange\2011_daily_avg_gas_exchange.csv'
dlm=',' firstobs=2 lrecl=10000;
input Obs block DOY co2 h2o$ ring cica Photo cond ci @@;
run;
 
****************************the ANOVAS below are for daily average values*************************************************;
proc mixed data =dailyavg_gasex11 order =data covtest ;
class block co2 h2o DOY;
model ci = co2|h2o|DOY / ddfm =kr;
random block ;
lsmeans co2|h2o|DOY/ pdiff slice=DOY;
ods output diffs = dailyavg_ci11;
ods output lsmeans= dailyavg_cilsmeans_11;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Gas Exchange\SAS\Gas Exchange SAS Output\dailyavg_ci_lsmeans_11.csv';
proc print data =dailyavg_cilsmeans_11;  run;ODS CSV close;

ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Gas Exchange\SAS\Gas Exchange SAS Output\dailyavg_ci_pdiff_11.csv';
proc print data =dailyavg_ci11; where DOY= _DOY; run;ODS CSV close;

proc mixed data =dailyavg_gasex11 order =data covtest ;
class block co2 h2o DOY;
model cond = co2|h2o|DOY / ddfm =kr;
random block ;
lsmeans co2|h2o|DOY/ pdiff slice=DOY;
ods output diffs = dailyavg_cond11;
ods output lsmeans= dailyavg_condlsmeans_11;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Gas Exchange\SAS\Gas Exchange SAS Output\dailyavg_cond_lsmeans_11.csv';
proc print data =dailyavg_condlsmeans_11;  run;ODS CSV close;

ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Gas Exchange\SAS\Gas Exchange SAS Output\dailyavg_cond_pdiff_11.csv';
proc print data =dailyavg_cond11; where DOY= _DOY; run;ODS CSV close;

proc mixed data =dailyavg_gasex11 order =data covtest ;
class block co2 h2o DOY;
model photo = co2|h2o|DOY / ddfm =kr;
random block ;
lsmeans co2|h2o|DOY/ pdiff slice=DOY;
ods output diffs = dailyavg_photo11;
ods output lsmeans= dailyavg_photolsmeans_11;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Gas Exchange\SAS\Gas Exchange SAS Output\dailyavg_photo_lsmeans_11.csv';
proc print data =dailyavg_photolsmeans_11;  run;ODS CSV close;

ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Gas Exchange\SAS\Gas Exchange SAS Output\dailyavg_photo_pdiff_11.csv';
proc print data =dailyavg_photo11; where DOY= _DOY; run;ODS CSV close;

proc mixed data =dailyavg_gasex11 order =data covtest ;
class block co2 h2o DOY;
model cica = co2|h2o|DOY / ddfm =kr;
random block ;
lsmeans co2|h2o|DOY/ pdiff slice=DOY;
ods output diffs = dailyavg_cica11;
ods output lsmeans= dailyavg_cicalsmeans_11;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Gas Exchange\SAS\Gas Exchange SAS Output\dailyavg_cica_lsmeans_11.csv';
proc print data =dailyavg_cicalsmeans_11;  run;ODS CSV close;

ODS CSV file='C:\Users\sbgray\Documents\Summer 2011\Gas Exchange\SAS\Gas Exchange SAS Output\dailyavg_cica_pdiff_11.csv';
proc print data =dailyavg_cica11; where DOY= _DOY; run;ODS CSV close;
