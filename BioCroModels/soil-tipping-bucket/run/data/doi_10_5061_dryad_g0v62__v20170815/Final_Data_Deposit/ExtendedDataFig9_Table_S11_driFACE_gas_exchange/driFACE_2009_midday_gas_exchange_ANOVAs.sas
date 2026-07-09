Data dailyavg_gasex09;
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\ExtendedDataFig9_Table_S11_driFACE_gas_exchange\2009_daily_avg_gas_exchange.csv'
dlm=',' firstobs=2 lrecl=10000;
input Obs block DOY co2 h2o$ ring cica Photo cond ci @@;
run;


******the ANOVAs below are for whole season data, using daily averages***************************************;

proc mixed data =dailyavg_gasex09 order =data covtest ;
class block co2 h2o DOY;
model cica = co2|h2o|DOY / ddfm =kr;
random block ;
lsmeans co2|DOY|h2o/ pdiff slice=DOY;
ods output diffs = dailyavg_cica;
ods output lsmeans= dailyavg_lsmeans;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2009\Gas Exchange\SAS\dailyavgcicalsmeans.csv';
proc print data =dailyavg_lsmeans;  run;ODS CSV close;

ODS CSV file='C:\Users\sbgray\Documents\Summer 2009\Gas Exchange\SAS\dailyavg_cica_pdiff.csv';
proc print data =dailyavg_cica; where DOY= _DOY; run;ODS CSV close;

proc mixed data =dailyavg_gasex09 order =data covtest ;
class block co2 h2o DOY;
model cond = co2|h2o|DOY / ddfm =kr;
random block ;
lsmeans co2|DOY|h2o/ pdiff slice=DOY;
ods output diffs = dailyavg_cond;
ods output lsmeans= dailyavg_lsmeans;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2009\Gas Exchange\SAS\dailyavgcondlsmeans.csv';
proc print data =dailyavg_lsmeans;  run;ODS CSV close;

ODS CSV file='C:\Users\sbgray\Documents\Summer 2009\Gas Exchange\SAS\dailyavg_cond_pdiff.csv';
proc print data =dailyavg_cond; where DOY= _DOY; run;ODS CSV close;

proc mixed data =dailyavg_gasex09 order =data covtest ;
class block co2 h2o DOY;
model ci = co2|h2o|DOY / ddfm =kr;
random block ;
lsmeans co2|DOY|h2o/ pdiff slice=DOY;
ods output diffs = dailyavg_ci;
ods output lsmeans= dailyavg_lsmeans;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2009\Gas Exchange\SAS\dailyavgcilsmeans.csv';
proc print data =dailyavg_lsmeans;  run;ODS CSV close;

ODS CSV file='C:\Users\sbgray\Documents\Summer 2009\Gas Exchange\SAS\dailyavg_ci_pdiff.csv';
proc print data =dailyavg_ci; where DOY= _DOY; run;ODS CSV close;

proc mixed data =dailyavg_gasex09 order =data covtest ;
class block co2 h2o DOY;
model Photo = co2|h2o|DOY / ddfm =kr;
random block ;
lsmeans co2|DOY|h2o/ pdiff slice=DOY;
ods output diffs = dailyavg_Photo;
ods output lsmeans= dailyavg_lsmeans;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2009\Gas Exchange\SAS\dailyavgPhotolsmeans.csv';
proc print data =dailyavg_lsmeans;  run;ODS CSV close;

ODS CSV file='C:\Users\sbgray\Documents\Summer 2009\Gas Exchange\SAS\dailyavg_Photo_pdiff.csv';
proc print data =dailyavg_Photo; where DOY= _DOY; run;ODS CSV close;
