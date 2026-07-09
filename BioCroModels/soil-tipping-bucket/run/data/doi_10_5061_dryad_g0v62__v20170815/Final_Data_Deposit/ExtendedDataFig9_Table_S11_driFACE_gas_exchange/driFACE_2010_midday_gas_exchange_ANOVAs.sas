Data dailyavg_gasex10;
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\ExtendedDataFig9_Table_S11_driFACE_gas_exchange\2010_daily_avg_gas_exchange.csv'
dlm=',' firstobs=2 lrecl=10000;
input Obs block DOY co2 h2o$ ring cica Photo cond ci @@;
run;

/***these anovas use the whole season's data in the form of daily averages*************************************************;*/;
proc mixed data =dailyavg_gasex10 order =data covtest ;
class block co2 h2o DOY;
model ci = co2|h2o|DOY / ddfm =kr;
random block ;
lsmeans co2|h2o|DOY/ pdiff slice=DOY;
ods output diffs = dailyavg_ci_pdiffs;
ods output lsmeans= dailyavg_ci_lsmeans;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Gas Exchange\SAS\Gas Exchange SAS Output\daily_avg_ci_lsmeans.csv';
proc print data =dailyavg_ci_lsmeans;  run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Gas Exchange\SAS\Gas Exchange SAS Output\daily_avg_ci_pdiff.csv';
proc print data =dailyavg_ci_pdiffs; where DOY= _DOY; run;ODS CSV close;

proc mixed data =dailyavg_gasex10 order =data covtest ;
class block co2 h2o DOY;
model cond = co2|h2o|DOY / ddfm =kr;
random block ;
lsmeans co2|h2o|DOY/ pdiff slice=DOY;
ods output diffs = dailyavg_cond_pdiffs;
ods output lsmeans= dailyavg_cond_lsmeans;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Gas Exchange\SAS\Gas Exchange SAS Output\daily_avg_cond_lsmeans.csv';
proc print data =dailyavg_cond_lsmeans;  run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Gas Exchange\SAS\Gas Exchange SAS Output\daily_avg_cond_pdiff.csv';
proc print data =dailyavg_cond_pdiffs; where DOY= _DOY; run;ODS CSV close;


proc mixed data =dailyavg_gasex10 order =data covtest ;
class block co2 h2o DOY;
model Photo = co2|h2o|DOY / ddfm =kr;
random block ;
lsmeans co2|h2o|DOY/ pdiff slice=DOY;
ods output diffs = dailyavg_Photo_pdiffs;
ods output lsmeans= dailyavg_Photo_lsmeans;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Gas Exchange\SAS\Gas Exchange SAS Output\daily_avg_Photo_lsmeans.csv';
proc print data =dailyavg_Photo_lsmeans;  run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Gas Exchange\SAS\Gas Exchange SAS Output\daily_avg_Photo_pdiff.csv';
proc print data =dailyavg_Photo_pdiffs; where DOY= _DOY; run;ODS CSV close;

proc mixed data =dailyavg_gasex10 order =data covtest ;
class block co2 h2o DOY;
model cica = co2|h2o|DOY / ddfm =kr;
random block ;
lsmeans co2|h2o|DOY/ pdiff slice=DOY;
ods output diffs = dailyavg_cica_pdiffs;
ods output lsmeans= dailyavg_cica_lsmeans;
ods listing exclude diffs;
run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Gas Exchange\SAS\Gas Exchange SAS Output\daily_avg_cica_lsmeans.csv';
proc print data =dailyavg_cica_lsmeans;  run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Gas Exchange\SAS\Gas Exchange SAS Output\daily_avg_cica_pdiff.csv';
proc print data =dailyavg_cica_pdiffs; where DOY= _DOY; run;ODS CSV close;
