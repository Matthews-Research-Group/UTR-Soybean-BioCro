**the following code infiles the monthly weather data csv file from the Illinois State Water Survey website***;
data month_avg;
infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\ExtendedDataFig5_Long_Term_Weather_Data\ISWS_Urbana_monthly_weather_1901thru2008.csv'
dlm=',' firstobs=2 lrecl=20000;
input Year$ Month$ Precip Mean_T;  
run;

**remove summary data lines and other month averages from the month_avg dataset;
data summer_month_avg;
set 
month_avg;
if Month = 'Total' then delete;
if Month = 'Winter' then delete;
if Month = 'Spring' then delete;
if Month = 'Summer' then delete;
if Month = 'Fall' then delete;
if Month = '1' then delete;
if Month = '2' then delete;
if Month = '3' then delete;
if Month = '4' then delete;
if Month = '5' then delete;
if Month = '10' then delete;
if Month = '11' then delete;
if Month = '12' then delete;
run;

**convert degrees F into degrees C and precipitation (inches) into precipitation (mm);
data summer_month_avg;
set 
summer_month_avg;
Temp_C=(Mean_T - 32)/1.8;
precip_mm = Precip*25.4;
run;

*calculate mean jun-sep temperature and total jun-sep precipitation for each year;
proc means mean noprint data=summer_month_avg;
class year ; var Temp_C ;
output mean = out = mean_temp;
run;
ODS CSV file='C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\Fig_S6_Long_Term_Weather_Data\Urbana_JunthruSep_temp_1901thru2008.csv';
proc print data = mean_temp; where _type_=1; run; ODS CSV close;

*calculate total jun-sep precipitation for each year;
proc means sum noprint data=summer_month_avg;
class year ; var precip_mm;
output sum = out = total_precip;
run;

ODS CSV file='C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\ExtendedDataFig5_Long_Term_Weather_Data\Urbana_JunthruSep_precip_1901thru2008.csv';
proc print data = total_precip; where _type_=1; run; ODS CSV close;


