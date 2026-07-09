options ls=170 ps=10000 pageno=1;
title 'driFACE 2010 Water Potential ';
Data driFACE2010wp;
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\Table_S12_water_potential\dri2010WPringmeans.csv' 
dlm=',' firstobs=2 lrecl=10000; input DOY block ring CO2$ H2O$ WP OP TP @;
run;
**********Water potential ANOVA;
proc mixed data =driFACE2010wp order =data;
class DOY block CO2 H2O ;
model WP = DOY|CO2|H2O  / ddfm =kr;
random block block*co2 ;
lsmeans DOY|CO2|H2O/ pdiff; 
ods output lsmeans = WPseason;
ods output diffs = WPseasondiffs
run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Water Potential\Figures\2010_WP_lsmeans.csv';
proc print data = WPseason; run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Water Potential\Figures\2010_WP_pdiffs.csv';
proc print data = WPseasondiffs; where DOY=_DOY; run;ODS CSV close;

********Osmotic potential ANOVA;
proc mixed data =driFACE2010wp order =data;
class DOY block CO2 H2O ;
model OP = DOY|CO2|H2O  / ddfm =kr;
random block block*co2 ;
lsmeans DOY|CO2|H2O/ pdiff; 
ods output lsmeans = OPseason;
ods output diffs = OPseasondiffs
run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Water Potential\Figures\2010_OP_lsmeans.csv';
proc print data = OPseason; run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Water Potential\Figures\2010_OP_pdiffs.csv';
proc print data = OPseasondiffs; where DOY=_DOY; run;ODS CSV close;

********Turgor pressure ANOVA;
proc mixed data =driFACE2010wp order =data;
class DOY block CO2 H2O ;
model TP = DOY|CO2|H2O  / ddfm =kr;
random block block*co2 ;
lsmeans DOY|CO2|H2O/ pdiff; 
ods output lsmeans = TPseason;
ods output diffs = TPseasondiffs
run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Water Potential\Figures\2010_TP_lsmeans.csv';
proc print data = TPseason; run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Water Potential\Figures\2010_TP_pdiffs.csv';
proc print data = TPseasondiffs; where DOY=_DOY; run;ODS CSV close;


