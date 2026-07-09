options ls=170 ps=10000 pageno=1;
title 'Soil VWC environment variables';

Data soyFACEVWCenvironment;
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\Fig_3_Fig_S5_soil_moisture_regressions\VWC_environment_data.csv' 
dlm=',' firstobs=2 lrecl=20000; 
input year$ airTpeakLAI percLAI precippeakLAI peakVWC75 ccVWC75 dailymaxdeltacanT gsclosedcan;
run;

*Fig. 3A: Regression of % effect of elevated CO2 on soil moisture vs. % effect of elevated CO2 on maximum LAI;
*2011 LAI data were unreliable due to very large variation between instruments, so LAI data for this year were not used;
data no2011;
set soyFACEVWCenvironment;
if year=2011 then delete;run;
proc reg data=no2011;
model ccVWC75=percLAI;
run;

*Fig 3B: Regression of % effect of elevated CO2 on soil moisture vs. seasonal average increase in daily maximum canopy surface
temperature induced by elevated CO2 treatment during the period of canopy closure (LAI>3);
*Canopy temperature data for 2008 were unreliable due to problems with sensor calibration, so 2008 data are excluded from analysis;
proc reg data=soyFACEVWCenvironment;
model ccVWC75=dailymaxdeltacanT;
run;

*Fig. 3C:Regression of % effect of elevated CO2 on soil moisture vs. precipitation (mm) and average air temperature (deg C)
during the period of canopy development (planting through peak LAI);
proc reg data=soyFACEVWCenvironment;
model peakVWC75=precippeakLAI airTpeakLAI /stb vif tol collin stb;
run;



