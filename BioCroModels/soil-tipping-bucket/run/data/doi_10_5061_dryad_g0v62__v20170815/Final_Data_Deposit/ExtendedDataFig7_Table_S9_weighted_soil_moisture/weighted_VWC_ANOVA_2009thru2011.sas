*Infile 2009 data;
Data DriFACE2009_all;
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\ExtendedDataFig7_Table_S8_weighted_soil_moisture\VWC_RLD_2009.csv' 
dlm=',' firstobs=2 lrecl=10000; input year DOY block ring CO2 H2O$ VWCRLD ;;

*ANOVA for 2009 soil moisture data, weighted for root distribution;
proc mixed data =driFACE2009_all order =data covtest;
class block DOY co2 h2o ;
model VWCRLD = co2|h2o|DOY/ ddfm =kr;
random block ;
repeated DOY / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|DOY / pdiff slice=DOY; 
ods output diffs = VWC_RLD_pdiffs09;
ods output lsmeans = VWC_RLD_lsmeans09;
run ;
ODS CSV file='C:\Users\sbgray\Documents\driFACE 2009-2011 Manuscript\SAS library\VWC_RLD_lsmeans09.csv';
proc print data =VWC_RLD_lsmeans09; run ; ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\driFACE 2009-2011 Manuscript\SAS library\VWC_RLD_pdiffs09.csv';
proc print data =VWC_RLD_pdiffs09; run ; ODS CSV close;

*Infile 2010 data;
Data DriFACE2010_all;
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\ExtendedDataFig7_Table_S8_weighted_soil_moisture\VWC_RLD_2010.csv' 
dlm=',' firstobs=2 lrecl=10000; input year DOY block ring CO2 H2O$ VWCRLD ;;

*ANOVA for 2010 soil moisture data, weighted for root distribution;
proc mixed data =driFACE2010_all order =data covtest;
class block DOY co2 h2o ;
model VWCRLD = co2|h2o|DOY/ ddfm =kr;
random block ;
repeated DOY / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|DOY / pdiff slice=DOY; 
ods output diffs = VWC_RLD_pdiffs10;
ods output lsmeans = VWC_RLD_lsmeans10;
run ;
ODS CSV file='C:\Users\sbgray\Documents\driFACE 2009-2011 Manuscript\SAS library\VWC_RLD_lsmeans10.csv';
proc print data =VWC_RLD_lsmeans10; run ; ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\driFACE 2009-2011 Manuscript\SAS library\VWC_RLD_pdiffs10.csv';
proc print data =VWC_RLD_pdiffs10; run ; ODS CSV close;

*Infile 2011 data;
Data DriFACE2011_all;
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\ExtendedDataFig7_Table_S8_weighted_soil_moisture\VWC_RLD_2011.csv' 
dlm=',' firstobs=2 lrecl=10000; input year DOY block ring CO2 H2O$ VWCRLD ;;

*ANOVA for 2011 soil moisture data, weighted for root distribution;
proc mixed data =driFACE2011_all order =data covtest;
class block DOY co2 h2o ;
model VWCRLD = co2|h2o|DOY/ ddfm =kr;
random block ;
repeated DOY / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|DOY / pdiff slice=DOY; 
ods output diffs = VWC_RLD_pdiffs11;
ods output lsmeans = VWC_RLD_lsmeans11;
run ;

ODS CSV file='C:\Users\sbgray\Documents\driFACE 2009-2011 Manuscript\SAS library\VWC_RLD_lsmeans11.csv';
proc print data =VWC_RLD_lsmeans11; run ; ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\driFACE 2009-2011 Manuscript\SAS library\VWC_RLD_pdiffs11.csv';
proc print data =VWC_RLD_pdiffs11; run ; ODS CSV close;
