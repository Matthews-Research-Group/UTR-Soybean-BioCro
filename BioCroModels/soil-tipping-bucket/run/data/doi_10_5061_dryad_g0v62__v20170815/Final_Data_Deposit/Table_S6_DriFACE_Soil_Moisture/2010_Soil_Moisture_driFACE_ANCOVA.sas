Data soystats2010_driFACE;
Infile 'C:\Users\sbgray\Documents\Combined Soil H2O and driFACE manuscript\Raw data and code\Final\Table_S6_DriFACE_Soil_Moisture\2010_driFACE_soil_moisture_ringmeans.csv'
dlm=',' firstobs=2 lrecl=10000; input Obs ring doy h2o$ co2$ block 
true_ten true_twenty true_thirty true_forty true_fifty true_sixty true_seventy 
start_ten start_twenty start_thirty start_forty start_fifty start_sixty start_seventy @@;
run;
proc sort data =soystats2010_driFACE;
by block doy co2 h2o;
run ;
proc mixed data =soystats2010_driFACE order =data covtest;
class block doy co2 h2o ;
model true_ten = co2|h2o|doy start_ten/ ddfm =kr;
random block block*co2 ;
repeated doy / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|doy/ pdiff slice=doy; 
ods output diffs = co2_h2oeffect;
ods output lsmeans = co2_h2o_doylsmeans;
run ;

proc print data =co2_h2o_doylsmeans; run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Soil Moisture\Final SAS Library\ten_lsmeans.csv';
proc print data =co2_h2o_doylsmeans; run; ODS CSV close;

proc print data =co2_h2oeffect; where doy= _doy; run;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Soil Moisture\Final SAS Library\ten_pdiffs.csv';
proc print data =co2_h2oeffect; where doy= _doy; run; ODS CSV close;

proc mixed data =soystats2010_driFACE order =data covtest;
class block doy co2 h2o ;
model true_twenty = co2|h2o|doy start_twenty/ ddfm =kr;
random block block*co2 ;
repeated doy / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|doy/ pdiff slice=doy; 
ods output diffs = co2_h2oeffect;
ods output lsmeans = co2_h2o_doylsmeans;
ods listing exclude diffs;
run ;

proc print data =co2_h2o_doylsmeans; run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Soil Moisture\Final SAS Library\twenty_lsmeans.csv';
proc print data =co2_h2o_doylsmeans; run; ODS CSV close;

proc print data =co2_h2oeffect; where doy= _doy; run;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Soil Moisture\Final SAS Library\twenty_pdiffs.csv';
proc print data =co2_h2oeffect; where doy= _doy; run; ODS CSV close;


proc mixed data =soystats2010_driFACE order =data covtest;
class block doy co2 h2o;
model true_thirty = co2|h2o|doy start_thirty/ ddfm =kr;
random block block*co2 ;
repeated doy / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|doy/ pdiff slice=doy; 
ods output diffs = co2_h2oeffect;
ods output lsmeans = co2_h2o_doylsmeans;
ods listing exclude diffs;
run ;

proc print data =co2_h2o_doylsmeans; run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Soil Moisture\Final SAS Library\thirty_lsmeans.csv';
proc print data =co2_h2o_doylsmeans; run; ODS CSV close;

proc print data =co2_h2oeffect; where doy= _doy; run;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Soil Moisture\Final SAS Library\thirty_pdiffs.csv';
proc print data =co2_h2oeffect; where doy= _doy; run; ODS CSV close;


proc mixed data =soystats2010_driFACE order =data covtest;
class block doy co2 h2o;
model true_forty = co2|h2o|doy start_forty/ ddfm =kr;
random block block*co2 ;
repeated doy / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|doy/ pdiff slice=doy; 
ods output diffs = co2_h2oeffect;
ods output lsmeans = co2_h2o_doylsmeans;
ods listing exclude diffs;
run ;

proc print data =co2_h2o_doylsmeans; run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Soil Moisture\Final SAS Library\forty_lsmeans.csv';
proc print data =co2_h2o_doylsmeans; run; ODS CSV close;

proc print data =co2_h2oeffect; where doy= _doy; run;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Soil Moisture\Final SAS Library\forty_pdiffs.csv';
proc print data =co2_h2oeffect; where doy= _doy; run; ODS CSV close;

proc mixed data =soystats2010_driFACE order =data covtest;
class block doy co2 h2o;
model true_fifty = co2|h2o|doy start_fifty/ ddfm =kr;
random block block*co2 ;
repeated doy / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|doy/ pdiff slice=doy; 
ods output diffs = co2_h2oeffect;
ods output lsmeans = co2_h2o_doylsmeans;
run ;

proc print data =co2_h2o_doylsmeans; run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Soil Moisture\Final SAS Library\fifty_lsmeans.csv';
proc print data =co2_h2o_doylsmeans; run; ODS CSV close;

proc print data =co2_h2oeffect; where doy= _doy; run;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Soil Moisture\Final SAS Library\fifty_pdiffs.csv';
proc print data =co2_h2oeffect; where doy= _doy; run; ODS CSV close;

proc mixed data =soystats2010_driFACE order =data covtest;
class block doy co2 h2o;
model true_sixty = co2|h2o|doy start_sixty/ ddfm =kr;
random block block*co2 ;
repeated doy / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|doy/ pdiff slice=doy; 
ods output diffs = co2_h2oeffect;
ods output lsmeans = co2_h2o_doylsmeans;
run ;

proc print data =co2_h2o_doylsmeans; run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Soil Moisture\Final SAS Library\sixty_lsmeans.csv';
proc print data =co2_h2o_doylsmeans; run; ODS CSV close;

proc print data =co2_h2oeffect; where doy= _doy; run;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Soil Moisture\Final SAS Library\sixty_pdiffs.csv';
proc print data =co2_h2oeffect; where doy= _doy; run; ODS CSV close;

proc mixed data =soystats2010_driFACE order =data covtest;
class block doy co2 h2o;
model true_seventy = co2|h2o|doy start_seventy/ ddfm =kr;
random block block*co2 ;
repeated doy / type =ar(1) subject =block*co2*h2o;
lsmeans co2|h2o|doy/ pdiff slice=doy; 
ods output diffs = co2_h2oeffect;
ods output lsmeans = co2_h2o_doylsmeans;
ods listing exclude diffs;
run ;

proc print data =co2_h2o_doylsmeans; run ;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Soil Moisture\Final SAS Library\seventy_lsmeans.csv';
proc print data =co2_h2o_doylsmeans; run; ODS CSV close;

proc print data =co2_h2oeffect; where doy= _doy; run;
ODS CSV file='C:\Users\sbgray\Documents\Summer 2010\Soil Moisture\Final SAS Library\seventy_pdiffs.csv';
proc print data =co2_h2oeffect; where doy= _doy; run; ODS CSV close;
