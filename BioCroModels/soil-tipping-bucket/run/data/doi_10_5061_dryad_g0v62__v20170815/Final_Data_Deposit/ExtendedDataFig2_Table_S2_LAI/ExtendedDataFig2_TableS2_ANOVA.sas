options ls=170 ps=10000 pageno=1;
title 'LAI 2004 through 2011';
Data LAI_2004_ringmeans;
Infile 'C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\2004 LAI ring means.csv' 
dlm=',' firstobs=2 lrecl=10000; input year CO2$ treatment$ ring block date$ DOY LAI ;;
run;
Data LAI_2005_ringmeans;
Infile 'C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\2005 LAI ring means.csv' 
dlm=',' firstobs=2 lrecl=10000; input year CO2$ treatment$ ring block date$ DOY LAI ;;
run;
Data LAI_2006_ringmeans;
Infile 'C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\2006 LAI ring means.csv' 
dlm=',' firstobs=2 lrecl=10000; input year CO2$ treatment$ ring block date$ DOY LAI ;;
run;
Data LAI_2007_ringmeans;
Infile 'C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\2007 LAI ring means.csv' 
dlm=',' firstobs=2 lrecl=10000; input year CO2$ treatment$ ring block date$ DOY LAI ;;
run;
Data LAI_2008_ringmeans;
Infile 'C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\2008 LAI ring means.csv' 
dlm=',' firstobs=2 lrecl=10000; input year CO2$ treatment$ ring block date$ DOY LAI ;;
run;
Data LAI_2009_ringmeans;
Infile 'C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\2009 LAI ring means.csv' 
dlm=',' firstobs=2 lrecl=10000; input year CO2$ DOY treatment$ block ring  LAI ;;
run; 
Data LAI_2010_ringmeans;
Infile 'C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\2010 LAI ring means.csv' 
dlm=',' firstobs=2 lrecl=10000; input year treatment$ date$ DOY time$  block ring  LAI ;;
run;
Data LAI_ringmeans_2004thru2011;
set LAI_2004_ringmeans LAI_2005_ringmeans LAI_2006_ringmeans LAI_2007_ringmeans LAI_2008_ringmeans LAI_2009_ringmeans
LAI_2010_ringmeans LAI_2011_ringmeans;
run;
data LAI_ringmeans_2004thru2011;
set LAI_ringmeans_2004thru2011;
if ring in (1 4 10 11 17 24 25 32) then CO2="ambient";
if ring in (3 5 14 15 20 21 29 28) then CO2="elevated";
if ring in (1 5) then block=1;
if ring in (3 4) then block=2;
if ring in (10 14) then block=3;
if ring in (11 15) then block=4;
if ring in (17 21) then block=5;
if ring in (20 24) then block=6;
if ring in (25 29) then block=7;
if ring in (28 32) then block=8;
drop treatment meter LAI_corr time;
run;
proc sort data =LAI_ringmeans_2004thru2011;
by year doy ring;
run ;
data LAI_ringmeans_2004;
set LAI_ringmeans_2004thru2011;
where year=2004;
run;
proc mixed data =LAI_ringmeans_2004 order =data covtest;
class block doy co2 ;
model LAI = co2|doy/ ddfm =kr;
random block ;
lsmeans co2|doy/ pdiff slice=doy; 
ods output diffs = pdiffs2004LAI;
ods output lsmeans = lsmeans2004LAI;
run ;

ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\SAS Output\2004LAIlsmeans.csv';
proc print data = lsmeans2004LAI;run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\SAS Output\2004LAIpdiffs.csv';
proc print data = pdiffs2004LAI; where DOY=_DOY; run; ODS CSV close;

data LAI_ringmeans_2005;
set LAI_ringmeans_2004thru2011;
where year=2005;
run;
proc mixed data =LAI_ringmeans_2005 order =data covtest;
class block doy co2 ;
model LAI = co2|doy/ ddfm =kr;
random block ;
lsmeans co2|doy/ pdiff slice=doy; 
ods output diffs = pdiffs2005LAI;
ods output lsmeans = lsmeans2005LAI;
run ;

ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\SAS Output\2005LAIlsmeans.csv';
proc print data = lsmeans2005LAI;run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\SAS Output\2005LAIpdiffs.csv';
proc print data = pdiffs2005LAI; where DOY=_DOY; run; ODS CSV close;

data LAI_ringmeans_2006;
set LAI_ringmeans_2004thru2011;
where year=2006;
run;
proc mixed data =LAI_ringmeans_2006 order =data covtest;
class block doy co2 ;
model LAI = co2|doy/ ddfm =kr;
random block ;
lsmeans co2|doy/ pdiff slice=doy; 
ods output diffs = pdiffs2006LAI;
ods output lsmeans = lsmeans2006LAI;
run ;

ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\SAS Output\2006LAIlsmeans.csv';
proc print data = lsmeans2006LAI;run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\SAS Output\2006LAIpdiffs.csv';
proc print data = pdiffs2006LAI; where DOY=_DOY; run; ODS CSV close;

data LAI_ringmeans_2007;
set LAI_ringmeans_2004thru2011;
where year=2007;
run;
proc mixed data =LAI_ringmeans_2007 order =data covtest;
class block doy co2 ;
model LAI = co2|doy/ ddfm =kr;
random block ;
lsmeans co2|doy/ pdiff slice=doy; 
ods output diffs = pdiffs2007LAI;
ods output lsmeans = lsmeans2007LAI;
run ;

ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\SAS Output\2007LAIlsmeans.csv';
proc print data = lsmeans2007LAI;run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\SAS Output\2007LAIpdiffs.csv';
proc print data = pdiffs2007LAI; where DOY=_DOY; run; ODS CSV close;

data LAI_ringmeans_2008;
set LAI_ringmeans_2004thru2011;
where year=2008;
run;
proc mixed data =LAI_ringmeans_2008 order =data covtest;
class block doy co2 ;
model LAI = co2|doy/ ddfm =kr;
random block ;
lsmeans co2|doy/ pdiff slice=doy; 
ods output diffs = pdiffs2008LAI;
ods output lsmeans = lsmeans2008LAI;
run ;

ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\SAS Output\2008LAIlsmeans.csv';
proc print data = lsmeans2008LAI;run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\SAS Output\2008LAIpdiffs.csv';
proc print data = pdiffs2008LAI; where DOY=_DOY; run; ODS CSV close;

data LAI_ringmeans_2009;
set LAI_ringmeans_2004thru2011;
where year=2009;
run;
proc mixed data =LAI_ringmeans_2009 order =data covtest;
class block doy co2 ;
model LAI = co2|doy/ ddfm =kr;
random block ;
lsmeans co2|doy/ pdiff slice=doy; 
ods output diffs = pdiffs2009LAI;
ods output lsmeans = lsmeans2009LAI;
run ;

ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\SAS Output\2009LAIlsmeans.csv';
proc print data = lsmeans2009LAI;run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\SAS Output\2009LAIpdiffs.csv';
proc print data = pdiffs2009LAI; where DOY=_DOY; run; ODS CSV close;


data LAI_ringmeans_2010;
set LAI_ringmeans_2004thru2011;
where year=2010;
run;
proc mixed data =LAI_ringmeans_2010 order =data covtest;
class block doy co2 ;
model LAI = co2|doy/ ddfm =kr;
random block ;
lsmeans co2|doy/ pdiff slice=doy; 
ods output diffs = pdiffs2010LAI;
ods output lsmeans = lsmeans2010LAI;
run ;

ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\SAS Output\2010LAIlsmeans.csv';
proc print data = lsmeans2010LAI;run;ODS CSV close;
ODS CSV file='C:\Users\sbgray\Documents\Soil Moisture Manuscript\Leaf Area Index\SAS Output\2010LAIpdiffs.csv';
proc print data = pdiffs2010LAI; where DOY=_DOY; run; ODS CSV close;


