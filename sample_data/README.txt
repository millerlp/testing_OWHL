The file 15032600_halfday.CSV contains pressure and temperature data from an early prototype of the Open Wave Height Logger. The unit was deployed Mar 26 2015 in approximately 10 meters of water off Hopkins Marine Station in Pacific Grove, CA (36.2 degrees latitude). Time stamps are in Pacific Daylight Time zone. 

The dataset contains measurements prior to the start of the deployment, and the submergence begins around 12pm in the data file (around roughly row 15,000). 

The file contains 1 extra header line at the top. Column headers are in row 2.

Columns
POSIXt	- elapsed seconds since 1970-01-01 00:00 in Pacific Daylight Time zone
DateTime - human-readable time stamp YYYY-MM-DD HH:MM:SS in Pacific Daylight Time zone
frac.second - representation of the fractional seconds * 100. So a value of 25 = 0.25 seconds, a value of 50 = 0.5 seconds, a value of 75 = 0.75 seconds. Used to add to the timestamp columns to create sub-second timestamps. Sampling rate was 4Hz
Pressure.mbar -  Absolute pressure in millibar from the MS5803-14 sensor 
TempC	- temperature from the MS5803-14 pressure sensor internal sensor, degrees C.  

NOTE - this prototype had a silicone tube over the pressure sensor, which imparted a substantial offset (approx +700mbar above ambient pressure), evidenced by the readings around 1700mbar at the start of the data file when the unit was still at the surface. Atmospheric pressure at sealevel on this day was 1019-1020 mbar. 