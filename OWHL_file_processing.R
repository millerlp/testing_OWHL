# OWHL_file_processing.R
# 
# Author: Luke Miller Jan 20, 2016
###############################################################################

require(oce)
require(RCurl)
myTimeZone = 'PST8PDT'
Sys.setenv(TZ = myTimeZone)

# The test file 15032600_halfday.CSV contains data from an early prototype
# deployed in approx 10m depth at Hopkins Marine Station in Pacific Grove CA
# on Mar 26, 2015. The sensor sat approximately 0.5m above the true sand bottom.
fname = "./testing_OWHL/sample_data/15032600_halfday.CSV"

# There is a little bit of mission info stored in the first row
missioninfo = scan(fname, what=character(),nlines = 1, sep = ',')

# To get the real headers, skip the first line
dat = read.csv(fname, skip = 1)
###########################
# Columns:
# POSIXt: elapsed seconds since 1970-01-01 00:00:00 (unix epoch) in whatever
#         timezone the sensor was set to during deployment
# DateTime: human-readable character date and time, in whatever timezone the
#         sensor was set to during deployment (Pacific Daylight Time for the
#         example file 15032600_halfday.CSV)
# frac.seconds: integer value representing the fractional seconds of the 
#         sample * 100. To get the actual fractional second, divide these values
#         by 100 (resulting in 0.0, 0.25, 0.50, 0.75)
# Pressure.mbar: temperature-compensated pressure in millibar, from the 
#         MS5803-14BA sensor. This is absolute pressure (including atmospheric)
# TempC:  internal temperature of the MS5803-14BA pressure sensor.
#########################

# Convert the DateTime column to a POSIXct object. The example data were
# collected with the clock set to Pacific Daylight Time
dat$DateTime = as.POSIXct(dat$DateTime, tz='PST8PDT') 
# Add on the fractional seconds for each row using the values in frac.seconds
dat$DateTime = dat$DateTime + (dat$frac.seconds/100)


###########
# Filter suspect time points. Suspect time points have ms values that aren't
# equal to 0, 25, 50, or 75. (there shouldn't be any in the sample data file)
# Make a matrix to hold the test results
filt = matrix(0,nrow = nrow(dat), ncol = 5)
filt[,1] = rep(0L,nrow(dat)) == dat[,'frac.seconds'] # compare every ms value to 0
filt[,2] = rep(25L,nrow(dat)) == dat[,'frac.seconds'] # compare every ms value to 25
filt[,3] = rep(50L,nrow(dat)) == dat[,'frac.seconds'] # compare every ms value to 50
filt[,4] = rep(75L,nrow(dat)) == dat[,'frac.seconds'] # compare every ms value to 75
# Sum every row of filt. If none of the tests above returned true, the result
# in the 5th column will be 0 (== FALSE). Otherwise it will be TRUE
filt[,5] = rowSums(filt[,1:4]) 
# Now find every row that has a suspect ms value
badrows = which(!(filt[,5]))
if (length(badrows)>0){
	# Remove any rows with suspect ms values
	dat = dat[-badrows,]	
	cat('Removed',length(badrows),'suspect rows from data.\n')
} else {
	cat('No troublesome data\n')
}

################################


################################
# Retrieve air pressure data from local tide station
station = 9413450 # 9413450 = Monterey, CA
# Generate start and end times for data retrieval based on the times in df
# Note that you have to replace any spaces in the timestamp with %20 codes
# so that the URL that is eventually generated has the correct html code for 
# a space. 
startT = strftime(dat$DateTime[1],format='%Y%m%d %H:%M')
startT = sub(" ","%20",startT)
endT = strftime(dat$DateTime[nrow(dat)],format='%Y%m%d %H:%M')
endT = sub(" ","%20",endT)
# Assemble a query, based on the API specifications described at 
# http://tidesandcurrents.noaa.gov/api/ 
query = 'http://tidesandcurrents.noaa.gov/api/datagetter?begin_date='
query = paste0(query,startT)
query = paste0(query,"&end_date=")
query = paste0(query,endT)
query = paste0(query,"&station=")
# You entered a station ID number at the top of this file, it gets used here
query = paste0(query,station)  
# Request air pressure
query = paste0(query,"&product=air_pressure&units=metric&time_zone=")
# Careful here, specify the correct desired time zone: gmt, lst, or lst_ldt
query = paste0(query,"lst_ldt")
query = paste0(query,"&format=csv")
# Send the complete URL to the CO-OPS website, requesting a csv file in return
# that can be parsed by the read.csv function.
ap = getURL(query) # Store as a text string
# Check the integrity of the response
if (regexpr("Bad", ap) != -1 | regexpr("Wrong", ap) != -1) {
	# Finding those text strings usually means the returned info is bad
	rm(ap)
	cat('Error retrieving air pressure data\n')
} else {
	# Otherwise the data are probably good, so convert them into a data frame
	airPressure = read.csv(textConnection(ap))
	# Convert the Date.Time column to a POSIXct value, using the correct time 
	# zone
	airPressure$Date.Time = as.POSIXct(airPressure$Date.Time, tz = myTimeZone)
	# Rename the pressure column to include units
	names(airPressure)[2] = "Pressure.mbar"
	# Drop the extra columns
	airPressure = airPressure[,-(c(3,4,5))]
}
meanAirPressure.mbar = mean(airPressure$Pressure.mbar, na.rm=TRUE)
# If the above won't work for you, the pressure at the site on Mar 26 2015 was
# approximately 1020 mbar. 
####################################################

###########################################
# NOTE: The following section is only necessary for this early test dataset
# from March 2015 that used a silicone tube over the sensor. Later units 
# should not have the same huge offset. 

# Generate a set of pressure readings that subtract off the extra pressure from
# the silicone tubing. The air pressure downloaded from a local weather station
# tells you the true sea level pressure, which can be compared to the readings 
# from the OWHL just before it was submerged in the ocean. Someday this should
# probably be converted into an automatic process or user-interactive process.
rawSeaLevPres = 1700  # OWHL Reading just before submergence
trueSeaLevPres = meanAirPressure.mbar # Obtained from local weather station
# Calculate a correction value based on the sea level air pressure from 
# the measured pressure in the OWHL prior to deployment.
#correc = rawSeaLevPres - slPres
correc = rawSeaLevPres - trueSeaLevPres
# Subtract off the correction value from all OWHL pressure readings, and 
# convert to bar. This should be the approximate absolute pressure, including 
# both atmospheric pressure and seawater pressure, without the offset introduced
# by the silcone tubing.
dat$Pressure.mbar = (dat$Pressure.mbar - correc)

# At this point all of the entries in dat$Pressure.mbar should be absolute
# pressure: the pressure due to atmospheric pressure at sea level plus 
# the pressure due to the seawater above the submerged OWHL, without the offset
# induced by the silicone tubing on the pressure port. 
###############################################################################
###############################################################################

# Convert pressure to gauge (seawater only) pressure by subtracting 
# the air pressure near the site on that date (downloaded earlier in the 
# script)
dat$swPressure.mbar = dat$Pressure.mbar - meanAirPressure.mbar

# Convert pressure to decibar so the oce swDepth function can be used
dat$swPressure.dbar = dat$swPressure.mbar / 100

mylatitude = 36.2 # Monterey, CA
# Use the oce function swDepth to estimate the height of seawater above the
# sensor.
dat$swDepth.m = swDepth(dat$swPressure.dbar, latitude = mylatitude)

# Data from the submerged portion of the data file begin around row 15,000
dat = dat[15000:nrow(dat),]

