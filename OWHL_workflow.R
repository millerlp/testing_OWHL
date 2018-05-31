# OWHL_workflow.R
# 
# Author: Luke Miller May 30 2018
###############################################################################

# To install the oceanwaves package for estimating wave height and period:
install.packages('devtools')
library(devtools)
install_github('millerlp/oceanwaves')


setwd("~/R_public/oceanwaves") # package must be the working directory
library(devtools)  # load devtools
devtools::document() # Regenerate documents, help files, namespace etc
build() # Generate a tar.gz file of the entire package
install() # Re-install oceanwaves after an update
load_all() # Actually load the functions in oceanwaves 
setwd("~/R_public/testing_OWHL") # Change back to the normal directory

library(zoo)

# Manually enter a sea level air pressure (mbar) for the deployment date in 
# order to correct for any residual pressure signal on the sensor. This could
# be replaced with a value from a local weather station for the time of 
# the initial deployment
initialSurfacePressure = 1014
# Set the R environment's time zone to UTC
myTimeZone = 'UTC'
Sys.setenv(TZ = myTimeZone)
# Define a surface reading for OWHL and start of subsurface deployment
# These are values specific to Kristen's first example data set
owhlSurfaceTime = as.POSIXct('2016-08-19 18:15', tz = myTimeZone)
owhlDeployStartTime = as.POSIXct('2016-08-19 19:15', tz = myTimeZone)
owhlDeployEndTime = as.POSIXct('2016-09-27 18:15', tz = myTimeZone)

# Assume we're starting with a directory of daily csv data files from the OWHL
filenames = dir('~/Dropbox/OWHL_misc/Deployment_Elsmore_Marguerite_201608_short/',
                pattern = "*.csv", full.names=TRUE)

# Create a text progress bar for the long-running data ingestion
pb = txtProgressBar(min=0,max = length(filenames), style = 3)
# Cycle through each daily data file and concatenate it onto the 'dat' data frame
for (i in 1:length(filenames)){
  # There is a little bit of mission info stored in the first row
  missioninfo = scan(filenames[i], what=character(),nlines = 1, sep = ',',
                     quiet = TRUE)
  # To get the real headers, skip the first line
  temp = read.csv(filenames[i], skip = 1)
  if (i == 1){
    dat = temp
  } else {
    # add on the subsequent files
    dat = rbind(dat,temp)
  }
  setTxtProgressBar(pb,i) # update text progress bar
}
close(pb) 	# shut down text progress bar

# Convert the DateTime column to a POSIXct object. For now we'll assume the 
# OWHL clock was set to UTC time zone. This could take a few minutes if there's
# lots of data
dat$DateTime = as.POSIXct(dat$DateTime, tz='UTC') 
# Add on the fractional seconds for each row using the values in frac.seconds
dat$DateTime = dat$DateTime + (dat$frac.seconds/100)

# Order all of the data based on the DateTime column
dat = dat[order(dat$DateTime),]

#################################
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
#########################################
# If there are brief bouts of missing data, no longer than 1 second, impute
# missing values to generate a continuous time series
alldiffs = diff(dat$DateTime)
gaps = which((alldiffs > 0.25))
# Find locations of any larger gaps
largegaps = which(alldiffs > 1.25)
# Remove any larger gaps from the gaps vector
if (length(largegaps)>0){
  matches = which(gaps == largegaps)
  if (length(matches > 0)){
    # Remove indices of large gaps
    gaps = gaps[-matches]
  }
}

# Go through each index value in gaps and impute the missing rows
if (length(gaps)> 0){
for (i in 1:length(gaps)){
  # Determine the length of this time gap
  timediff = dat[gaps[i]+1,'POSIXt'] - dat[gaps[i],'POSIXt']
  # If 2 timestamps are 2 seconds apart, then there must be a missing 
  # second of data in between them
  if (timediff == 2){
    # One second of data is missing
    if (i == 1){
      # Set aside the first chunk of data
      dat1 = dat[1:gaps[i],]
    } else {
      # Get the chunk from the last gap to the current gap and add it to dat1
      tempdat = dat[(gaps[i-1]+1):gaps[i],]
      dat1 = rbind(dat1,tempdat)
    }
    # Calculate how many rows to insert
    timesteps = (timediff-1) * 4 # 4 samples per second
    # Generate a set of time stamps for the missing second as a data frame
    tempdataframe = data.frame(POSIXt = rep(dat$POSIXt[gaps[i]]+1,timesteps),
                               DateTime = rep(as.POSIXct(dat$POSIXt[gaps[i]]+1,
                                                         tz=myTimeZone,
                                                         origin = '1970-1-1'),timesteps) +
                                 c(0,0.25,0.5,0.75),
                               frac.seconds = c(0,25,50,75),
                               Pressure.mbar = rep(NA,timesteps),
                               TempC = rep(NA,timesteps) )
    # Add on the newly created chunk of time stamps to 'dat1'
    dat1 = rbind(dat1,tempdataframe)
  }
}
  # Finish by adding the final rows of 'dat' onto the expanded 'dat1' data frame
  tempdat = dat[(gaps[i]+1):nrow(dat),]
  dat1 = rbind(dat1,tempdat)
  
  # dat1 should now have no 1-second gaps in the timestamps, but could have
  # missing pressure and temperature values where timestamps had to be inserted
  # Use the na.approx function from the zoo package to do linear interpolation
  # of the missing pressure and temperature values for the 1 second gaps
  pressure = zoo(dat1$Pressure.mbar,order.by = dat1$DateTime)
  pressure2 = as.numeric(na.approx(pressure))
  dat1$Pressure.mbar = pressure2
  temperature = zoo(dat1$TempC, order.by = dat1$DateTime)
  temperature2 = as.numeric(na.approx(temperature))
  dat1$TempC = temperature2
  
  # Now that dat1 should have all the short gaps filled, rewrite it to dat
  dat = dat1
}

#################################################################
# Grab an average surface pressure reading for the OWHL unit
surfaceIndx = which.min(abs(dat$DateTime - owhlSurfaceTime))

# Grab 10 seconds of surface pressure readings. This will be used to determine
# any pressure offset in the sensor, when compared to local sea surface air
# pressure at the time of deployment.
surfaceOWHLPress = dat$Pressure.mbar[surfaceIndx:(surfaceIndx+40)]
surfaceOWHLPress = round(mean(surfaceOWHLPress),dig=1)
# Generate a set of pressure readings that subtract off any extra pressure from
# the sensor setup. The air pressure downloaded from a local weather station
# tells you the true sea level pressure, which can be compared to the readings 
# from the OWHL just before it was submerged in the ocean. Someday this should
# probably be converted into an automatic process or user-interactive process.
rawSeaLevPres = surfaceOWHLPress  # OWHL Reading just before submergence
trueSeaLevPres = initialSurfacePressure # Obtained from local wave buoy
# Calculate a correction value based on the sea level air pressure from 
# the measured pressure in the OWHL prior to deployment.
correc = rawSeaLevPres - trueSeaLevPres  # units of millibars
# Subtract off the correction value from all OWHL pressure readings. 
# This should be the approximate absolute pressure, including 
# both atmospheric pressure and seawater pressure, without the offset introduced
# by any of the protective tubing attached to the sensor. 
dat$Pressure.mbar.corr = (dat$Pressure.mbar - correc)

# Subtract off mean air pressure from all pressure readings to convert from
# absolute pressure to gauge pressure (pressure due to seawater only). It would
# be nice to do this correction on an hourly interval since we have the buoy
# data to do it (but need to align the owhl/buoy data time values to accomplish
# this). 
dat$swPressure.mbar = dat$Pressure.mbar.corr - initialSurfacePressure

# Convert pressure to decibar so the oce package swDepth function can be used
dat$swPressure.dbar = dat$swPressure.mbar / 100

mylatitude = 33.72 # Los Angeles
# Use the oce function swDepth to estimate the height of seawater above the
# sensor.
dat$swDepth.m = oce::swDepth(dat$swPressure.dbar, latitude = mylatitude)

# Remove all surface data prior to deployment and after end
dat = dat[dat$DateTime >= owhlDeployStartTime,]
dat = dat[dat$DateTime <= owhlDeployEndTime,]

################################################################################
# Seawater depth values can be input into the oceanwaves::prCorr function to  
# correct for the depth attenuation of the pressure signal at the benthos. 
###### 
# Inputs
# pt	- A vector of surface height values (meters) derived from the 
#			original pressure sensor time series (uncorrected for depth
#			attenuation effects)
# Fs	- Sampling frequency (Hz). Normally 4 Hz for OWHL logger
# zpt	- Height of pressure sensor above seabed (meters)
# M		- Length of time series segments that will be used in the 
#			detrending and attenuation correction operations. 512 samples
#			is the default. Should be an even number.
# Corr_lim - [min max] frequency for attenuation correction (Hz, 
#                optional, default [0.05 0.33])


dat$swDepthcorr.m = NA
# In some cases the pressure data were collected in chunks of time rather 
# than continuously (i.e. from a Seabird, rather than from OWHL). 
# First determine the boundaries of each contiguous chunk by looking for 
# time steps between rows that are greater than the normal 0.25 second step
# If there are not breaks due to missing data, the if/else statement below
# will handle that case and instead break the continuous time series into 15
# minute chunks
bounds = which(diff(dat$DateTime) > 0.25)
if (length(bounds) == 0){
  stepSize = 15 * 60 * 4 # 15 minutes x 60 seconds x 4 Hz
  bounds = seq(from = 0, to = nrow(dat), by = stepSize)
} else {
  bounds = c(0,bounds) # Add on the first row index as well
}

cat("Applying pressure attenuation correction\n")
pb = txtProgressBar(min=0,max = length(bounds), style = 3)
for (i in 2:length(bounds)){
  # Apply pr_corr function to the swDepth data in the 15-minute chunks
  # The segment ((bounds[i-1])+1) extracts the previous row index from bounds
  # which is the end of the chunk, so you want to add 1 row to the value so
  # that you arrive at the start of the next good chunk. The value of 
  # bounds[i] is the end of the good chunk. 
  dat$swDepthcorr.m[((bounds[i-1])+1):bounds[i]] = 
    oceanwaves::prCorr(dat$swDepth.m[((bounds[i-1])+1):bounds[i]], Fs = 4, 
                       zpt = 0.1, M = 512)
  setTxtProgressBar(pb,i) # update text progress bar
}
close(pb)	# shut down text progress bar

# If prCorr works correctly, the values in swDepthcorr.m should be sea 
# surface height in meters, after accounting for the depth attenuation of the
# OWHL pressure signal that is unavoidable in bottom-mounted pressure 
# transducers (as opposed to surface buoys). Generally these corrected sea
# surface height values will be slightly more extreme than the raw sea water
# depth calculated from the pressure signal. 


# Go through each 15 minute time chunk and calculate the significant wave
# height of the de-trended data (removing tide signal). Calculate the peak
# wave period as well. 
# A 15 minute chunk at 4Hz should be 3600 rows
pb = txtProgressBar(min=0,max = length(bounds), style = 3)
for (i in 2:length(bounds)){
  # Extract the chunk of contiguous data
  tempdat = dat[((bounds[i-1])+1):bounds[i],]
  # Extract the date time in the middle of the chunk
  halfwidth = (bounds[i]-bounds[i-1])/2
  timeindex = bounds[i-1] + halfwidth
  tempDateTime = dat[timeindex,'DateTime']
  
  # Using the zero-crossing function to estimate sig. wave height. 
  zcstats = oceanwaves::waveStatsZC(tempdat$swDepthcorr.m, Fs = 4)
  # Extract the peak period
  Period.zerocross = zcstats$Tsig
  # Extract the significant wave height
  Hsig.zerocross = zcstats$Hsig
  ####
  # Now use the spectral analysis method instead
  spstats = oceanwaves::waveStatsSP(tempdat$swDepthcorr.m, Fs = 4)
  # Extract the significant wave height estimate
  Hsig.spectral = spstats$Hm0 
  # Extract the peak period estimate
  Period.spectral = spstats$Tp
  
  # Stick the wave heights and periods in an output data frame 'results'
  if (i == 2){
    results = data.frame(DateTime = tempDateTime,
                         Period.zerocross = Period.zerocross,
                         Hsig.zerocross = Hsig.zerocross,
                         Period.spectral = Period.spectral,
                         Hsig.spectral = Hsig.spectral)
  } else {
    tempresults = data.frame(DateTime = tempDateTime,
                             Period.zerocross = Period.zerocross,
                             Hsig.zerocross = Hsig.zerocross,
                             Period.spectral = Period.spectral,
                             Hsig.spectral = Hsig.spectral)
    results = rbind(results,tempresults)
  }
  setTxtProgressBar(pb,i) # update text progress bar
  
}
close(pb)

# Show the results data frame, showing 2 types of estimates of significant
# wave height and peak period. 
results

