# pr_corr.R
################################################################################
# Originally in the MATLAB script pr_corr.m
# From http://neumeier.perso.ch/matlab/waves.html
# % written by Urs Neumeier, 2003-2004
# % Modified from the Pcorr3.m function written by T. Mason, SOC, January 1997

# Adapted to R by Luke Miller et al January 2016
# Should probably determine the license status from U. Neumeier some day

## NOTE JAN 2016 THE FOLLOWING CODE IS NOT TESTED AT ALL

###############################################################################


################################################################################
################################################################################
#############################################
# function y<-wavenumL(f,h);
# % y<-wavenum(f,h): FUNCTION for the calculation of the wavenumber.
# %                   The dispersion relation is solved using a 
# %                   polynomial approximation.
# %                   f, wave frequency; f<-1/T.
# %                   h, water depth (in m).
# %
# %       George Voulgaris, SUDO, 1992
wavenumL <- function(f, h) {
	w <- 2*pi*f
	dum1 <- (w^2)*h/9.81
	dum2 <- dum1 + 
			(1.0+0.6522*dum1 + 0.4622*dum1^2 + 0.0864*dum1^4 + 
				0.0675*dum1^5)^(-1);
	dum3 <- sqrt(9.81*h*dum2^(-1)) / f
	y <- 2*pi*dum3^(-1)
	y # return y
}


################################################################################

# function H=pr_corr(pt,h,Fs,zpt,M,Corr_lim)
# % Input: PT = sea surface elevation time series (detrended)
# %        H = mean water depth (m)
# %        Fs = sampling frequency (Hz)
# %        Zpt = height of PT from seabed (m)
# %        M = length of segment (optional, default: 512)
# %        Corr_lim = [min max] frequency for attenuation correction (Hz, 
#                                     %    optional, default [0.05 0.33])                     

# % This function corrects a detrended sea surface time series for depth 
# % attenuation.
# % Frequency correction from 0.05 to 0.33 Hz. The program checks the length 
# % of the input and zero pads as necessary.  After correction, the output 
# % vector is truncated to the same length as the input.  
# %
# %    SURFACE = PR_CORR(PT,[],Fs,Zpt,M)
# %
# % Alternatively, if the second argument (H) is an empty matrix, PT should be 
# % the sea-surface above bottom time-series. In this case, each segment of PT 
# % will be linearly detrended, corrected for attenuation, and the linear trend 
# % added back before constructing the result SURFACE.
# %


pr_corr <- function(pt, Fs, zpt, M = 512, Corr_lim = c(0.05, 0.33) ){
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
	library(signal) # for hanning() function
	
	# normally the maximum attenuation correction should not be higher than 5
	max_attenuation_correction <- 5
	# % higher value to process boat-waves for Jean Ellis
	# % max_attenuation_correction <- 20; 
	
	# mininum frequency, below which no correction is applied (0.05) = 20 s period
	min_frequency <- Corr_lim[1]
	# maximum frequency, above which no correction is applied (0.33) = 3 s period
	max_frequency <- Corr_lim[2]
	
	H_with_NaN <- pt	# make a copy of the input depth values
	notNA <- which(!is.nan(pt)) # get the indices of good data (not NA)
	pt<-pt[notNA] # make a reduced copy with only the good rows of pt
	
	# do_detrend=isempty(h); # little h is the mean water depth. MATLAB version
	do_detrend <- TRUE # Need to insert code to check for presence of  
	# an argument 'h', h and skip the detrending if it
	# was provided. LPM note: h is not currently an argument
	
	m <- length(pt)	# m = number of samples (all the NAs are excluded)
	
	Noverlap <- M/2	# length of overlap of segments. M is an argument to the
	# function and has a default of 512
	N <- (ceiling(m/M)) * M # length of array, will be 
							# zero-padded to nearest multiple of M
	
	f <- c(NA, seq(1,(M/2),by = 1) * Fs/M)  
	# f will be a vector of frequencies, with a NA in the 1st position. 
	
	### LPM: I think we can initially ignore this section below, since we need to 
	### detrend the data. Skip down further.
# if ~do_detrend
#   K = wavenumL(f,h);  % calculates wave number for each frequency using 
#   % function defined below
#   Kpt=cosh(K*zpt)./cosh(K*h); % correction factor of spectrum for pressure
#   if isempty(pt), H=Kpt;end
#   Kpt(f<min_frequency | f>max_frequency) = 1;	% attenuation (0.05-0.33Hz only)
#   Kpt(Kpt < 1/max_attenuation_correction) = 1/max_attenuation_correction;
#   % correction factor never higher than max_attenuation_correction
#   % linear decrease of correction above f>max_frequency
#   fb_max=max(find(f<=max_frequency));
#   fKlin=[fb_max:min(fb_max+fix(length(K)/10),length(K))];
#   Kpt(fKlin)=(length(fKlin):-1:1)*(Kpt(fb_max)-1)/length(fKlin)  + 1; 
#   Kpt(1)=1;
#   Kpt(M:-1:M/2+2)=Kpt(2:M/2);			% second half of series is symetric  
#   if isempty(pt)
#   H=[f(2:end) 1./H(2:length(f)) 1./Kpt(2:length(f))];
#   fprintf([' The returned matrix contains in column 1 the frequencies,\n',...
#            ' in column 2 the theoretical correction factor, and in column 3 \n'...,
#            ' the effectively used correction factor, for a water depth of\n',...
#            ' %g m and PT position of %g m above the sea bed.\n'],h,zpt);
#   return
#   end
#   else
#     x=(1:M)';
#   end  
	
	x <- seq(1, M, by = 1) # vector of indices from 1 to M
	
	H <- vector(mode="numeric", length = N) # Make a vector of zeros, length N
	
	overlap_window <- hanning(M); # R version of MATLAB's hann() function,
								# returns a L-point symmetric Hann window that
								# is used as an array of coefficients to 
								# combine overlapping segments
	
	# Mirror the overlap values in the 2nd half of overlap_window
	overlap_window[((M/2)+1):length(overlap_window)] <- 1 - overlap_window[1:(M/2)]
	
# Step through segments of pt vector, detrend the data if needed
	for (q in seq(1, N-Noverlap,by = Noverlap)){
		o <- min(c(q+M-1, m))
		ptseg <- pt[q:o]
		seg_len <- length(ptseg)
		
		if (do_detrend){
			# Fit a simple linear regression through the segment of sea surface
			# height data and extract the intercept + slope
			trend <- coef(lm(ptseg ~  x[1:seg_len]))
			
			# Calculate a depth h at the midpoint of the segment of 
			# data using the regression coefficients stored in 'trend' 
			# (i.e. intercept and slope) 
			h <- trend[1] + ( trend[2] * ( (seg_len+1)/2 ) )  
			
			# Remove the linear trend from the segment of sea surface heights
			ptseg <- ptseg - (trend[1] + (trend[2]* x[1:seg_len])) 
			
			# Calculate the wave number for each frequency in f, using the
			# mean sea surface height h for this segment of data
			K <- wavenumL(f,h);  # wavenumL function defined separately
			
			# Correction factor of spectrum for pressure
			Kpt <- cosh(K*zpt) / cosh(K*h) 
			# Recall that zpt is supplied as an argument to the function, 
			# giving the height of the pressure transducer above the seabed.
			
			# Set values for frequencies outside the desired range to 
			# 1 so that no attenuation is applied 
			# (i.e. attenuate 0.05-0.33Hz only)
			Kpt[f < min_frequency | f > max_frequency] <- 1 	
			
			Kpt[Kpt < (1/max_attenuation_correction)] <- 1 / max_attenuation_correction	
			#   % correction factor never higher than max_attenuation_correction
			#   % linear decrease of correction above f>max_frequency
			
			fb_max <- max(which(f <= max_frequency))
			
			# matlab fix() function rounds toward zero, R trunc() function 
			# should substitute
			# The min() function will return the smaller of 
			# fb_max+fix(length(K)/10) or length(K)
			mymin <- min(c(fb_max+trunc(length(K)/10), length(K)))
			fKlin <- seq(fb_max,mymin, by = 1) 
			
			Kpt[fKlin] <- ((seq(length(fKlin),1,by= -1)) * 
						(Kpt[fb_max]-1)/length(fKlin)) + 1
		
			Kpt[1] <- 1 # Overwrites that NA that's been hanging out in row 1
			
			# 2nd half of series is symmetric
			Kpt[seq(M,((M/2)+2),by=-1)] <- Kpt[2:(M/2)] 	
		}  # End of if(do_detrend) statement
		
		if (seg_len < M){
			ptseg[(seg_len+1):M] <- 0 # pad ptseg with zeros if its length is 
									  # shorter than M
		}# end of if(seg_len < M)
		
		P <- fft(ptseg) # Calculate spectrum
		Pcor <- P / Kpt # Apply correction factor
		# The R fft(x, inverse=TRUE) function returns unnormalized series, so it
		# is necessary to divide by length(x)
		# The Re() function keeps only the real part of the result
		Hseg <- Re(fft(Pcor, inverse = TRUE) / length(Pcor)) 

		# Keep only the part that is the same length as the  original segment.
		Hseg <- Hseg[1:seg_len] 
		
		if (do_detrend){
			# Add linear regression values back on to the detrended data so that
			# they are once again trended.
			Hseg <- Hseg + (trend[1] + trend[2]*x[1:seg_len])
			# Hseg should now contain sea surface height values (meters)
		}
		
		# Insert the corrected sea surface heights back into the vector H
		# using the coefficients in overlap_window to taper the ends of the
		# the segment Hseg
		H[q:o] <- H[q:o] + Hseg * overlap_window[1:seg_len]
		
		# Handle the 1st segment in H 
		if (q == 1){
			H[1: min(c(Noverlap,seg_len))] <- Hseg[1: min(c(Noverlap,seg_len))]
		}
		
		# Deal with the tail end of H
		if ( ( (q+M) >= N) & (seg_len > Noverlap)){
			H[(q+Noverlap):o] <- Hseg[(Noverlap+1):length(Hseg)]
		}
		
	}   # End of for (q in seq(1,N-Noverlap,by=Noverlap)) loop
		
	H <- H[1:m] # truncate H back to original length m
	
	# Place the new corrected heights back into the
	# the positions in the original vector that had
	# good data. 
	H_with_NaN[notNA] <- H; 
	# And rename that to be H
	H <- H_with_NaN 
	
	H # return H with all of the depth-corrected surface heights and the 
	# missing values (NA) re-inserted. Should be the same length as the input
	# vector pt. Units should be the same as the original input values
	# in the vector pt (usually meters).
	
}  # end of p_corr function

###############################################################################
###############################################################################
# Function: detrend
# Simple function fits a straight line to a vector of values, and uses the
# regression coefficients to subtract off the linear trend from the values. 
detrend <- function(pt){
	# Determine length of pt
	seg_len <- length(pt)
	# Make a sequence of indices
	x <- seq(1,seg_len, by = 1)
	# Fit a simple linear regression through the segment of sea surface
	# height data and extract the intercept + slope
	trend <- coef(lm(pt ~  x[1:seg_len]))
	
	# Calculate a depth h at the midpoint of the segment of 
	# data using the regression coefficients stored in 'trend' 
	# (i.e. intercept and slope) 
	h <- trend[1] + ( trend[2] * ( (seg_len+1)/2 ) )  
	
	# Remove the linear trend from the segment of sea surface heights
	pt <- pt - (trend[1] + (trend[2]* x[1:seg_len])) 
	# Return the detrended vector
	pt
}


###############################################################################
###############################################################################
## Function: wavesp 
## **********CURRENTLY NONFUNCTIONAL + UNTESTED****************

## Spectral wave parameters from pressure transducer data (PT) after 
## correction of the pressure attenuation (see function pr_corr)
#
## Ported to R from MATLAB by Luke Miller 2016
## % written by Urs Neumeier, 2003
## % modified from Coast07.m written by T Mason, UPCE, 1999 (psd modified 7/99)


wavesp <- function(PT, Zpt, Fs){
	# Inputs
	# PT = vector of pressure transducer values, converted to seawater height
	#		(height of sea surface above pressure transducer, units = meters)
	# 
	
	min_frequency = 0.05;	#	% mininum frequency, below which no correction is applied (0.05)
	max_frequency = 0.33;		#% maximum frequency, above which no correction is applied (0.33)	
	
	#%*************************************************************
	#		% Prepare data for spectral analysis
			
	m <- length(PT);		# Get length of surface height record
	# matlab fix() function rounds toward zero, R trunc() function 
	# should substitute
#	M = fix(m/Noseg/2)*2;	#% length of segment
	M <- trunc(m/Noseg/2)*2;	#% length of segment
	h <- mean(PT);			#% mean water depth
	PT <- detrend(PT); # MATLAB detrend() function computes the least-squares 
	# fit of a straight line (or composite line for piecewise linear trends) to 
	# the data and subtracts the resulting function from the data. The 
	# detrend() function defined separately in this source file does the
	# equivalent operation. 
	
	#%***********************************************************%
		#	% Calculation of spetral density 
	myspec = spec.pgram(PT, detrend = FALSE)
	p = myspec$spec
	Snf = p
	# Produce the vector of frequencies in Hz from the spectrum object
	F = myspec$freq * 4
	
	# The spectrum function in MATLAB is doing a routine where it steps through
	# the input data (PT) in windows of length M (900 in Noseg=4) and 
	# calculating the spectrum using overlapping windows (overlap =450, total 
	# of 7 windows in a 3600 sample chunk when M = 900). The question is whether
	# it is worthwhile to implement something similar in R
	
#	[P,F] <- spectrum(PT,M,M/2,[],Fs);  # % gives spectrum from 0 to Fnyq
	# The above in MATLAB feeds an empty vector [] to the window argument,
	# which results in the spectrum() function using hanning(M) to generate
	# the hanning window. But the hanning() function in MATLAB does not 
	# included the zero-weighted samples at the start and end of the output
	# window, while the hann() function does. 
	
#	p <- P(2:end,1)*2/Fs;	#% normalization so that psd=cov(7/99), and remove first						
#	F(1) <- [];		#% element of P and F, which contains no useful information
			
#	% frequency range over which the spectrum is integrated for calculation of 
# the moments
	integmin=which.min(F[F >= 0])	 #% this influences Hm0 and other wave parameters
	integmax=which.max(F[F <= (max_frequency*1.5 )])
	
	# TODO: I stopped here, since the matlab and R routines are deviating 
# quite a bit at this point, and numbers aren't coming up similar at all. Take
# a look back at the spectrum stuff above, since that's where they really 
# deviate.
	df = F[1];	#	% bandwidth (Hz)
	moment = numeric(length = 7)
	for (i in -2:4){	#		% calculation of moments of spectrum
		moment[i+3]=sum( (F[integmin:integmax]^i) * Snf[integmin:integmax])*df
	}
	#%	fprintf('m%g  =  %g\n',i,moment(i+3));

m0 = moment[0+3]
	
}



