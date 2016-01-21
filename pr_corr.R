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
# function y=wavenumL(f,h);
# % y=wavenum(f,h): FUNCTION for the calculation of the wavenumber.
# %                   The dispertion relation is solved using a 
# %                   polynomial approximation.
# %                   f, wave frequency; f=1/T.
# %                   h, water depth (in m).
# %
# %       George Voulgaris, SUDO, 1992
# echo off
wavenumL = function(f,h) {
	# f=f(:);
	w = 2*pi*f
	dum1 = (w^2)*h/9.81
	dum2 = dum1 + 
			(1.0+0.6522*dum1 + 0.4622*dum1^2 + 0.0864*dum1^4 + 0.0675*dum1^5)^(-1);
	dum3 = sqrt(9.81*h*(dum2^(-1)))/f
	y = 2*pi*(dum3^(-1))
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

# % This function corrects a detrended sea surface time series for depth attenuation.
# % Frequency correction from 0.05 to 0.33 Hz. The programme checks the length of the 
# % input and zero pads as necessary.  After correction, the output vector is truncated 
# % to the same length as the input.  
# %
# %    SURFACE = PR_CORR(PT,[],Fs,Zpt,M)
# %
# % Alternatively, if the second argument (H) is an empty matrix, PT should be the 
# % sea-surface above bottom time-series. In this case, each segment of PT will be 
# % linearly detrended, corrected for attenuation, and the linear trend added back before 
# % constructing the result SURFACE.
# %


pr_corr = function(pt, Fs, zpt, M = 512, Corr_lim = c(0.05, 0.33) ){
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
	
	
	# normally the maximum attenuation correction should not be higher than 5
	max_attenuation_correction = 5
	# % higher value to process boat-waves for Jean Ellis
	# % max_attenuation_correction = 20; 
	
	# mininum frequency, below which no correction is applied (0.05) = 20 s period
	min_frequency = Corr_lim[1]
	# maximum frequency, above which no correction is applied (0.33) = 3 s period
	max_frequency = Corr_lim[2]
	
	H_with_NaN = pt;	# make a copy of the input depth values
	not_NaN = which(!isnan(pt)); # get the indices of good data (not NA)
	pt=pt[not_NaN]; # make a reduced copy with only the good rows of pt
	pt_dimension=length(pt); # get length of the smaller pt
#	pt = as.matrix(pt, ncol = 1) # make sure pt is a column array
	# pt=pt(:);					% assures column array MATLAB specific
	
	# do_detrend=isempty(h); # little h is the mean water depth. MATLAB version
	do_detrend = TRUE # Need to insert code to check for presence of  
	# an argument 'h', h and skip the detrending if it
	# was provided. LPM note: h is not currently an argument
	
	m = length(pt)	# m = number of samples (all the NAs are excluded)
	
	Noverlap = M/2	# length of overlap of segments. M is an argument to the
	# function and has a default of 512
	N = (ceiling(m/M)) * M # length of array zero-padded to nearest multiple of M
	
	# f = [NaN; (1:M/2)'*Fs/M];          % frequencies column vector
	f = c(NA, seq(1,(M/2),by=1) * Fs/M)  # LPM: not sure why the NA needs to be included
	f = as.matrix(f, ncol = 1)
	
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
	
	### LPM: I think we need the x=(1:M)' for the next section, since that's the
	### only thing that happens in the if-else statement above if you need to 
	### detrend the data. It is just a vector of indices from 1 to M
	x = seq(1,M, by = 1)
	
	# H = zeros(N,1);	% pre-allocate vector # Making a N x 1 matrix of zeros
	H = matrix(0,nrow = N, ncol = 1) # Make a N x 1 matrix of zeros
	
#   overlap_window=hann(M);				% coefficient array to combine overlapping segments
	## w = hann(L) returns an L-point symmetric Hann window in the column vector w
	library(signal) # for hanning() function
	overlap_window = hanning(M); ## CHECK THAT THIS IS A SENSIBLE SUBSTITUTE
	
#   overlap_window(M/2+1:end)=1-overlap_window(1:M/2);
	overlap_window[((M/2)+1):nrow(overlap_window)] = 1 - overlap_window[1:(M/2)]
	
	
#   %overlap_window=triang(M);  # LPM: this was commented out in the original
#   
#   h_dif=[];  # LPM: can be ignored for now I think
	
# Step through segments of pt vector, detrend the data if needed, and do 
# some other stuff
	for (q in seq(1,N-Noverlap,by=Noverlap)){
		o = min(c(q+M-1, m))
		ptseq = pt[q:o]
		seg_len = length(ptseg)
		
		if (do_detrend){
			# trend=polyfit(x(1:seg_len),ptseg,1);	% calculates trend
			# MATLAB HELP: p = polyfit(x,y,n) returns the coefficients for a polynomial 
			# p(x) of degree n that is a best fit (in a least-squares sense) for the 
			# data in y. The coefficients in p are in descending powers, and the length 
			# of p is n+1
			# In the above case, it is a simple linear fit (degree = 1)
			# R equivalent: trend = coef(lm(y~x)), but coefficients are listed in 
			# ascending order, opposite of matlab
			trend = coef(lm(ptseg ~  x[1:seg_len]))
			
			# h=polyval(trend,(seg_len+1)/2);		    % mean water depth of segment
			# The above is calculating a depth h at the midpoint of the segment of 
			# data using the regression coefficients stored in 'trend' 
			# (i.e. intercept and slope) 
			h = trend[1] + ( trend[2] * ( (seg_len+1)/2 ) )  
			
			
			#   ptseg=ptseg-polyval(trend,x(1:seg_len));% remove trend
			ptseg = ptseg - (trend[1] + (trend[2]* x[1:seg_len])) # remove trend
			
			K = wavenumL(f,h); # % calculates wave number for each frequency using 
			# function defined earlier
			
			#   Kpt=cosh(K*zpt)./cosh(K*h);	 % correction factor of spectrum for pressure
			Kpt = cosh(K*zpt) / cosh(K*h) # correction factor of spectrum for pressure
			# recall that zpt is supplied as an argument to the function, giving the
			# height of the pressure transducer above the seabed.
			
			Kpt[f < min_frequency | f > max_frequency] = 1 # attenuation (0.05-0.33Hz only)	
			#   Kpt(Kpt < 1/max_attenuation_correction) = 1/max_attenuation_correction;
			Kpt[Kpt < (1/max_attenuation_correction)] = 1 / max_attenuation_correction	
			
#   % correction factor never higher than max_attenuation_correction
#   % linear decrease of correction above f>max_frequency
			
#   fb_max=max(find(f<=max_frequency));
			fb_max = which.max(f <= max_frequency)
			
#   fKlin=[fb_max:min(fb_max+fix(length(K)/10),length(K))];
			# matlab fix() function rounds toward zero, R trunc() function should 
			# substitute
			# The min() function will return the smaller of fb_max+fix(length(K)/10) or 
			# length(K)
			mymin = min(c(fb_max+trunc(length(K)/10), length(K)))
			fKlin = seq(fb_max,mymin, by = 1) # CHECK THIS RESULT
			
#   Kpt(fKlin)=(length(fKlin):-1:1)*(Kpt(fb_max)-1)/length(fKlin)  + 1;
			Kpt[fKlin] = ((seq(fKlin,1,by= -1)) * (Kpt[fb_max]-1)/length(fKlin)) + 1
			
#   Kpt(1)=1;
			Kpt[1] = 1
			
#   Kpt(M:-1:M/2+2)=Kpt(2:M/2);			% second half of series is symetric
			Kpt[seq(M,((M/2)+2),by=-1)] = Kpt[2:(M/2)] # 2nd half of series is symmetric	
		}  # End of if(do_detrend) statement
		if (seg_len < M){
			#   ptseg(M,1)=0;				% zero-pads to nearest length M 
			# LPM note: I may not be following what's happening in the above statement
			# since I think that statement would just drop a 0 in position (M,1) of
			# ptseg, which isn't the same as padding out to length M
		}# end of if(seg_len < M)
		
		P = fft(ptseq) # calculate spectrum. Need to check that R returns a format 
		# that works with the next operation. The R fft() function
		# returns the unnormalized univariate Fourier transform
		Pcor = P / Kpt # apply correction factor
#   Hseg = real(ifft(Pcor));	% corrected PT time series
		Hseg = Re(fft(Pcor, inverse = TRUE) / length(Pcor)) # CHECK ME CAREFULLY
		# The R fft(x, inverse=TRUE) function returns unnormalized series, so it
		# is necessary to divide by length(x)
		# The Re() function keeps only the real part of the result
		
#   Hseg=Hseg(1:seg_len);		% same length than original segment
		Hseg = Hseg[1:seg_len] # Keep only the part that is the same length as the
		# original segment.
		
		if (do_detrend){
#   Hseg=Hseg + polyval(trend,x(1:seg_len))% add linear regression values if data were detrended
			# Add linear regression values back on to the detrended data so that
			# they are once again trended.
			Hseg = Hseg + (trend[1] + trend[2]*x[1:seg_len])	
		}
#   H(q:o) = H(q:o) + Hseg.*overlap_window(1:seg_len);
		H[g:o] = H[g:o] + Hseg * overlap_window[1:seg_len]
		
		if (q == 1){
#   H(1:min(Noverlap,seg_len)) = Hseg(1:min(Noverlap,seg_len));
			H[1: min(c(Noverlap,seg_len))] = Hseg[1: min(c(Noverlap,seg_len))]
		}
		
		if ( ( (q+M) >= N) & (seg_len > Noverlap)){
#   H(q+Noverlap:o) = Hseg(Noverlap+1:end);
			H[(q+Noverlap):o] = Hseg[(Noverlap+1):length(Hseg)]
		}
		
		
	}   # End of for (q in seq(1,N-Noverlap,by=Noverlap)) loop
	
#   H = H(1:m);	
	H = H[1:m]
	
	##########################################
	## TODO: Need to implement the following in R
#   H = reshape(H,pt_dimension); % convert H into an array with dimensions 
	# % contained in pt_dimension
#   H_with_NaN(not_NaN) = H;	% add again NaN value that were present in pt
#   H = H_with_NaN;
	
	H # return H with all of the depth-corrected surface heights and the 
	# missing values (NA) re-inserted. Should be the same length as the input
	# vector pt.
	
}  # end of p_corr function
###############################################################################




