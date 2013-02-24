import numpy as np
import inout as io
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as rpyn

#x=io.readtext('txtfiles/temp.sfc.4ysl.nino.rsc.txt')
#y=io.readtext('txtfiles/temp.sfc.4ysl.land.rsc.txt')

def xcorr(ts1,ts2=None,maxlag=72,freq=12):
	""" Replicates the matlab function xcorr using rpy2 (which may need to be installed)
	Input: 		ts1 and ts2 are timeseries, 1D arrays. If only one array is entered the the auto-correlation will be calculated.
						maxlag... the maximum lag values
						freq: For monthly data use freq=12, for annual I guess freq = 1
	Output:	cor_out:  The correlation values at each lag (numpy array)
						lags:			An array with the values of the lags. 
	To match the matlab output and make it easier to plot cross and auto-correlations together, the auto-correlations are 'mirrored' for <0.
	"""

	#define R functions
	rts=robjects.r['ts']		# R function used to create timeseries
	rccf=robjects.r['ccf']	# R function to calculate cross-correlations
	racf=robjects.r['acf']	# R function to calculate auto-correlations

	#Convert python array to an R vector (Floatvector), then an R timeseries
	ts1_r=rts(robjects.FloatVector(ts1),frequency=freq)

	if ts2==None:
		# Use autocorrelation if there's only one timeseries
		acf_ts1 = racf(ts1_r,lag_max=maxlag,plot=False)
		ac_ts1=rpyn.ri2numpy(acf_ts1[0])[:,0,0]										# Converts R array back to numpy array
		cor_out=np.concatenate([ac_ts1[::-1],ac_ts1[1:maxlag+1]])	# This mirrors the positive values of the auto-corr

	elif ts2!=None:
		ts2_r=rts(robjects.FloatVector(ts2),frequency=freq)
		ccf_ts12 = rccf(ts1_r,ts2_r,lag_max=maxlag,plot=False)
		cc_ts12=rpyn.ri2numpy(ccf_ts12[0])[:,0,0]
		cor_out=cc_ts12

	lags=np.concatenate([np.linspace(-maxlag,-1,num=maxlag),np.linspace(0,maxlag,num=maxlag+1)])

	return cor_out, lags

# http://stat.ethz.ch/R-manual/R-patched/library/stats/html/acf.html
# http://oceansciencehack.blogspot.com.au/2010/05/r_15.html






