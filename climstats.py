import numpy as np
import scipy as sp
from scipy import stats
import inout as io
import rstats as rs

# Correlation
def cormap(txt_ts,ncarray,lag=0):
	"""input: txt_file... a text file string
	          ncfilein... an ncfile - string
	"""
	tdim = ncarray.shape[0]
	if ncarray.ndim==4:
		map_ts=ncarray[:,0,:,:]
	elif ncarray.ndim==3:
		map_ts=ncarray
	
	t,lat,lon = map_ts.shape
	cor_map = np.zeros((lat,lon))
	for i in xrange(lat):
		for j in xrange(lon):
			#slope, intercept, r_val, p_val, stderr = stats.linregress(map_ts[:,i,j],txt_ts[:,col])
			r_mat = np.corrcoef(map_ts[lag:tdim,i,j],txt_ts[0:tdim-lag,col])
			cor_map[i,j] = r_mat[1,0]  
			out_cor_map = np.nan_to_num(cor_map)
	return out_cor_map

# Regression
"""Maps the regression values (or other values) between a timeseries, e.g. NINO3 and an array, e.g. Hadisst
txt_ts is a text file of the timeseries, and col is the column to use (if there are more than one, otherwise leave blank)
statvar is the statistical variable: slope=0, intercept=1, r_val=2, p_val=3, stderr=4
Lag is in MONTHS 
"""
def regmap(txt_ts,data,lag=0,statvar=0):
	"""input:		txt_ts = 1D timeseries array 
							data   = 3D or 4D data array
	statvar=0 for regression, 2 for correlation
	"""
	#txt_ts = io.readtext(txt_file)
	#ncarray = io.readnc(ncfilein,ncvar)[0]
	### Note: program changed, originally the input was a string of a text or nc file, now it's an array
	###	Also removed 'col=0' from input, so txt array has to be 1D
	tdim=data.shape[0]
	if data.ndim==4:
		map_ts=data[:,0,:,:]
	elif data.ndim==3:
		map_ts=data
	
	t,lat,lon = map_ts.shape
	reg_map = np.zeros((lat,lon))
	for i in xrange(lat):
		for j in xrange(lon):
			#slope, intercept, r_val, p_val, stderr = stats.linregress(map_ts[:,i,j],txt_ts[:,col])
			#stattup = stats.linregress(map_ts[lag:1706,i,j],txt_ts[0:1706-lag,col])
			stattup = stats.linregress(txt_ts[0:tdim-lag],map_ts[lag:tdim,i,j])
			#stattup = stats.linregress(txt_ts[11+lag:1695,col],map_ts[11:1695-lag,i,j])
			reg_map[i,j] = stattup[0].copy()
			reg_II = np.nan_to_num(reg_map)
			out_reg = np.ma.masked_outside(reg_II,-1e9,1e9)
	return out_reg



def regmap_box(temp,var,stat=0,lag=0):
	"""
	Input a temperature array [time,lat,lon] and array of variable in question
	e.g. pressure [time,lat,lon].
	Outputs an array [lat,lon] of regression values, where each grid box is the temperature regessed with the variale
	stat=0: use 0 for regression, 2 for correlatino
	"""
	t,lat,lon = temp.shape
	reg_map = np.zeros((lat,lon))
	for i in xrange(lat):
		for j in xrange(lon):
			#slope, intercept, r_val, p_val, stderr = stats.linregress(map_ts[:,i,j],txt_ts[:,col])
			#stattup = stats.linregress(map_ts[lag:1706,i,j],txt_ts[0:1706-lag,col])
			stattup = stats.linregress(temp[lag:t,i,j],var[0:t-lag,i,j])
			#stattup = stats.linregress(txt_ts[11+lag:1695,col],map_ts[11:1695-lag,i,j])
			reg_map[i,j] = stattup[stat]
			reg_II = np.nan_to_num(reg_map)
			out_reg = np.ma.masked_outside(reg_II,-100000,100000)
	return out_reg



def crosscormap(txt_ts,data,maxlag=24):
	"""input:		txt_ts = 1D timeseries array (monthly data is best)
							data   = 3D or 4D data array
	This program will caculate the maximum cross-correlation, and lag time of that max, between a given timeseries and each gridbox
	"""
	tdim=data.shape[0]
	if data.ndim==4:
		map_ts=data[:,0,:,:]
	elif data.ndim==3:
		map_ts=data

	t,lat,lon = map_ts.shape
	maxcor_map = np.zeros((lat,lon))
	maxcorlag_map = np.zeros((lat,lon))
	for i in xrange(lat):
		for j in xrange(lon):
			xcv, lags=rs.xcorr(txt_ts[0:tdim],ts2=map_ts[0:tdim,i,j],maxlag=maxlag)
			maxcor=abs(xcv).max()
			maxcorlag=np.array([abs(np.argmax(xcv)-maxlag),abs(np.argmin(xcv)-maxlag)]).min()
			maxcor_map[i,j] = np.nan_to_num(maxcor.copy())
			maxcorlag_map[i,j] = np.nan_to_num(maxcorlag.copy())
			#out_reg = np.ma.masked_outside(reg_II,-1e-9,1e9)
	return maxcor_map, maxcorlag_map







