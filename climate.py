import numpy as np
import scipy as sp
from numpy import ma
from scipy import stats
import inout as io
#import rstats as rs

def masklatlon(lonmin=0,lonmax=360,latmin=-90,latmax=90,t=None,had=None):
	""" Create a mask (1s = True, 0s = False) for a given latitude and longitude 
	Currently only suitable for 73x96 size grid"""

	if had==None:
		lon = np.linspace(0,356.25,96)
		lat = np.linspace(-90,90,73)
	elif had=='hadisst':
		lon = np.linspace(-179.5,179.5,360)
		lat = np.linspace(-89.5,89.5,180)
	elif had=='hadcrut':
		lon = np.linspace(-177.5,177.5,72)
		lat = np.linspace(-87.5,87.5,36)
	elif had=='rghadcrut':
		lon = np.linspace(0,355,72)
		lat = np.linspace(-87.5,87.5,36)
	elif had=='cmip3':
		lon = np.linspace(0,357.5,144)
		lat = np.linspace(-88.75,88.75,72)
	elif had=='ncep':
		lon = np.linspace(0,357.5,144)
		lat = np.linspace(90,-90,73)
	elif had=='plev':
		lon = np.linspace(1.875,358.125,96)
		lat = np.linspace(-88.75,88.75,72)

	LON,LAT = np.meshgrid(lon,lat)
	mask=np.array((LON>=lonmin)*(LON<=lonmax)*(LAT>=latmin)*(LAT<=latmax),dtype=int) # Use <,> or <=, >=. greater or equal seems better ??
	if t==None:
		maskout=mask
	else:
		maskout = np.repeat(np.expand_dims(mask,axis=0),t,0)
	return maskout

def maskdata(data,mask):
  """ Create a masked array using the weighted data
  reversing the mask (~) """
  out_data=np.ma.masked_array(data,~mask)
  #out_data=maskeddata.mask*maskeddata.data
  return out_data

def add_anom(climatology, anomaly, scale=1, t=12):
  # Adds an anomaly to a climatology
  #repdim_anom = np.repeat(np.expand_dims(anomaly,axis=0),12,0)
  dims=anomaly.ndim
  cldims=climatology.ndim
  if dims==cldims:
    if dims==3:
      ancil=np.zeros([t,73,96])
      ancil=climatology+anomaly*scale
    elif dims==4:
      ancil=np.zeros([t,1,73,96])
      ancil=climatology+anomaly*scale
  elif dims==3 and cldims==4:
    ancil=np.zeros([t,1,73,96])
    exdanom=np.expand_dims(anomaly,axis=1)
    ancil = climatology + exdanom*scale
  elif dims==4 and cldims==3:
    ancil=np.zeros([t,73,96])
    anom2=anomaly[:,0,:,:]
    ancil = climatology + anom2*scale
    
  return ancil

def landocmask(icenc,lsmasknc):
	""" generate an icefree ocean mask, and spits out a land mask for good measure
	Input: nc file of ice concentration, land-sea mask (as strings)
	Output: ice free ocean mask (1s and 0s)
	To use: ocean = ma.masked_where(icexp<1,tempsfc)
	example: 
	"""
	lsmask=io.readnc(lsmasknc,'lsm')[0][0,0,:,:]
	#Create an icefree ocean mask
	ice=io.readnc(icenc,'iceconc')[0]
	icesum=ice.sum(axis=0)
	icelsm=icesum+lsmask
	ifmask = ma.masked_greater(icelsm,0)+1
	icefree = ifmask.data*~ifmask.mask
	sansice = icefree[0,:,:]
	#ocean = ma.masked_where(icexp<1,tempsfc)
	return lsmask, sansice

def weightave(data, mask=None,had=None):
	""" Assume the earth is a perfect sphere (close enough!), 
	calculate the weighting for a given latitude and create a
	weighting array, then use that to weight the data """
	if data.ndim == 4:
		data3 = data[:,0,:,:]
	elif data.ndim ==3:
		data3 = data.copy()
	#if mask != None:
	#	if mask.ndim == 4:
	#		mask3 = mask[:,0,:,:]
	#	elif mask.ndim ==3:
	#		mask3 = mask.copy()
	  
	# From Thom Chubb; create a grid of lat/lon values, only works for 96x73
	if had==None:
		lon = np.linspace(0,356.25,96)
		lat = np.linspace(-90,90,73)
	elif had=='hadisst':
		lon = np.linspace(-179.5,179.5,360)
		lat = np.linspace(-89.5,89.5,180)
	elif had=='hadcrut':
		lon = np.linspace(-177.5,177.5,72)
		lat = np.linspace(-87.5,87.5,36)
	elif had=='rghadcrut':
		lon = np.linspace(0,355,72)
		lat = np.linspace(-87.5,87.5,36)
	elif had=='cmip3':
		lon = np.linspace(0,357.5,144)
		lat = np.linspace(-88.75,88.75,72)
	elif had=='ncep':
		lon = np.linspace(0,357.5,144)
		lat = np.linspace(90,-90,73)
	elif had=='plev':
		lon = np.linspace(1.875,358.125,96)
		lat = np.linspace(-88.75,88.75,72)

	LON,LAT = np.meshgrid(lon,lat)
	cosLAT = np.cos(LAT*np.pi/180)
	coslat = np.cos(lat*np.pi/180)
	#clext = np.repeat(np.expand_dims(coslat,axis=0),tdim,0)
	tdim = data3.shape[0]
	
	if had==None:
		lat_dist=2.5*111.319
		lon_dist=3.75*111.319
	elif had=='hadisst':
		lat_dist=1*111.319
		lon_dist=1*111.319
	elif had == 'hadcrut':
		lat_dist=5*111.319
		lon_dist=5*111.319
	elif had == 'rghadcrut':
		lat_dist=5*111.319
		lon_dist=5*111.319
	elif had == 'cmip3':
		lat_dist=2.5*111.319
		lon_dist=2.5*111.319
	elif had == 'ncep':
		lat_dist=2.5*111.319
		lon_dist=2.5*111.319
	elif had=='plev':
		lat_dist=2.5*111.319
		lon_dist=3.75*111.319
	area=np.zeros(data3[0,:,:].shape)
	for i in xrange(len(lat)):
		area[i,:] = lat_dist*lon_dist*coslat[i]
	if mask != None:
		area_mask = ma.masked_where(mask<1.0,area)
	elif mask == None:
		area_mask = area

	w_mean=np.zeros((tdim))
	for k in xrange(tdim):
		w_mean[k] = ma.average(data3[k,:,:],weights=area_mask)

	#merid_mean = np.sum(data3*area,axis=1)/np.sum(area,axis=1)
	#merid_area = area_mask.mean(axis=1)
	#w_mean = np.sum(merid_mean*merid_area,axis=-1)/np.sum(merid_area,axis=-1)
	#out_data = clext * data
	return w_mean





def region_ave(data,genmask=None,region=None,latlon=[-90,90,0,360],lomask=None,had=None):
	""" Gives the weighted average for specified region, or for specified latlon box of land or ocean
	Input: 
	data = [time, (height), lat, lon]
	region = 'nino3', 'tropics',
	latlon = [latmin, latmax, lonmin, lonmax]
	landsea = 'land', 'sea'
	lsmask = 1s for Land, 0s for Ocean
	ocmask = 1s for Ocean, 0s elsewhere (icefree ocean)
	Output: 1D array

	For all masks, 1s will be kept data, 0s excluded
	"""

	if data.ndim == 4:
		data0 = data[:,0,:,:]
	elif data.ndim == 3:
		data0 = data.copy()
	tdim = data.shape[0]
	xys = data0.shape[1:]
	#print(xys)
	if region=='nino3':
		if had==None:
			latlon = [-5,5,210,270] #for data from 0 - 360 degress
		else:
			latlon = [-5,5,-150,-90] # for data from -180 - 180 degrees
	
	if region=='tropics':
		latlon = [-23.5,23.5,0,360]
	
	latlonmask = masklatlon(latmin=latlon[0],latmax=latlon[1],lonmin=latlon[2],lonmax=latlon[3],had=had)
	#print('latlon: '+str(latlonmask.shape))
	#print('tdim: '+str(tdim))
	onesies=np.ones(xys)
	if lomask!=None:
		#maskout = np.repeat(np.expand_dims(mask,axis=0),t,0)
		data_mask = ma.masked_where(lomask*latlonmask<1,onesies)
	else:
		data_mask = ma.masked_where(latlonmask<1,onesies)
	if genmask!=None:
		data_mask = ma.masked_where(genmask<1,onesies)

	masked_ave = weightave(data0,mask=data_mask,had=had)

	return masked_ave, data_mask


#===================================================================================



def troposave(data, startlev=0, endlev=11,weights=False):
	# pressure in Pa
	plev = np.array([100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000,\
			20000, 15000, 10000, 7000, 5000, 3000, 2000])
	# geometr. height in m
	height = np.array([110.9, 762.1, 1457.6, 3013.6, 4209.2, 5579.3, 7193.6, 9177.2, 10379.9,\
		11805.9, 13637.6, 16221, 18495.3, 20643, 23938.4, 26592.3])
	# temperature in C
	temperature = np.array([14.28, 10.05, 5.53, -4.58, -12.34, -21.23, -31.71, -44.57, -52.36,\
		-56.5, -56.5, -56.5, -56.5, -55.92, -52.65, -50.02])
	# density in kg/cm
	density = np.array([1.212, 1.1378, 1.0625, 0.90796, 0.80142, 0.69142, 0.57713, 0.4572, 0.39445,\
		0.32159, 0.24119, 0.16079, 0.11256, 8.0184e-2, 4.7397e-2, 3.1223e-2])
	# Calulate a mass weighted tropospheric mean


	meanT = (data[:,startlev+1:endlev,::]+data[:,startlev:endlev-1,::])*0.5

	meanRho = (density[startlev+1:endlev]+density[startlev:endlev-1])*0.5

	diffH = np.diff(height[startlev:endlev])
	massw = meanRho * diffH
	w_ave = np.average(meanT,axis=1,weights=massw)
	if weights == True:
		outdata = [w_ave, massw]
	else:
		outdata = w_ave
	return outdata














