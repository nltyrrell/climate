import numpy as np
import scipy as sp
from scipy import stats
import Scientific.IO.NetCDF as S
#import matplotlib.pyplot as plt
#import mpl_toolkits.basemap as bm

"""
@author Nicholas Tyrrell
@brief Get a weighted area average

Input:  NetCDF file, one 2D variable, area
Output: Text file timeseries

"""
def algorhythm_ave(lonmin,lonmax,latmin,latmax,areaname):
  """ What would Al Gore's band be called? Algorhythm.
  This just calls all the functions in order """
  nc_data, ncobj=readnetcdf('ncfiles/had.rsc.nc','sst')
  #hadisst.nc','sst')
  w_data = weighting(nc_data)
  print w_data.shape
  maskext = masklatlon(lonmin,lonmax,latmin,latmax)
  m_data = maskdata(w_data,maskext)
  ave_data = weightedave(m_data)
  print ave_data.shape
  writetxt(ave_data,areaname)
  return ave_data

# Import netcdf file as array
def readnetcdf(nc_file,variable):
  """ Function the reads in an nc file and outputs a numpy array """
  nc_fileobj = S.NetCDFFile(nc_file,mode='r')
  nc_data = nc_fileobj.variables[variable].getValue()
  return nc_data, nc_fileobj

"""
i)    Create weighting - coslat
      Weight data
ii)   Create mask (lat x lon)
      extend mask (time x lat x lon)
      Mask data
iii)  Get mean of masked/weighted data

iv)   Output timeseries as text
"""
# Set up lat/lon grid  -> ideally this would be automated from the input dimensions
lon = np.linspace(0,356.25,96)
lat = np.linspace(-90,90,73)
LON,LAT = np.meshgrid(lon,lat)
tdim = 1706

def weighting(data):
  """ Assume the earth is a perfect sphere (close enough!), 
  calculate the weighting for a given latitude and create a
  weighting array, then use that to weight the data """
  coslat = np.cos(LAT*np.pi/180)
  clext = np.repeat(np.expand_dims(coslat,axis=0),tdim,0)
  out_data = clext * data
  return out_data

def masklatlon(lonmin,lonmax,latmin,latmax):
  """ Create a mask (TRUE/FALSE vals) for a given latitude/longitude """
  mask=(LON>=lonmin)*(LON<=lonmax)*(LAT>=latmin)*(LAT<=latmax) # Use <,> or <=, >= ??
  maskext = np.repeat(np.expand_dims(mask,axis=0),tdim,0)
  return maskext

def maskdata(data,mask):
  """ Create a masked array using the weighted data
  reversing the mask (~) """
  out_data=np.ma.masked_array(data,~mask)
  #out_data=maskeddata.mask*maskeddata.data
  return out_data

def weightedave(data):
  """ Get the mean of the masked, weighted area. Easy! """
  ave = data.mean(axis=1).mean(axis=-1)
  return ave

def writetxt(ave_1D,area):
  outputstr = ['\n']*len(ave_1D)
  for i in xrange(len(ave_1D)):
    outputstr[i] = str(ave_1D[i]) + outputstr[i] #- Convert to string and add newlines
  fileout = open('nino3_had_' + area + '.txt', 'w') 
  fileout.writelines(outputstr) 
  fileout.close() 
  return













""" Shit I didn't use"""
## Program from Thom to find x,y for ANY given lat/lon
#def lonlatcalc(lon_in,lat_in):
#  dist=(LON-lon_in)**2+(LAT-lat_in)**2
#  jj,ii = np.where(dist==dist.min())
#  return ii, jj
## Import timeseries
#def readtext(txt_file):
#  # Function reads in correctly formatted text file of SST anomalies, outputs numpy array
#  fileobj = open(txt_file, 'r')
#  outputstr = fileobj.readlines()
#  fileobj.close()
#  textarray = np.array([[float(n) for n in line.split()] for line in outputstr])
#  return textarray
#
##function from Shannon
#from pyproj import Geod
#def get_dlon(lat):
#  g = Geod(ellps="WGS84")
#  return g.inv(1,lat,0,lat)[2]
#
#def get_dlat(lat):
#  g = Geod(ellps="WGS84")
#  return g.inv(0,lat,0,lat+1)[2]
#
#
## Regression
#"""Maps the regression values (or other values) between a timeseries, e.g. NINO3 and an array, e.g. Hadisst
#txt_ts is a text file of the timeseries, and col is the column to use (if there are more than one, otherwise leave blank)
#statvar is the statistical variable: slope=1, intercept=2, r_val=3, p_val=4, stderr=5
#"""
#def regmap(txt_ts,map_ts,statvar=0,col=0):
#  t,lat,lon = map_ts.shape
#  reg_map = np.zeros((lat,lon))
#  for i in xrange(lat):
#    for j in xrange(lon):
#      #slope, intercept, r_val, p_val, stderr = stats.linregress(map_ts[:,i,j],txt_ts[:,col])
#      stattup = stats.linregress(map_ts[:,i,j],txt_ts[:,col])
#      reg_map[i,j] = stattup[statvar]
#  out_reg = np.nan_to_num(reg_map)
#  return out_reg


## Export result map
#def writenc(nc_fileout,var):
#  nc_out = S.NetCDFFile(nc_fileout,mode='a')
#  var=nc_out.variables['temp']
#  data_out=var.getValue()
#  ancil2=np.expand_dims(ancil,axis=1)
#  data_out[:,:,:,:]=ancil2[:,:,:,:]
#  var[:]=data_out
#  nc_out.close()
#  return ancil2


## Display result map
#
#def makemap(nc_obj,mapvar,color_code,v_min=None,v_max=None):
#  """
#  Map a climate variable using basemap and matplotlib
#  for color_code, need to do plt.cm.hot_r etc
#  """
#  
#  if v_min == None:
#    v_min=mapvar.mean()-2*mapvar.std()
#    print v_min
#  elif v_max == None:
#    v_max=mapvar.mean()+2*mapvar.std()
#    print v_max
#
#  lon=nc_obj.variables['lon'].getValue()
#  lon_units = nc_obj.variables['lon'].units
#  lat = nc_obj.variables['lat'].getValue()
#  lat_units = nc_obj.variables['lat'].units
#  #- Make 2-D longitude and latitude arrays:
#  [lon2d, lat2d] = np.meshgrid(lon, lat)
#
#  #- Set up map:
#  plt.clf()
#  mapproj = bm.Basemap(projection='cyl', llcrnrlat=-90.0, llcrnrlon=0.0, urcrnrlat=90.0, urcrnrlon=360.0)
#  mapproj.drawcoastlines()
#  mapproj.drawparallels(np.array([-90, -45, 0, 45, 90]), labels=[1,0,0,0])
#  mapproj.drawmeridians(np.array([0, 90, 180, 270, 360]), labels=[0,0,0,1])
#  lonall, latall = mapproj(lon2d, lat2d)
#
#  themap = plt.pcolormesh(lonall, latall, mapvar, vmin=v_min, vmax=v_max, cmap=color_code)  #plt.cm.hot_r)
#  #themapc = plt.contour(lonall, latall, mapvar, 10, color=(0,0,0))
#  #plt.clabel(themapc, fontsize=12)
#  themap.cmap.set_over((0.0,0.0,0.0))
#  themap.cmap.set_under((0.9,0.9,0.9))
#  plt.title('DJF Mean Temperature response, [K]', fontsize=12)
#  plt.axis([0, 360, -90, 90])
#  #plt.xlabel('Longitude [' + lon_units + ']')
#  #plt.ylabel('Latitude [' + lat_units + ']')
#  plt.colorbar(themap, orientation='horizontal', extend='both', ticks=[-200,-100,-50,-20,-10,-5,-2,-1,0,1,2,5,10,20,50,100,200])
#  plt.savefig('temp.test.eps')
#  plt.draw()
#  #plt.show()
#  return

