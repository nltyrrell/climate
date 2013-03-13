import Scientific.IO.NetCDF as S
import numpy as np
import matplotlib.pyplot as plt
import sys
import mpl_toolkits.basemap as bm
import scipy as sc
from scipy import stats
from mycmaps import jetwhite

def worldmap(mapvar,color_code=jetwhite(n=256),vmin=None,vmax=None, name='test', title='title',lons='longitude',lats='latitude',clf=None,cbaror='Horizontal'):
	"""
	Map a climate variable using basemap and matplotlib
	for color_code: plt.cm.hot_r etc
	"""

	latlon=mapvar.shape

	if vmin == None:
		vmin=np.float('%8.0e'%(np.around(mapvar.mean(),decimals=0)))-np.float('%8.0e'%(2*mapvar.std()))
		print('vmin = '+str(vmin))
		print('mean = '+str(np.around(mapvar.mean(),decimals=0)))
	if vmax == None:
		vmax=np.float('%8.0e'%(np.around(mapvar.mean(),decimals=0)))+np.float('%8.0e'%(2*mapvar.std()))
		print('vmax = '+str(vmax))

#	mapvarII, lon=bm.shiftgrid(0., mapvar, np.linspace(0,356.25,latlon[1]), start=True)
	lon=np.linspace(-1.875,354.375,96)  #latlon[1])
	lat=np.linspace(-91,90,latlon[0])
	lon_units = 'degrees east'
	lat_units = 'degrees north'
	#- Make 2-D longitude and latitude arrays:
	[lon2d, lat2d] = np.meshgrid(lon, lat)

	#- Set up map:
	if clf==None:
		plt.clf()
	mapproj = bm.Basemap(projection='cyl', llcrnrlat=-90.0, llcrnrlon=-1.875, urcrnrlat=90.0, urcrnrlon=354.375)
	mapproj.drawcoastlines()
	mapproj.drawparallels(np.array([-90, -45, 0, 45, 90]), labels=[1,0,0,0])
	mapproj.drawmeridians(np.array([0, 90, 180, 270, 360]), labels=[0,0,0,1])
	lonall, latall = mapproj(lon2d, lat2d)

	plt.ion()
	themap = plt.pcolormesh(lonall, latall, mapvar, vmin=vmin, vmax=vmax, cmap=color_code)  #plt.cm.hot_r)
	#themap = plt.contourf(lonall, latall, mapvar, np.linspace(vmin,vmax,21), cmap=color_code, extend="both")
	#plt.clabel(themapc, fontsize=12)
	themap.cmap.set_over((0.0,0.0,0.0))
	themap.cmap.set_under((0.7,0.7,0.7))
	plt.title(title, fontsize=14)
	plt.axis([-2, 355, -90, 90])
	#plt.xlabel('Longitude [' + lon_units + ']')
	#plt.ylabel('Latitude [' + lat_units + ']')
	if cbaror!=None:
		plt.colorbar(themap, orientation='horizontal', extend='both', ticks=np.linspace(vmin,vmax,11)) #orientation=cbaror,
	plt.savefig(name + '.eps')
	plt.show()
	return

def makemap(nc_obj,mapvar,color_code=jetwhite(n=256),vmin=None,vmax=None, name='test', title='title',lons='longitude',lats='latitude'):
  """
  Map a climate variable using basemap and matplotlib
  for color_code, need to do plt.cm.hot_r etc
  """
  
  if vmin == None:
    vmin=np.floor(mapvar.mean()-2*mapvar.std())
    print vmin
  elif vmax == None:
    vmax=np.floor(mapvar.mean()+2*mapvar.std())
    print vmax

  lon=nc_obj.variables[lons].getValue()
  lon_units = nc_obj.variables[lons].units
  lat = nc_obj.variables[lats].getValue()
  lat_units = nc_obj.variables[lats].units
  #- Make 2-D longitude and latitude arrays:
  [lon2d, lat2d] = np.meshgrid(lon, lat)

  #- Set up map:
  plt.clf()
  mapproj = bm.Basemap(projection='cyl', llcrnrlat=-90.0, llcrnrlon=0.0, urcrnrlat=90.0, urcrnrlon=360.0)
  mapproj.drawcoastlines()
  mapproj.drawparallels(np.array([-90, -45, 0, 45, 90]), labels=[1,0,0,0])
  mapproj.drawmeridians(np.array([0, 90, 180, 270, 360]), labels=[0,0,0,1])
  lonall, latall = mapproj(lon2d, lat2d)

  plt.ion()
  themap = plt.pcolormesh(lonall, latall, mapvar, vmin=vmin, vmax=vmax, cmap=color_code)  #plt.cm.hot_r)
  #themap = plt.contourf(lonall, latall, mapvar, np.linspace(vmin,vmax,21), cmap=color_code, extend="both")
  #plt.clabel(themapc, fontsize=12)
  themap.cmap.set_over((0.0,0.0,0.0))
  themap.cmap.set_under((0.7,0.7,0.7))
  plt.title(title, fontsize=14)
  plt.axis([0, 357, -90, 90])
  #plt.xlabel('Longitude [' + lon_units + ']')
  #plt.ylabel('Latitude [' + lat_units + ']')
  plt.colorbar(themap, orientation='horizontal', extend='both', ticks=np.linspace(vmin,vmax,11))
  plt.savefig(name + '.eps')
  plt.show()
  return
  #raise SystemExit

def linsimplot(var):
	plt.ion()
	plt.clf()
	plt.plot(var)
	plt.grid()

def simplot(var,clf=True):
	plt.ion()
	if clf==True:
		plt.clf()
	plt.pcolormesh(var)
	plt.colorbar()
	return
