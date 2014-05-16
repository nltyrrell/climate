import Scientific.IO.NetCDF as S
import numpy as np
import matplotlib.pyplot as plt
import sys
import mpl_toolkits.basemap as bm
import scipy as sc
from scipy import stats
from mycmaps import jetwhite
from matplotlib.colors import LogNorm

def worldshift(mapvar,color_code=jetwhite(n=256),vmin=None,vmax=None, name='test', title='title',lons='longitude',lats='latitude',clf=None,cbaror='Horizontal'):
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
	lon_ax=np.linspace(0,355,latlon[1])
	dataout, lon_axis = bm.shiftgrid(1.0, mapvar,lon_ax,start=True)
	lat_axis=np.linspace(-90,90,latlon[0])
	lon_units = 'degrees east'
	lat_units = 'degrees north'
	#- Make 2-D longitude and latitude arrays:
	[lon2d, lat2d] = np.meshgrid(lon_axis, lat_axis)

	#- Set up map:
	if clf==None:
		plt.clf()
	mapproj = bm.Basemap(projection='cyl', llcrnrlat=-90.0, llcrnrlon=0, urcrnrlat=90.0, urcrnrlon=360)
	mapproj.drawcoastlines()
	mapproj.drawparallels(np.array([-90, -45, 0, 45, 90]), labels=[1,0,0,0])
	mapproj.drawmeridians(np.array([0, 90, 180, 270, 360]), labels=[0,0,0,1])
	lonall, latall = mapproj(lon2d, lat2d)

	plt.ion()
	themap = plt.pcolormesh(lonall, latall, dataout, vmin=vmin, vmax=vmax, cmap=color_code)  #plt.cm.hot_r)
	#themap = plt.contourf(lonall, latall, mapvar, np.linspace(vmin,vmax,21), cmap=color_code, extend="both")
	#plt.clabel(themapc, fontsize=12)
	themap.cmap.set_over((0.0,0.0,0.0))
	themap.cmap.set_under((0.7,0.7,0.7))
	plt.title(title, fontsize=16)
	plt.axis([0, 355, -90, 90])
	#plt.xlabel('Longitude [' + lon_units + ']')
	#plt.ylabel('Latitude [' + lat_units + ']')
	if cbaror!=None:
		plt.colorbar(themap, orientation='horizontal', extend='both', ticks=np.linspace(vmin,vmax,11)) #orientation=cbaror,
	plt.savefig(name + '.eps')
	plt.show()
	return

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
	lon=np.linspace(-1.875,354.375,latlon[1])  #latlon[1])
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
	#themap = plt.pcolormesh(lonall, latall, mapvar, vmin=vmin, vmax=vmax, cmap=color_code)  #plt.cm.hot_r)
	themap = plt.contourf(lonall, latall, mapvar, np.linspace(vmin,vmax,21), cmap=color_code, extend="both")
	#plt.clabel(themapc, fontsize=12)
	themap.cmap.set_over((0.0,0.0,0.0))
	themap.cmap.set_under((0.7,0.7,0.7))
	plt.title(title, fontsize=16)
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
  mapproj.drawparallels(np.array([-90,-60,-30, 0,30, 60, 90]), labels=[1,0,0,0])
  mapproj.drawmeridians(np.array([0, 90, 180,270,360]), labels=[0,0,0,1])
  lonall, latall = mapproj(lon2d, lat2d)

  plt.ion()
  #themap = plt.pcolormesh(lonall, latall, mapvar, vmin=vmin, vmax=vmax, cmap=color_code)  #plt.cm.hot_r)
  themap = plt.contourf(lonall, latall, mapvar, np.linspace(vmin,vmax,21), cmap=color_code, extend="both")
  #plt.clabel(themapc, fontsize=12)
  themap.cmap.set_over((0.0,0.0,0.0))
  themap.cmap.set_under((0.7,0.7,0.7))
  plt.title(title, fontsize=14)
  plt.axis([0, 356, -90, 90])
  #plt.xlabel('Longitude [' + lon_units + ']')
  #plt.ylabel('Latitude [' + lat_units + ']')
  plt.colorbar(themap, orientation='horizontal', extend='both', ticks=np.linspace(vmin,vmax,11))
  plt.savefig(name + '.eps')
  plt.show()
  return
  #raise SystemExit

def hatchmap(nc_obj,mapvar,insig=np.array([1,2]),color_code=jetwhite(n=256),vmin=None,vmax=None, name='test', title='title',lons='longitude',lats='latitude',cbar=True):
	"""
	Map a climate variable using basemap and matplotlib
	for color_code, need to do plt.cm.hot_r etc
	Hatches out insignificant or other bits given by 'insig', which is p value masked at greater 0.05
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
	plt.figure(figsize=(9,4))
	mapproj = bm.Basemap(projection='cyl', llcrnrlat=-90.0, llcrnrlon=0.0, urcrnrlat=90.0, urcrnrlon=360.0)
	mapproj.drawcoastlines()
	mapproj.drawparallels(np.array([-90,-60,-30, 0,30, 60, 90]), labels=[1,0,0,0])
	mapproj.drawmeridians(np.array([0, 90, 180,270,360]), labels=[0,0,0,1])
	lonall, latall = mapproj(lon2d, lat2d)

	#themap = plt.pcolormesh(lonall, latall, mapvar, vmin=vmin, vmax=vmax, cmap=color_code)  #plt.cm.hot_r)
	themap = plt.contourf(lonall, latall, mapvar, np.linspace(vmin,vmax,21), cmap=color_code, extend="both")
	#plt.clabel(themapc, fontsize=12)
	themap.cmap.set_over((0.0,0.0,0.0))
	themap.cmap.set_under((0.7,0.7,0.7))
	plt.title(title, fontsize=19)
	plt.axis([0, 356, -90, 90])
	#plt.xlabel('Longitude [' + lon_units + ']')
	#plt.ylabel('Latitude [' + lat_units + ']')
	if cbar == True:
		plt.colorbar(themap, orientation='vertical', extend='both', ticks=np.linspace(vmin,vmax,11))
	if insig.shape == (73,96):
		themap = plt.contourf(lonall,latall, insig.mask.astype(int), hatches=['','','.',''],colors="none")
	plt.savefig(name + '.eps')
	plt.show()
	return
	#raise SystemExit

def linsimplot(var, clf=True):
	plt.ion()
	if clf==True:
		plt.clf()
	plt.plot(var)
	plt.grid()

def simplot(var,clf=True,clim=False):
	dims = len(var.shape)
	if dims == 4:
		var = var[0,0,::]
	elif dims == 3:
		var = var[0,::]
	else:
		var=var
	plt.ion()
	if clf==True:
		plt.clf()
	if clim!=False:
		plt.pcolormesh(var,cmap=jetwhite(n=256)).set_clim(vmin=-clim,vmax=clim)
	if clim==False:
		plt.pcolormesh(var)
	if clf==True:
		plt.colorbar()
	return

def hadmap(mapvar,color_code=jetwhite(n=256),vmin=None,vmax=None, name='test', title='title',lons='longitude',lats='latitude',clf=None,cbaror='Horizontal'):
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
	lon_ax=np.linspace(-180,180,latlon[1])
	dataout, lon_axis = bm.shiftgrid(0, mapvar,lon_ax,start=True)
	lat_axis=np.linspace(-90,90,latlon[0])
	lon_units = 'degrees east'
	lat_units = 'degrees north'
	#- Make 2-D longitude and latitude arrays:
	[lon2d, lat2d] = np.meshgrid(lon_axis, lat_axis)

	#- Set up map:
	if clf==None:
		plt.clf()
	mapproj = bm.Basemap(projection='cyl', llcrnrlat=-90.0, llcrnrlon=0, urcrnrlat=90.0, urcrnrlon=360)
	mapproj.drawcoastlines()
	mapproj.drawparallels(np.array([-90, -45, 0, 45, 90]), labels=[1,0,0,0])
	mapproj.drawmeridians(np.array([0, 90, 180, 270, 360]), labels=[0,0,0,1])
	lonall, latall = mapproj(lon2d, lat2d)

	plt.ion()
	#themap = plt.pcolormesh(lonall, latall, dataout, vmin=vmin, vmax=vmax, cmap=color_code)  #plt.cm.hot_r)
	themap = plt.contourf(lonall, latall, dataout, np.linspace(vmin,vmax,21), cmap=color_code, extend="both")
	#plt.clabel(themapc, fontsize=12)
	themap.cmap.set_over((0.0,0.0,0.0))
	themap.cmap.set_under((0.7,0.7,0.7))
	plt.title(title, fontsize=16)
	plt.axis([0, 360, -90, 90])
	#plt.xlabel('Longitude [' + lon_units + ']')
	#plt.ylabel('Latitude [' + lat_units + ']')
	if cbaror!=None:
		plt.colorbar(themap, orientation='horizontal', extend='both', ticks=np.linspace(vmin,vmax,11)) #orientation=cbaror,
	plt.savefig(name + '.eps')
	plt.show()
	return

#---------------------------------------- Log colorbar ------------------


def logmap(nc_obj,mapvar,color_code=jetwhite(n=256),vmin=None,vmax=None, name='test', title='title',lons='longitude',lats='latitude'):
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
  mapproj.drawparallels(np.array([-90,-60,-30, 0,30, 60, 90]), labels=[1,0,0,0])
  mapproj.drawmeridians(np.array([0, 90, 180,270,360]), labels=[0,0,0,1])
  lonall, latall = mapproj(lon2d, lat2d)

  plt.ion()
  #themap = plt.pcolormesh(lonall, latall, mapvar, vmin=vmin, vmax=vmax, cmap=color_code)  #plt.cm.hot_r)
  themap = plt.contourf(lonall, latall, mapvar, norm = LogNorm(vmin=vmin,vmax=vmax), cmap=color_code)
  #plt.clabel(themapc, fontsize=12)
  themap.cmap.set_over((0.0,0.0,0.0))
  themap.cmap.set_under((0.7,0.7,0.7))
  plt.title(title, fontsize=14)
  plt.axis([0, 356, -90, 90])
  #plt.xlabel('Longitude [' + lon_units + ']')
  #plt.ylabel('Latitude [' + lat_units + ']')
  plt.colorbar(themap, orientation='horizontal')
  plt.savefig(name + '.eps')
  plt.show()
  return
  #raise SystemExit

