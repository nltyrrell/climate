import Scientific.IO.NetCDF as S
import numpy as n
import matplotlib.pyplot as plt
import sys
import mpl_toolkits.basemap as bm

clim = S.NetCDFFile('./ncfiles/test.temp.nc', mode='r')
gpot = S.NetCDFFile('./ncfiles/test.gp.nc',mode='r')


# puts the data into a t x h x lat x lon array
gpot_data = gpot.variables['ht'].getValue()
clim_data = clim.variables['temp'].getValue()

# Define constants
R=287
g=9.81

# height in m, masked and temp too
hei = n.ma.masked_where(gpot_data<=0,gpot_data)
temp=n.ma.masked_where(clim_data<=250,clim_data)

# layer thicknes in km
dhei = n.diff(hei,axis=1)/1000

# diffence in temp
dtemp = n.diff(temp,axis=1)

# laps rate in K/1000m
laps = dtemp/dhei
lapsm= n.ma.masked_where(laps>5,laps)

Ltmean = n.mean(laps,axis=0)
ltm03 = n.mean(Ltmean[3:4,:,:],axis=0)
lon=clim.variables['longitude'].getValue()
lon_units = clim.variables['longitude'].units
lat = clim.variables['latitude'].getValue()
lat_units = clim.variables['latitude'].units
#[lonall, latall] = n.meshgrid(lon, lat)
#- Make 2-D longitude and latitude arrays:

[lon2d, lat2d] = n.meshgrid(lon, lat)


#- Set up map:

mapproj = bm.Basemap(projection='cyl', 
                     llcrnrlat=-90.0, llcrnrlon=0.0,
                     urcrnrlat=90.0, urcrnrlon=360.0)
mapproj.drawcoastlines()
mapproj.drawparallels(n.array([-90, -45, 0, 45, 90]), labels=[1,0,0,0])
mapproj.drawmeridians(n.array([0, 90, 180, 270, 360]), labels=[0,0,0,1])
lonall, latall = mapproj(lon2d, lat2d)

plt.clf()
mymapf = plt.contourf(lonall, latall, ltm03, n.linspace(-10,0,21), cmap=plt.cm.jet)
#mymap = plt.contour(lonall, latall, ltm03, n.linspace(-10,0,6), colors='k')
#mymapf.set_clim(0, 100)
#mymap.set_clim(0, 100)
#plt.clabel(mymap, fontsize=12)
plt.axis([0, 360, -90, 90])
#plt.xlabel('Longitude [' + lon_units + ']')
#plt.ylabel('Latitude [' + lat_units + ']')
plt.colorbar(mymapf, orientation='horizontal')
plt.savefig('exercise-T-contour.png')
#plt.show()

#raise SystemExit
sys.exit("done")
