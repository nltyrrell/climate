import cdms2
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid

#read in data
fin = cdms2.open('~/project/access/clim/ncfiles/temp.sfc.clim.nc')
wsp = fin('temp')

# make cylindrical basemap 
m = Basemap(llcrnrlon=-180, llcrnrlat=-90,  
            urcrnrlon=180, urcrnrlat=90,
            projection='cyl')
#            resolution='c', projection='cyl')

# shift the longitude grid to be -180 to 180
datout, lon_axis = shiftgrid(185., wsp[0,0, ::]-wsp[1,0,::], wsp.getLongitude()[:], start=False)
#datout = wsp[12,0,:,:]-wsp[13,0,:,:]
lat_axis = wsp.getLatitude()[:]
#lon_axis = wsp.getLongitude()[:]

# create figure, add axes
fig1 = plt.figure(figsize=(8,10))
ax = fig1.add_axes([0.1,0.1,0.8,0.8])

# make 2-d grid of lons, lats
# compute native x,y coordinates of grid
lons, lats = numpy.meshgrid(lon_axis, lat_axis)
x, y = m(lons, lats)

# set desired contour levels.
clevs = numpy.arange(-10, 10, 2)
parallels = numpy.arange(-90., 90, 30.)
meridians = numpy.arange(0., 360., 30.)

# create the contour plot
cs1 = m.pcolormesh(x, y, datout)
#cs1 = m.contourf(x, y, datout, clevs)

# draw coastlines, parallels, meridians.
m.drawcoastlines(linewidth=1.5)
labels = [0, 1, 0, 1] #designate sides for labels
m.drawparallels(parallels, labels=labels, fontsize=8, linewidth=0.5)
m.drawmeridians(meridians, labels=labels, fontsize=8, linewidth=0.5)

# add colorbar
cb = m.colorbar(cs1, "bottom", pad="10%")
cb.set_label('m s-1')

# add title
ax.set_title('250 hPa wind speed, January 1979')

plt.show()
