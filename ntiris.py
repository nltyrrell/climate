import iris
import numpy as np
import iris.coord_categorisation
import iris.plot as iplt
import iris.quickplot as qplt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys



def bc(aux_cube, dim_cube, comparable_coord):
    if type(aux_cube.coord(comparable_coord)) is not iris.coords.AuxCoord:
        raise TypeError
    if type(dim_cube.coord(comparable_coord)) is not iris.coords.AuxCoord:
        raise TypeError
    
    aux_cube_coord = aux_cube.coord(comparable_coord)
    dim_cube_coord = dim_cube.coord(comparable_coord)
    aux_cube_dim, = aux_cube.coord_dims(aux_cube_coord)
    dim_cube_dim, = aux_cube.coord_dims(dim_cube_coord)
    
    s_aux_cube = [slice(None)]*len(aux_cube.shape)
    s_aux_cube[aux_cube_dim] = 0
    s_dim_cube = [slice(None)]*len(dim_cube.shape)
    s_dim_cube[dim_cube_dim] = 0
    a = aux_cube[tuple(s_aux_cube)]
    a.attributes = None
    a.cell_methods = None
    b = dim_cube[tuple(s_dim_cube)]
    b.attributes = None
    b.cell_methods = None
    if not a.is_compatible(b):
        iris.util.describe_diff(a, b)
        raise RuntimeError("Cubes are not compatible")
    
    ind = []
    for p in aux_cube.coord(comparable_coord).points:
        i = np.where(dim_cube.coord(comparable_coord).points == p)
        ind.append(i[0][0])

    s = [slice(None)]*len(dim_cube.shape)
    s[dim_cube_dim] = ind
    new_data = dim_cube.data[tuple(s)]
    new_cube = aux_cube.copy()
    new_cube.data = new_data
    new_cube.history = "%s comparable to %s in terms of %s" % (dim_cube.name(),
                                                               aux_cube.name(),
                                                               comparable_coord)
    
    return new_cube

def remove_seascyc(cube, time_name='time'):
    """
    Remove seasonal cycle from montly timeseries
    Input: cube, time_name='time'
    Output: cube_rsc
    """
    try:
        iris.coord_categorisation.add_month_number(cube, time_name, 'month_number')
    except:
        pass
    else:
        print "month_number added as coord"
    cube_mean = cube[:,:,::].collapsed(time_name,iris.analysis.MEAN)
    cube_anom = cube-cube_mean
    cube_mon_mean = cube_anom.aggregated_by('month_number', iris.analysis.MEAN)
    seasonal_cycle = bc(cube_anom, cube_mon_mean, 'month_number')
    cube_rsc = cube_anom - seasonal_cycle
    return cube_rsc

def enscyc_ag(cube):
    ens = np.tile(np.linspace(1,48,48),24)
    #trim cube
    cube = cube[0:ens.shape[0]]
    cube.coord('month_number').long_name = '48_months'
    cube.coord('48_months').points = ens

    m48 = cube.aggregated_by('48_months',iris.analysis.MEAN)
    return m48
    

def composite_m48(cube_name, ncfile_path='/home/nicholat/project/pacemaker/ncfiles/',notanom=False):
    cube = iris.load_cube(ncfile_path+cube_name)
    try:
        cube.coord('t').standard_name='time'
    except:
        pass
    else:
        print "t coord changed to time"

    if notanom:
        cube_rsc = cube
        iris.coord_categorisation.add_month_number(cube_rsc, 'time', 'month_number')
    else:
        cube_rsc = remove_seascyc(cube) 
    cube_m48 = enscyc_ag(cube_rsc)
    cube_m48.long_name = cube.long_name
    if notanom:
        new_name = cube_name[:-2]+'m48.abs.nc'
    else:
        new_name = cube_name[:-2]+'m48.nc'
    iris.save(cube_m48,ncfile_path+new_name)
    return cube_m48, cube_rsc, cube

# u_cube = composite_m48('u.thlev.4ysl.fix.nc')
# rhum_cube, rhum_rsc, cube = composite_m48('rhum.plv.4ysl.nc',notanom=True)
# 
# sys.exit('all done')
# temp_cube, temp_rsc, cube = composite_m48('temp.plv.4ysl.nc',notanom=True)
# gpht_cube = composite_m48('gpht.plv.4ysl.nc')
# rh_cube = composite_m48('rhum.plv.4ysl.nc')
# p_cube = composite_m48('pres.sfc.4ysl.nc')
# v_cube = composite_m48('v.plev.4ysl.nc')
# u_cube = composite_m48('u.plev.4ysl.nc')
# w_cube = composite_m48('w.thlev.4ysl.fix.nc')
# composite_m48('lwflux.clsky.sfc.4ysl.nc')
# composite_m48('lhf.sfc.4ysl.nc')
# composite_m48('dlwr.sfc.4ysl.nc')
# composite_m48('dswr.sfc.4ysl.nc')


def regmean(cube,loni,lonf,lati,latf):
    """ Define a region and get the area weighted mean
    Input:  cube, lon_i, lon_f, lat_i, lat_f
    Output: cube_reg, cube_regmean
    """
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()
    region = iris.Constraint(longitude=lambda l: (loni <= l <= lonf), latitude = lambda l: (lati <= l <= latf))
    cube_reg = cube.extract(region)
    grid_areas = iris.analysis.cartography.area_weights(cube_reg)
    cube_regmean = cube_reg.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=grid_areas)
#     print 'mean sfc '+str(cube_regmean[0].data)
#     print 'mean p=10 '+str(cube_regmean[10].data)
    return cube_reg, cube_regmean


def nino3_plot(cube):
    """
    Plots the nino 3 timeseries for cube, surface temp wold make sense
    Input: Iris cube (sfc temp)
    Output: plot of nino3 timeseries, and cube of same thing
    """
    try:
        cube.coord('t').standard_name='time'
    except:
        pass
    else:
        print "t coord changed to time"

    cube_rsc = remove_seascyc(cube) 

    loni = 210; lonf = 270; lati = -5; latf = 5
# 			latlon = [-5,5,210,270] #for data from 0 - 360 degress
    nino3, nino3_mean = regmean(cube_rsc,loni,lonf,lati,latf)
    plt.ion()
    plt.clf()
    qplt.plot(nino3_mean[:,0])
    plt.title('NINO3 timeseries')

    return nino3, nino3_mean


def linregts(cube1,cube2,name1='name1',name2='name2',ncfile_path='/home/nicholat/project/pacemaker/ncfiles/'):
    """
    Calculate the regression between a 2d map and a 1d timeseries
    i.e. regress v onto T_sfc_Aus
    Input: cube_map, cube_ts (time coord only), name1, name2
    Output: reg_cube, cor_cube
    """

    linreg_map = np.zeros(cube1.shape[-2::])
    cor_map = np.zeros(cube1.shape[-2::])
    #     cor_map = np.zeros(copy_cube.shape)
    #     print('Linreg/Cor map for '+name1+' and '+name2)

    regress_ts = cube2.data
    for nlat, lat in enumerate(cube1.coord('latitude')):
        for nlon, lon in enumerate(cube1.coord('longitude')):
            # get the sfc temp timeseries at each lat, lon.
            var1 = cube1.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data
    #             var2 = cube2.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data

            linreg = stats.linregress(var1,regress_ts)
            linreg_map[nlat,nlon] = linreg[0]
            cor_map[nlat,nlon] = linreg[2]

        

    if cube1.ndim==4:
        reg_cube = cube1[0,0,::].copy()
    elif cube1.ndim==3:
        reg_cube = cube1[0,::].copy()
    reg_cube.data[:] = linreg_map
    reg_cube.long_name = 'Lin Regression '+ name1 +' '+ name2
    reg_cube.units = 'no_unit'
    reg_cube.attributes['title'] = 'Lin Regression '+ name1 +' '+ name2
    reg_cube.attributes['name'] = 'reg'
    try:
        reg_cube.remove_coord('surface')
    except:
        pass
    try:
        reg_cube.remove_coord('time')
    except:
        pass
    iris.save(reg_cube,ncfile_path+'lreg.4ysl.'+name1+'.'+name2+'.nc')
    if cube1.ndim==4:
        cor_cube = cube1[0,0,::].copy()
    elif cube1.ndim==3:
        cor_cube = cube1[0,::].copy()
    cor_cube.data[:] = cor_map
    cor_cube.long_name = 'Correlation '+ name1 +' '+ name2
    cor_cube.units = 'no_unit'
    cor_cube.attributes['title'] = 'Correlation '+ name1 +' '+ name2
    cor_cube.attributes['name'] = 'r_val'
    try:
        cor_cube.remove_coord('surface')
    except:
        pass
    try:
        cor_cube.remove_coord('time')
    except:
        pass
    iris.save(cor_cube,ncfile_path+'cor.4ysl.'+name1+'.'+name2+'.nc')

    return reg_cube, cor_cube


