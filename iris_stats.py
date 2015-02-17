# (C) British Crown Copyright 2013, Met Office
#
# This file is part of Iris.
#
# Iris is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Iris is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Iris.  If not, see <http://www.gnu.org/licenses/>.
"""
Statistical operations between cubes.

"""

import numpy as np

import iris
import scipy.stats as stats

def _get_calc_view(cube_a, cube_b, corr_coords):
    """
    This function takes two cubes and returns cubes which are
    flattened so that efficient comparisons can be performed
    between the two.

    Args:

    * cube_a:
        First cube of data
    * cube_b:
        Second cube of data, compatible with cube_a
    * corr_coords:
        Names of the dimension coordinates over which
        to calculate correlations.

    Returns:

    * reshaped_a/reshaped_b:
        The data arrays of cube_a/cube_b reshaped
        so that the dimensions to be compared are
        flattened into the 0th dimension of the
        return array and other dimensions are
        preserved.
    * res_ind:
        The indices of the dimensions that we
        are not comparing, in terms of cube_a/cube_b

    """

    # Following lists to be filled with:
    # indices of dimension we are not comparing
    res_ind = []
    # indices of dimensions we are comparing
    slice_ind = []
    for i, c in enumerate(cube_a.dim_coords):
        if not c.name() in corr_coords:
            res_ind.append(i)
        else:
            slice_ind.append(i)

    # sanitise input
    dim_coord_names = [c.name() for c in cube_a.dim_coords]
    if corr_coords is None:
        corr_coords = dim_coord_names

    if ([c.name() for c in cube_a.dim_coords] !=
            [c.name() for c in cube_b.dim_coords]):
        raise ValueError("Cubes are incompatible.")

    for c in corr_coords:
        if c not in dim_coord_names:
            raise ValueError("%s coord "
                             "does not exist in cube." % c)

    # Reshape data to be data to correlate in 0th dim and
    # other grid points in 1st dim.
    # Transpose to group the correlation data dims before the
    # grid point dims.
    data_a = cube_a.data.view()
    data_b = cube_b.data.view()
    dim_i_len = np.prod(np.array(cube_a.shape)[slice_ind])
    dim_j_len = np.prod(np.array(cube_a.shape)[res_ind])
    reshaped_a = data_a.transpose(slice_ind+res_ind)\
                       .reshape(dim_i_len, dim_j_len)
    reshaped_b = data_b.transpose(slice_ind+res_ind)\
                       .reshape(dim_i_len, dim_j_len)

    return reshaped_a, reshaped_b, res_ind


def pearsonr(cube_a, cube_b, corr_coords=None):
    """
    Calculates the n-D Pearson's r correlation
    cube over the dimensions associated with the
    given coordinates.

    Returns a cube of the correlation between the two
    cubes along the dimensions of the given
    coordinates, at each point in the remaining
    dimensions of the cubes.

    For example providing two time/altitude/latitude/longitude
    cubes and corr_coords of 'latitude' and 'longitude' will result
    in a time/altitude cube describing the latitude/longitude
    (i.e. pattern) correlation at each time/altitude point.

    Args:

    * cube_a, cube_b (cubes):
        Between which the correlation field will be calculated.
        Cubes should be the same shape and have the
        same dimension coordinates.
    * corr_coords (list of str):
        The cube coordinate names over which to calculate
        correlations. If no names are provided then
        correlation will be calculated over all cube
        dimensions.

    Returns:
        Cube of correlations.

    Reference:
        http://www.statsoft.com/textbook/glosp.html#Pearson%20Correlation

    """

    # If no coords passed then set to all coords of cube.
    if corr_coords is None:
        corr_coords = [c.name() for c in cube_a.dim_coords]

    vec_a, vec_b, res_ind = _get_calc_view(cube_a,
                                           cube_b,
                                           corr_coords)

    sa = vec_a - np.mean(vec_a, 0)
    sb = vec_b - np.mean(vec_b, 0)
    flat_corrs = np.sum((sa*sb), 0)/np.sqrt(np.sum(sa**2, 0)*np.sum(sb**2, 0))

    corrs = flat_corrs.reshape([cube_a.shape[i] for i in res_ind])

    # Construct cube to hold correlation results.
    corrs_cube = iris.cube.Cube(corrs)
    corrs_cube.long_name = "Pearson's r"
    corrs_cube.units = "1"
    for i, dim in enumerate(res_ind):
        c = cube_a.dim_coords[dim]
        corrs_cube.add_dim_coord(c, i)

    return corrs_cube, flat_corrs, sa, sb

def linregress(cube_a, cube_b, corr_coords=None):
    """
    NOT WORKING
    might not be possible, or easy, without a 
    loop for calc linear regression
    16/01/2014

    Calculates linear regression between two cubes
    
    Started with 'pearsonr' and used scipy.stats
    to calc regressions instead

    Returns a cube of the linear regression between the two
    cubes along the dimensions of the given
    coordinates, at each point in the remaining
    dimensions of the cubes.

    For example providing two time/altitude/latitude/longitude
    cubes and corr_coords of 'latitude' and 'longitude' will result
    in a time/altitude cube describing the latitude/longitude
    (i.e. pattern) correlation at each time/altitude point.

    Args:

    * cube_a, cube_b (cubes):
        Between which the correlation field will be calculated.
        Cubes should be the same shape and have the
        same dimension coordinates.
    * corr_coords (list of str):
        The cube coordinate names over which to calculate
        correlations. If no names are provided then
        correlation will be calculated over all cube
        dimensions.

    Returns:
        Cube of regression coefficients, flat_corrs, sa, sb

    """

    # If no coords passed then set to all coords of cube.
    if corr_coords is None:
        corr_coords = [c.name() for c in cube_a.dim_coords]

    vec_a, vec_b, res_ind = _get_calc_view(cube_a,
                                           cube_b,
                                           corr_coords)

    sa = vec_a - np.mean(vec_a, 0)
    sb = vec_b - np.mean(vec_b, 0)
    n = len(sa[:,0])
    TINY = 1.0e-20
    xmean = np.mean(sa,None)
    ymean = np.mean(sb,None)
    ssxm, ssxym, ssyxm, ssym = np.cov(sa[:,0], sb[:,0], bias=1).flat
    r_num = ssxym
    r_den = np.sqrt(ssxm*ssym)
    if r_den == 0.0:
        r = 0.0
    else:
        r = r_num / r_den
        # test for numerical error propagation
        if (r > 1.0):
            r = 1.0
        elif (r < -1.0):
            r = -1.0
    # z = 0.5*log((1.0+r+TINY)/(1.0-r+TINY))
    df = n-2
    t = r*np.sqrt(df/((1.0-r+TINY)*(1.0+r+TINY)))
#     prob = distributions.t.sf(np.abs(t),df)*2
    slope = r_num / ssxm
    intercept = ymean - slope*xmean
    sterrest = np.sqrt((1-r*r)*ssym / ssxm / df)
    return slope, intercept, r, sterrest

#     # following line was from pearsonr function
#     #     flat_corrs = np.sum((sa*sb), 0)/np.sqrt(np.sum(sa**2, 0)*np.sum(sb**2, 0))
# 
#     corrs = flat_corrs.reshape([cube_a.shape[i] for i in res_ind])
# 
#     # Construct cube to hold correlation results.
#     corrs_cube = iris.cube.Cube(corrs)
#     corrs_cube.long_name = "Pearson's r"
#     corrs_cube.units = "1"
#     for i, dim in enumerate(res_ind):
#         c = cube_a.dim_coords[dim]
#         corrs_cube.add_dim_coord(c, i)
# 
#     return corrs_cube, flat_corrs, sa, sb


