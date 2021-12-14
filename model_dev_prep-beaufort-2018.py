#!/usr/bin/env python
# coding: utf-8


import numpy as np
import os
import pandas as pd
import xarray as xr

import matplotlib as mpl
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
# from matplotlib.colors import DivergingNorm
# import matplotlib.patches as patches
#%matplotlib inline
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import seaborn as sns; sns.set(color_codes=False)
import cmocean
# from windrose import WindroseAxes
import warnings

import cftime
import datetime


plt.rcParams['xtick.bottom'] = True # keep my tick marks
plt.rcParams['ytick.left'] = True
plt.rcParams['font.size'] = 18
# plt.rcParams["font.family"] = "Arial"
# plt.rcParams['figure.figsize'] = 12,8
# mpl.rcParams['figure.dpi'] = 300 # activate for presentation quality

# from sklearn.metrics.pairwise import haversine_distances

## this is a dummy grid that has the right conventions
grid = xr.open_dataset(
    '/glade/work/vcooper/grid_ref/sithick_SImon_CESM2_piControl_r1i1p1f1_gn_110001-120012.nc')
latmin = 72
latmax = 79
lonmin = 195
lonmax = 230
beau_mask = (
    (grid.coords['lon'] > lonmin) 
    & (grid.coords['lon'] < lonmax)
    & (grid.coords['lat'] > latmin)
    & (grid.coords['lat'] < latmax)
)

## circle boundary for plotting
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpl.path.Path(verts * radius + center)

## keep cice 2018 output because of detailed grid
mdir_adj = '/glade/u/home/vcooper/work/BGEP_vtc/adjfiles/' # location of adjusted files, changed dims to lat and lon
cice18 = xr.open_dataset(mdir_adj + 'cicefsdww3i.cice.h1.0086.nc') # 2018 rerun
cice18.TLAT[:] = grid.lat # hack to fix polar stereo
cice18.TLON[:] = grid.lon
tarea = cice18.tarea.rename(
    {'TLON':'longitude','TLAT':'latitude'}).reset_coords(names=['ULON','ULAT'],drop=True)

import copy
import dask
from glob import glob
latslice=slice(300,384)

#%%time
###### Coupled Model Update May 2021 #####
## loading in and concat data 
ww2012 = xr.open_dataset('/glade/scratch/vcooper/waveice_analysis/cesm23iws1tsks.ww3.hi.2012.nc')
# ww2013 = xr.open_dataset('/glade/scratch/vcooper/waveice_analysis/cesm23iws1tsks.ww3.hi.2013.nc')
# ww2014 = xr.open_dataset('/glade/scratch/vcooper/waveice_analysis/cesm23iws1tsks.ww3.hi.2014.nc')
# ww2015 = xr.open_dataset('/glade/scratch/vcooper/waveice_analysis/cesm23iws1tsks.ww3.hi.2015.nc')
# ww2016 = xr.open_dataset('/glade/scratch/vcooper/waveice_analysis/cesm23iws1tsks.ww3.hi.2016.nc')
# ww2017 = xr.open_dataset('/glade/scratch/vcooper/waveice_analysis/cesm23iws1tsks.ww3.hi.2017.nc')
# ww2018 = xr.open_dataset('/glade/scratch/vcooper/waveice_analysis/cesm23iws1tsks.ww3.hi.2018.nc')
# ww2019 = xr.open_dataset('/glade/scratch/vcooper/waveice_analysis/cesm23iws1tsks.ww3.hi.2019.nc')

ww_dict = {'ww2012' : ww2012}
#            'ww2013' : ww2013,
#            'ww2014' : ww2014,
#            'ww2015' : ww2015,
#            'ww2016' : ww2016,
#            'ww2017' : ww2017,
#            'ww2018' : ww2018,
#            'ww2019' : ww2019}

temp_f = xr.open_dataset('./ww1719ef_beau_cat.nc')

for key,val in ww_dict.items():
    print(key)
#     val = val.rename({'UAX': 'uwnd',
#                 'UAY': 'vwnd',
#                 'ICE': 'ice',
#                 'HS':  'hs',
#                 'T02': 't02',
#                 'T0M1':'t0m1',
#                 'T01': 't01',
#                 'FP0': 'fp',
#                 'THM': 'dir',
#                 'EF':  'ef',
#                 'FREQ':'f',
#                 'NX':  'ni',
#                 'NY':  'nj'})
    
    val['latitude'] = (['nj','ni'],grid.lat.values)
    val['longitude'] = (['nj','ni'],grid.lon.values)
    val = val.set_coords(['time','latitude','longitude'])
#     val['f'] = temp_f.f
    val['FREQ'] = temp_f.f.values
#     val = val.sel(nj=latslice)
#     val['tarea'] = (['nj','ni'],cice18.tarea.sel(nj=latslice).values)
    val['tarea'] = (['nj','ni'],cice18.tarea.values)
    val.coords['mask'] = (('nj','ni'), beau_mask.values)
#     val = val.where(val.mask > 0, drop=True) ## only keep if dropping to beau
    ww_dict[key] = val

##########################################

from sklearn.metrics.pairwise import haversine_distances

## Updated version to 15% ice concentration threshold
def icedistance(iceconc_input):
    # turn icefracs to numpy array
    icefracsnp = iceconc_input.values
    lats = iceconc_input.TLAT.values # cice version
    lons = iceconc_input.TLON.values # cice version
#     lats = iceconc_input.latitude.values # wavewatch version
#     lons = iceconc_input.longitude.values # wavewatch version


    # create array to hold the distances
    distances = icefracsnp.copy() # same size array as the evaluated data
    distances -= distances # make zeros or nan; we will keep these values for cells that don't need a calc

    
    ##### GET OPEN WATER -> WATER/ICE EDGE LOCATIONS #####
    
    # get all open water locations except at edge of domain to avoid computation breaking
    icefracsnp_noborder = icefracsnp[1:-1,1:-1] # exclude borders for open water checking neighbors
    locations_openw = np.transpose(np.where(icefracsnp_noborder<0.15))
    locations_openw += 1 # adjust indices for the border exclusion

    # create 4 arrays, each represents the offset of open water location in coords by 1 unit
    latp1 = np.append(locations_openw[:,0]+1,locations_openw[:,1]).reshape(locations_openw.shape,order='F')
    latm1 = np.append(locations_openw[:,0]-1,locations_openw[:,1]).reshape(locations_openw.shape,order='F')
    lonp1 = np.append(locations_openw[:,0],locations_openw[:,1]+1).reshape(locations_openw.shape,order='F')
    lonm1 = np.append(locations_openw[:,0],locations_openw[:,1]-1).reshape(locations_openw.shape,order='F')

    # get max icefrac of 4 neighbor cells at each open water cell
    iceneighbormax = np.nanmax(np.stack((icefracsnp[latp1[:,0],latp1[:,1]],
                                         icefracsnp[lonm1[:,0],lonm1[:,1]],
                                         icefracsnp[lonp1[:,0],lonp1[:,1]],
                                         icefracsnp[latm1[:,0],latm1[:,1]])),axis=0)

    # get index of the open water cells with ice neighbor>15% # these are values for which we will calc distance
    wateredge = locations_openw[np.where(iceneighbormax>0.15)]
    wateredgeT = wateredge.T
    wateredgelatlon = np.array([[lats[wateredgeT[0],wateredgeT[1]]],
                                [lons[wateredgeT[0],wateredgeT[1]]]]).squeeze().T # Nx2 matrix of lat,lon
    
    ##### CALCULATION OF DISTANCES #####
    
    ### addition to limit to beaufort
    latmin = 72
    latmax = 79
    lonmin = 195
    lonmax = 230
    #################################
    
    # get all cell locations with ice > 15% 
    # and in beaufort region
    #icewhere = np.where(icefracsnp>0.15)
    icewhere = np.where((icefracsnp>0.15)
                        & (lons > lonmin) 
                        & (lons < lonmax)
                        & (lats > latmin)
                        & (lats < latmax))
    
    icecells = np.transpose(icewhere) # index by array position
    icelatlon = np.array([[lats[icewhere]],
                          [lons[icewhere]]]).squeeze().T # Nx2 matrix of lat,lon

    # calculate minimum distance
    mindist = haversine_distances(np.deg2rad(icelatlon),
                                  np.deg2rad(wateredgelatlon)).min(axis=1)*6371000/1000 # x by Radius-earth for km
    
    icecellsT = np.transpose(icecells) # transpose for vectorized indexing
    distances[icecellsT[0],icecellsT[1]] = mindist # put mindist into each grid point

    return(distances)


#%%time
## BOOKMARK
## calculate distance in ice for beaufort region

## load one time step of just the ice concentrations
## SET THIS
casename = 'cesm23iw_dtice100'
ctrlname = 'cesm23iws1tsks_bittest1'
monsel = np.arange(1,3)

###### WW3 ########
def preprocess(ds):
    ds = ds[['ICE']]
    ds = ds.expand_dims(time = [datetime.datetime.now()]) ## dummy time
    return ds

## experiment ww3
rundir = '/glade/scratch/vcooper/' + casename + '/run/'
rundir_cp = copy.copy(rundir)

year = str(2018)

flist = casename + '.ww3.hi.'+year+'-*.nc'

dfiles = sorted(glob(rundir + flist))
N = len(dfiles)
dfiles = dfiles[0:N]
dfiles_n = len(dfiles)

dfiles_np = np.array(dfiles)

print(len(dfiles), len(ww2012.time[0:N]))

mfds_temp = xr.open_mfdataset(dfiles,
                              combine='by_coords',
                              coords='minimal',
                              concat_dim='time',
                              compat='override',
                              preprocess=preprocess,
                              parallel=False)
mfds_temp['time'] = ww2012.time[0:N]
exp_ww = mfds_temp.load()

## add coordinates
# exp_ww = exp_ww.assign_coords({'FREQ':temp_f.f.values})
exp_ww['TLON'] = (exp_ww.ICE[0].dims, grid.lon.values)
exp_ww['TLAT'] = (exp_ww.ICE[0].dims, grid.lat.values)

exp_ww = exp_ww.set_coords(['TLAT','TLON'])

## initialize distances array
ice_conc_data = exp_ww.ICE
distances = np.zeros(ice_conc_data.shape)

with warnings.catch_warnings():
    warnings.simplefilter("ignore") ## ignore warnings from div by zero leading to nans

    ## loop through time, calculating distance from ice edge
    for i,da in enumerate(ice_conc_data):
        distances[i] = icedistance(da)

## convert to xarray
distances_da = xr.DataArray(distances, dims=ice_conc_data.dims,coords=ice_conc_data.coords).rename('dist')

## create beaufort subset of distances

latmin = 72
latmax = 79
lonmin = 195
lonmax = 230
beau_mask_xr = xr.where((distances_da.TLON > lonmin) 
                      & (distances_da.TLON < lonmax)
                      & (distances_da.TLAT > latmin)
                      & (distances_da.TLAT < latmax),1,0)

distances_b = distances_da.where(beau_mask_xr,drop=True)


## load other beaufort data variables

## BOOKMARK
## calculate distance in ice for beaufort region

## load one time step of just the ice concentrations
## SET THIS
casename = 'cesm23iw_dtice100'
ctrlname = 'cesm23iws1tsks_bittest1'

###### WW3 ########
def preprocess(ds):
    ds = ds[['HS','ICE','EF']].where(beau_mask_xr,drop=True)
    ds = ds.expand_dims(time = [datetime.datetime.now()]) ## dummy time
    return ds

## experiment ww3
rundir = '/glade/scratch/vcooper/' + casename + '/run/'
rundir_cp = copy.copy(rundir)
flist = casename + '.ww3.hi.'+year+'-*.nc'

dfiles = sorted(glob(rundir + flist))[0:N]
dfiles_n = len(dfiles)

dfiles_np = np.array(dfiles)

print(len(dfiles), len(ww2012.time[0:N]))

mfds_temp = xr.open_mfdataset(dfiles,
                              combine='by_coords',
                              coords='minimal',
                              concat_dim='time',
                              compat='override',
                              preprocess=preprocess,
                              parallel=False)
mfds_temp['time'] = ww2012.time[0:N]
exp_ww = mfds_temp.load()

## add coordinates
exp_ww = exp_ww.assign_coords({'FREQ':temp_f.f.values})
exp_ww['TLON'] = (distances_b[0].dims, distances_b.TLON.values)
exp_ww['TLAT'] = (distances_b[0].dims, distances_b.TLAT.values)

exp_ww = exp_ww.set_coords(['TLAT','TLON'])
exp_ww['dist'] = distances_b

new_filename = '/glade/scratch/vcooper/waveice_analysis/model_dev/'+casename+'.ww3.hi.beaufort.'+year+'.nc'
print('saving to ', new_filename)

exp_ww.to_netcdf(path=new_filename)
print('finished saving')
