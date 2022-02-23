#!/usr/bin/env python
# coding: utf-8

import numpy as np
import os
import pandas as pd
import xarray as xr

import cftime
import datetime



## this is a dummy grid that has the right conventions
grid = xr.open_dataset(
    '/glade/work/vcooper/grid_ref/sithick_SImon_CESM2_piControl_r1i1p1f1_gn_110001-120012.nc').isel(time=0)
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
## keep cice 2018 output because of detailed grid
mdir_adj = '/glade/u/home/vcooper/work/BGEP_vtc/adjfiles/' # location of adjusted files, changed dims to lat and lon
cice18 = xr.open_dataset(mdir_adj + 'cicefsdww3i.cice.h1.0086.nc') # 2018 rerun
cice18.TLAT[:] = grid.lat # hack to fix polar stereo
cice18.TLON[:] = grid.lon
tarea = cice18.tarea.rename(
    {'TLON':'longitude','TLAT':'latitude'}).reset_coords(names=['ULON','ULAT'],drop=True)

latslice=slice(300,384)

ww2018 = xr.open_dataset(
    '~/uwas0070/vcooper/waveice/scratch_analysis/cesm23iws1tsks.ww3.hi.2018.nc').ICE

temp_f = xr.open_dataset(
    '/glade/u/home/vcooper/analysis/waveice/ww1719ef_beau_cat.nc').f

scratchpath = '/glade/scratch/vcooper/'

## change these
case = 'ic4m3rad1'
ncasefile = '/cesm23iw_ic4m3'
wsel = 'w5nlDlnU'
wexp = '/windinpt_exp'

#fsel = scratchpath + 'cesm23iw_ic4m1/run/hourly_' + wsel + '/cesm23iw_ic4m1.ww3.ti.2018.EF.nc'
#fsel = scratchpath + 'cesm23iw_dtice100_w-nl-hf/run/hourly_t2t4-' + wsel + '/cesm23iw_dtice100_w-nl-hf.ww3.ti.2018.EF.nc'
fsel = scratchpath + 'cesm23iw_ic4m3/run/hourly_' + wsel + '/cesm23iw_ic4m3.ww3.ti.2018.EF.nc'
# fsel = scratchpath + 'cesm23iw_ic4m4/run/hourly_' + wsel + '/cesm23iw_ic4m4.ww3.ti.2018.EF.nc'
# fsel = scratchpath + 'cesm23iw_ic4m5/run/hourly_' + wsel + '/cesm23iw_ic4m5.ww3.ti.2018.EF.nc'
# fsel = scratchpath + 'cesm23iw_ic4m7/run/hourly_' + wsel + '/cesm23iw_ic4m7.ww3.ti.2018.EF.nc'
# fsel = scratchpath + 'cesm23iw_fsd-t2t4/run/hourly_' + wsel + '/cesm23iw_fsd-t2t4.ww3.ti.2018.EF.nc'

tempww = xr.open_dataset(fsel).isel(time=slice(0,8760))
tempww['time'] = ('time',ww2018.time)
tempww = tempww.set_coords('time').sel(time=slice('2018-07-01','2018-07-31'))
tempww['f'] = ('FREQ', temp_f)
tempww = tempww.set_coords('f')

npath = '/glade/campaign/univ/uwas0070/vcooper/waveice/zenodo_modeloutput/minfiles/'

newf = npath + case + wexp + ncasefile + '.ww3.ti.2018-07.EF.nc'
print('saving')
tempww.to_netcdf(newf)
print('finished ' + newf)

