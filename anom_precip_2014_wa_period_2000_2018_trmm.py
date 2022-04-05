#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 08:47:30 2019

@author: sowm
"""





import os
import salem
import matplotlib
import numpy as np
import xarray as xr
import netCDF4 as nc4
import numpy.ma as ma
import geopandas as gpd
import matplotlib as mpl
import cartopy.crs as ccrs
from pandas import read_csv
import matplotlib.pyplot as plt
import cartopy.mpl.ticker as cmt
import matplotlib.colors as colors
import scipy.integrate as integrate
import matplotlib.gridspec as gridspec
from datetime import datetime, timedelta
import cartopy.io.shapereader as shpreader
from netCDF4 import Dataset, num2date, date2num
from mpl_toolkits.axes_grid1 import make_axes_locatable


plt.rc('text', usetex=True)
try:
    from calc_dryspells import *
except:
    print("Module 'spells' not imported.")
    
try:
    from mask_array_using_shapefile import *
except:
    print("Module 'mask' not imported.")

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))  
#=======================================================
coord = read_csv('/home/sowm/magatte/These/data/coord_stations.csv')
coord= coord.set_index('station')
#- Reading data for a timeslice, latitude, and longitude:
fname = '/home/sowm/Africa.shp'
fig_dir = '/home/sowm/magatte/These/calculs/python/spells/figures/cartes/'
filename = 'trmm_3b42_dly_anom_2014_2000_2018_wa.nc'
FILEPATH = '/home/sowm/magatte/These/data/data_cnrm/ppt_trmm_3b42/west_africa_extract/'
inputnc = FILEPATH+filename

#os.system('cdo -setmissval,-99 -ifthen mask_senegal_chirps.nc '+inputnc+ ' chirps_dly_2018_18w11w_11n18n_masked.nc')



proj = ccrs.PlateCarree()
fig,axes = plt.subplots(2, 3, subplot_kw={'projection': proj}, figsize=(17,10), gridspec_kw={'wspace':0.15,'hspace':-0.3, 'top':0.95, 'bottom':0.05, 'left':0.05, 'right':0.92}, facecolor='none')
ax=axes.flat
#=======================================================


#fig.subplots_adjust(hspace=0.01, wspace=0.01)

levels = np.linspace(-5, 5, 17)

#=======================================================
mon_len = [31,30,31,31,30,31]
months = ['May', 'Jun', 'July', 'August', 'September', 'October']
d1 = int(date2num(datetime(2014, 5, 1, 0, 0), units='days since 2014-01-01', calendar='gregorian'))   # 2014 MJJASO
d2 = int(date2num(datetime(2014, 10, 31, 0, 0), units='days since 2014-01-01', calendar='gregorian')) # period
d_1 = int(date2num(datetime(2000, 1, 1, 0, 0), units='days since 1981-01-01', calendar='gregorian'))   # 2000-2015
d_2 = int(date2num(datetime(2015, 12, 31, 0, 0), units='days since 1981-01-01', calendar='gregorian')) # period

nc = nc4.Dataset(FILEPATH + filename, mode='r')
dly_anom = nc['pr'][:]
dly_anom_mai = ma.mean(dly_anom[120:150,:,:], axis=0)
dly_anom_jun = ma.mean(dly_anom[151:180,:,:], axis=0)
dly_anom_jul = ma.mean(dly_anom[181:211,:,:], axis=0)
dly_anom_aug = ma.mean(dly_anom[212:242,:,:], axis=0)
dly_anom_sep = ma.mean(dly_anom[243:273,:,:], axis=0)
dly_anom_oct = ma.mean(dly_anom[274:303,:,:], axis=0)
data = [dly_anom_mai, dly_anom_jun, dly_anom_jul, dly_anom_aug, dly_anom_sep, dly_anom_oct]
# =============================================================================
# nc = nc4.Dataset(FILEPATH + 'chirps_dly_2014_18w11w_11n18n_masked.nc', mode='r')
# mjjaso_pr = nc['precip'][d1:d2]
# mjjaso_pr_cum  = np.sum(mjjaso_pr, axis=0)
# mjjaso_pr_mean = np.mean(mjjaso_pr, axis=0)
# nc1= nc4.Dataset(inputnc, mode='r')
# an_pr = nc1['precip'][d_1:d_2]  # 2000-2015 extraction
# an_pr_yrs = []
# for yr in range(2000,2016):
#     ds = int(date2num(datetime(yr, 5, 1, 0, 0), units='days since 2000-01-01', calendar='gregorian'))   # yr MJJASO
#     de = int(date2num(datetime(yr, 10, 31, 0, 0), units='days since 2000-01-01', calendar='gregorian')) # period
#     an_pr_yr = an_pr[ds:de,:,:]
#     an_pr_yrs.append(an_pr_yr)
# an_pr_yrs_clim = ma.mean(an_pr_yrs, axis=0)
# anom_mjjaso = mjjaso_pr - an_pr_yrs_clim
# =============================================================================
lon = nc['longitude'][:]
lat = nc['latitude'][:] 
#nc.close()
#- Make 2-D longitude and latitude arrays:
[lon2d, lat2d] = np.meshgrid(lon, lat)
adm1_shapes = list(shpreader.Reader(fname).geometries())

cmap = mpl.cm.bwr_r
# =============================================================================
for i in range(6):
    
    ax[i].coastlines(resolution='10m')
    ax[i].add_geometries(adm1_shapes, ccrs.PlateCarree(), edgecolor='black', facecolor='none', linewidth=2)
    ax[i].set_global()
    ax[i].set_xticks(np.round(np.linspace(-20., 20., 11),1), crs=ccrs.PlateCarree())
    ax[i].tick_params(labelsize=10)
    ax[i].set_yticks(np.round(np.linspace(0., 25., 10),1), crs=ccrs.PlateCarree())
    lon_formatter = cmt.LongitudeFormatter(zero_direction_label=True)
    lat_formatter = cmt.LatitudeFormatter()
    ax[i].xaxis.set_major_formatter(lon_formatter)
    ax[i].yaxis.set_major_formatter(lat_formatter)
    ax[i].set_extent([-20., 20., 0., 25.], ccrs.PlateCarree())
    ax[i].text(0., 26., months[i], fontweight='bold', fontsize=14)
    
    #============   stations repesentation
    #for j in range(len(coord.index)):
        #ax[j].plot(coord['lon'][j], coord['lat'][j], 'ko',markersize=5) 
        #plt.legend(loc='best', ncol=2)
        #ax[j].text(coord['lon'][j],coord['lat'][j]+0.1, list(coord.index)[j][:3], fontsize=9, ha="center", va="center")  # :1 fo selecting a sub set of string in list (e.g K)
    
    #=====================================
    
    
    #mymapf = ax[i].contourf( lon2d, lat2d, mjjaso_pr_cum, levels, transform=ccrs.PlateCarree(), cmap=cmap)
    mymapf = ax[i].contourf( lon2d, lat2d, data[i], levels, transform=ccrs.PlateCarree(), cmap=cmap, extend = 'both')#,
#divider = make_axes_locatable(ax)
#ax_cb = divider.new_horizontal(size="5%", pad=0.08, axes_class=plt.Axes)
#fig.add_axes(ax_cb)
#plt.colorbar(mymapf, cax=ax_cb)
#plt.savefig(fig_dir+'Percentage_of_dry_days_from_'+mon+'_1981-2010_period_in_CHIRPS_over_Senegal.eps')
#fig.colorbar(mymapf, ax=axes.ravel().tolist())
#plt.rcParams['axes.labelweight'] = 'bold'    
cbar_ax = fig.add_axes([0.94, 0.2, 0.02, 0.6])
cbar = fig.colorbar(mymapf, cax=cbar_ax, pad=0.05)
cbar.set_label('(mm)', labelpad=-36, y=1.06, rotation=0)
for lb in cbar.ax.yaxis.get_ticklabels():  #  change ticklabels fontsize
    lb.set_weight("bold")
    lb.set_fontsize(12)   

########================== colorbar parameters  =========================================================================
    
fig.get_axes()[0].annotate('Precipitation anomaly for MJJASO 2014 over 2000-2018 period in TRMM over west africa', (0.5, 0.88),xycoords='figure fraction', ha='center',fontsize=14, fontweight='bold')
plt.savefig(fig_dir+'Precipitation_anomaly_for_mjjaso_2014_over_2000_2018_period_TRMM_over_west_africa', dpi = 300)
plt.show()



