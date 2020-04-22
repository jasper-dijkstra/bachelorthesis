# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 10:13:52 2020

@author: jaspd
"""

import h5py
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


    
#%%


hdffile = r'C:\Users\jaspd\Documents\Python\00_bachelorthesis\bachelorthesis\GFED_fire_emissions_mask\GFED4.1s_2018_Beta.hdf5'

saving_fig_dir = r'C:\Users\jaspd\Desktop\THESIS_WORKINGDIR\plume_figures\fig_GEFD_emissions_10_13_C_2018_max1.png'

"""
# Defining boundaries
lat_min = -50
lat_max = 0
lon_min = 100
lon_max = 160
"""

lon_min = -179.875
lon_max = 179.875
lat_min = -89.875
lat_max = 89.875

# define a box
box_lat = [lat_min, lat_max]
box_lon = [lon_min, lon_max]


#%%

f = h5py.File(hdffile, mode='r')

# this is the full path within the HDF file
name = '/emissions/10/C'
daily_frac = '/emissions/10/daily_fraction/day_13'
# get the geolocation data
latitude = f['/lat'][:]
longitude = f['/lon'][:]
#lat_index = np.logical_and(latitude > box_lat[0], latitude < box_lat[1])
#lon_index = np.logical_and(longitude > box_lon[0], longitude < box_lon[1])
#box_index = np.logical_and(lat_index, lon_index)
#data = f[name][box_index]
data = f[name][:]
daily = data = f[daily_frac][:]
data = np.multiply(data, daily)
#data[data > 1] = 1
data[data == 0] = np.NaN # Set all zero values to 'nan'

#data = np.random.random((720,1440))

#%%

# Deciding on the nlon_t and nlat_t
field_t = data
nlon_t = len(field_t[0])
nlat_t = len(field_t)

# make new masked array of valid averages
count_t = np.ones((len(field_t), len(field_t[0])))
mask = (count_t == 0)
field_mt = ma.array(field_t, mask=mask)


# plot with cartopy
fig = plt.figure(figsize=(10,6))
ax = plt.axes(projection=ccrs.PlateCarree())

# Deciding what data will be plotted

cs = plt.pcolormesh(longitude, latitude, field_mt, cmap='rainbow', transform=ccrs.PlateCarree())
cbaxes = fig.add_axes([0.2, 0.1, 0.6, 0.03]) 
cb = plt.colorbar(cs, cax = cbaxes, orientation = 'horizontal' )
cb.set_label('g C m^-2 month^-1')

ax.coastlines()
plt.ioff() # Preventing figures from appearing as pop-up

# Saving figure
fig.savefig(saving_fig_dir, bbox_inches='tight')
#print('figure saved at: {}'.format(saving_path))
plt.close()