# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 14:28:08 2020

@author: Jasper Dijkstra

Open the CAMS database, to get the atmospheric background

Source: Inness et al. (2019), http://www.atmos-chem-phys.net/19/3515/2019/

https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-reanalysis-eac4?tab=form

"""


from netCDF4 import Dataset
import netCDF4
import numpy as np
import os

# Local imports
import raster_tools as raster


def GetCorrectOutput(fid, day, month):
    """
    Temporary function to open CAMS netCDF file supplied by Sander
    """

    # initialize time coversion
    time = fid.variables['time']
    date = netCDF4.num2date(time[:], time.units, time.calendar)
    months = np.zeros(date.shape)
    days = np.zeros(date.shape)
    for i in range(len(date)):
        months[i] = date[i].month
        days[i] = date[i].day
    
    # Only get the required indices
    indices = np.where((days == day) & (months == month))
    CO_level = fid.variables['co'][indices[0],:,:,:]
    
    return CO_level



def CAMSWeights():
    """
    Get list with weights for levels in CAMS data, to compute weighted average
    """
    weights = [237.000,353.732,516.378,716.195,944.860,1194.27,1457.65,1727.73,1998.78,\
               2264.98,2521.33,2763.06,2986.61,3188.15,3365.06,3514.64,3635.34,3725.89,\
                   3785.58,3814.08,3811.81,3779.81,3719.42,3632.63,3522.12,3390.72,3242.46,\
                       3081.20,2911.81,2739.41,2570.01,2393.17,2210.68,2024.48,1836.72,\
                           1649.57,1466.02,1288.61,1111.31,946.417,763.870,616.533,497.615,\
                               401.634,324.166,261.640,211.174,170.442,137.567, 111.033,\
                                   90.6462,75.4761,63.7167,54.1947,46.1010,38.8463,31.9892,\
                                       25.2225,18.4253,20.0000]
    return weights


def OpenCAMS(path, bbox, xres, yres, day, month):
    
    # EDGAR longitudes range from 0-360, while TROPOMI does -180-180
    # Therefore apply correction to bbox    
    cams_bbox = [bbox[0], bbox[1], bbox[2] + 180, bbox[3] + 180]
    
    # Open Dataset
    fid = Dataset(path)
    
    # Retrieve latitude and longitudes
    lon = fid.variables['longitude'][:]
    lat = fid.variables['latitude'][:]
    
    # Get CO level, and take the mean in all dimensions, until lat/lon dimensions are left
    #CO_level = fid.variables['co'][:] # 4-dimensions: [hours, model_level, lat, lon]
    CO_level = GetCorrectOutput(fid, day, month)
    
    # Take mean of time and atmospheric height dimensions, so only lat/lon (2D) are left
    CO_level = np.mean(CO_level, axis=0)
    weights = np.array(CAMSWeights()).astype(int) # Getting weighted average
    CO_level = np.average(CO_level, weights=weights, axis=0) # CO level in kg/kg
    
    # Put longitudes in correct order: [180 - 360] block to [-180 - 0]
    west = np.squeeze(CO_level[:,np.where(lon > 180)])
    east = np.squeeze(CO_level[:,np.where(lon <= 180)])
    CO_level = np.concatenate((west,east),axis=1)
    
    # Close dataset
    fid.close()
    
    # Get indices of latitudes and logitudes where desired extent is True
    ilon = np.where((lon >= cams_bbox[2]) & (lon <= cams_bbox[3]))
    ilat = np.where((lat >= cams_bbox[0]) & (lat <= cams_bbox[1]))
    
    # Clip lon, lat and CO_emissions to desired extent
    lon = lon[ilon[0]]
    lat = lat[ilat[0]]
    CO_level = CO_level[ilat[0], :]
    CO_level = CO_level[:, ilon[0]]
    
    # Resample to desired extent
    CO_level = raster.ResampleArray(bbox = bbox, array = CO_level, lon_resolution = xres, lat_resolution = yres)
    
    # Convert from kg/kg to ppb (assuming 1 ppb = 1 nano mole CO / mole dry air)
    ppb_factor = (28.01 * 1e9) / 28.9674 # Molar mass CO / Molar mass dry air
    CO_level = CO_level * ppb_factor
    
    # Turn the data around
    CO_level = np.flipud(CO_level)
    
    return CO_level



def FetchCams(CAMS_path, dates, bbox, xres, yres):
    """
    Fetch CAMS data (CO background concentration), from downloaded files.
    
    In the future it will be possible to download the required data with an api
    
    For now, please download the files required.
    If file does not exist, algorithm will proceed without CAMS

    Parameters
    ----------
    CAMS_path : path to .nc file
    day : day of TROPOMI overpass
    month : month of TROPOMI overpass
    year : year of TROPOMI overpass
    bbox : boundaries [lat_min, lat_max, lon_min, lon_max]
    xres : resolution in lon direction
    yres : resolution in lat direction

    Returns
    -------
    CO_level : array
        2D ndarray with CO background concentration.

    """
    # Dates
    day = dates[0]
    month = dates[1]
    year = dates[2]
    
    # TODO: API retrieval
    print('NOTE: In the future, CAMS data will be downloaded automatically with the cds api!')
    
    # For now, just open file without API
    # target file
    #file = os.path.join(CAMS_path + os.path.sep + rf'EAC4_m{month}_d{day}_y{year}.nc')
    file = os.path.join(CAMS_path + os.path.sep + rf'download.nc')
    
    if not os.path.isfile(file):
        print('CAMS file was not found!')
        print('Proceeding algorithm without CAMS')
        return
    
    CO_level = OpenCAMS(file, bbox, xres, yres, day, month)
    
    return CO_level


# =============================================================================
## Required for API
#import certifi
#import urllib3
#import cdsapi
#import yaml
#
# # Disable warnings
# urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
# 
# def Config_credentials_copernicus(store):
#     
#     assert store in ['ads', 'cds'], 'Credentials not recognised!'
#     
#     parent_dir = os.path.dirname(os.getcwd())
#     with open(rf'{parent_dir}\credentials\{str(store)}\.cdsapirc', 'r') as f:
#         credentials = yaml.safe_load(f)
# 
#     return credentials
# 
# 
# 
# def DownloadCAMS(timescope, basepath):
#     """
#     This function downloads ERA5 data making use of the CDS API from Copernicus,
#     as netCDF4 files for the time- and aerial scope defined in daily_data_dict
#     
#     Parameters
#     ----------
#     zone_extents : output of GetTimeZoneExtent()
#     timescope : list, tuple
#         time scope in the form of: [day, month, year]
#     pressure_level : int
#         Atmospheric pressure level (hPa) at which data has to be downloaded
#     basepath : string
#         Path to folder where meteodata will be stored
#         
#         
#     Returns
#     -------
#     downloaded .nc file in: (<working_dir>/store_meteo_data/<subdir>).
# 
#     """
#      
#     # redefining earlier results
#     year = str(timescope[0])
#     month = str(timescope[1])
#     day = str(timescope[2])
#     
#     # Setting output directory
#     #curr_directory = os.getcwd()
#     out_dir = os.path.join(basepath, rf'02_Store_CAMS_data//')
#     ut.DefineAndCreateDirectory(out_dir)
#     filename = f'EAC4_m{month}_d{day}_y{year}_api.nc'
#     
#     if os.path.isfile(os.path.join(out_dir, filename)) == True:
#         print('CAMS data already exists!')
#         print('No new attempt to download the data will be undertaken!')
#         return os.path.join(out_dir, filename)
#     
#     # Request parameters
#     timelist = ['00:00', '03:00', '06:00', '09:00', '12:00', '15:00', '18:00', '21:00']
#     model_level_list = list(map(str, range(1,61)))
#     pressure_level_list = ['1', '2', '3', '5', '7', '10', '20', '30', '50', '70', \
#                            '100', '150', '200', '250', '300', '400', '500', '600', \
#                                '700', '800', '850', '900', '925', '950', '1000']
#     
#     # Actual request:
#     credentials = Config_credentials_copernicus(store = 'ads')
#     c = cdsapi.Client(url=credentials['url'], key=credentials['key'])
#     
# 
#     http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
#     http.request('GET', credentials['url'])
#     c.retrieve(
#     'cams-global-reanalysis-eac4',
#     {
#         'variable': 'carbon_monoxide',
#         'pressure_level': pressure_level_list,
#         'model_level': model_level_list,
#         'date': f'{year}-{month}-{day}/{year}-{month}-{str(int(day)+1)}',
#         'time': timelist,
#         'format': 'netcdf',
#     },
#     os.path.join(out_dir, filename))
#     
#     
#     return os.path.join(out_dir, filename)
# =============================================================================