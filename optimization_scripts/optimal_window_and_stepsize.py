# -*- coding: utf-8 -*-
"""
Created on Tue Jun 2 10:00:00 2020

@author: Jasper Dijkstra

Script to obtain and visualise the optimal window and step size for moving window,
for a TROPOMI resolution of 10 by 10 kilometers, between October and November 10, 2018

Script derived from main script.

"""

#--------------------
# IMPORTING FUNCTIONS
#--------------------
import os
import sys
import copy
from datetime import datetime
import numpy as np
import numpy.ma as ma
import time
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# Local imports
sys.path.append("..") # get one directory up
import handling_input as inpt
import handling_GFED as GFED
import handling_EDGAR as EDGAR
import handling_output_textfiles as txt
import handling_output_figures as figs
import masking_functions as mask
import raster_tools as raster
import utilities as ut

# ========================================================
# FUNCTIONS
# ========================================================

def GetTotalLandCells():
    lat_min = -50
    lat_max = 0
    lon_min = 100
    lon_max = 160
    boundaries = [lat_min, lat_max, lon_min, lon_max]
    
    # Data resolution
    latres = 10 # km
    lonres = 10 # km
    
    # Boundaries and resolution in array resolution
    target_lon = int(abs((lon_max-lon_min)/(lonres/110)))
    target_lat = int(abs((lat_max-lat_min)/(latres/110)))
    
    # Initialize array from data
    arr = np.ones((target_lat, target_lon))
    
    # Count amount of arrays that are land
    land = mask.land_sea_mask(arr, boundaries)
    total_land_cells = int(len(land[land != 0]))
    
    return total_land_cells


def Plot3D(X, Y, RMSE):
    
    # Do not open figure
    plt.ioff()
    
    # Initialize Figure
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    # Plot the surface.
    surf = ax.plot_surface(X, Y, RMSE, cmap='rainbow', linewidth=0, antialiased=False)
    
    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.01f'))
    
    ax.xaxis.grid(True, zorder=0)
    ax.yaxis.grid(True, zorder=0)
    
    # Set Labels and Distances
    ax.set_xlabel('window step size')
    ax.xaxis.labelpad=10
    ax.set_ylabel('window length and width')
    ax.yaxis.labelpad=10
    ax.set_zlabel('RMSE')
    ax.zaxis.labelpad=10
    
    # Set image view (elevation and azimuth)
    ax.view_init(elev=47, azim=134)
      
    # Save the image
    plt.savefig(os.path.join(os.getcwd() + rf'//3d_rmse_elev47_azim134.png'), bbox_inches=None, dpi=1200) #'tight'
    
    plt.cla()
    plt.clf()
    
    return


def Plot2D(X, Y, RMSE):
    
    # Do not open figure
    plt.ioff()
    
    # Initialize Figure
    fig = plt.figure()
    ax = fig.gca()
    
    # Plot the surface.
    cs = ax.pcolormesh(X, Y, RMSE, norm=colors.PowerNorm(gamma=0.5), cmap='rainbow', linewidth=0, antialiased=False)
    
    # Customize axes
    #ax.xaxis.grid(True, zorder=0)
    #ax.yaxis.grid(True, zorder=0)
    
    # Set Labels, axticks and Distances
    #plt.xticks(X[0])
    ax.set_xlabel('window step size')
    ax.xaxis.labelpad=10
    #plt.yticks(Y[:,0])
    ax.set_ylabel('window length and width')
    ax.yaxis.labelpad=10
    
    # Set colorbar
    cb = plt.colorbar(cs, orientation = 'vertical', fraction=0.046, pad=0.04)
    cb.set_label('RMSE')
      
    # Save the image
    plt.savefig(os.path.join(os.getcwd() + rf'//2d_rmse_pcolormesh.png'), bbox_inches=None, dpi=1200) #'tight'
    
    plt.cla()
    plt.clf()
    
    return


# ========================================================
# PARAMETERS
# ========================================================

# == Boundary Conditions and resolutions == 
# Set the Target Boundaries (degrees)
lon_min = 100 
lon_max = 160 
lat_min = -50 
lat_max = 0 

# Setting the approximate target resolution
lonres = 10 # km
latres = 10 # km

# == Desired Behaviour ==
# Apply operations:
apply_land_sea_mask = True  # filter out all TROPOMI data above the ocean
apply_overlap_filter = False # filter out plumes that occur for more than 3 consecutive satellite overpasses
    

# == Directories ==
# Main Directory
basepath = ut.DefineAndCreateDirectory(r'C:\Users\jaspd\Documents\Python\00_bachelorthesis\bachelorthesis\THESIS_WORKINGDIR\\')

# Directory where GFED files are stored
GFED_path = os.path.join(basepath, '01_GFED_hdf5_files' + os.path.sep) # Path to gfed .hdf5 files
EDGAR_path = os.path.join(basepath, '02_EDGAR_files' + os.path.sep) # Path to EDGAR .nc files

#%%
# ========================================================
# INITIALIZATION
# ========================================================
# Store all defined extents in a list:
boundaries = [lat_min, lat_max, lon_min, lon_max]

# Calculating resolution based on lonres, latres
target_lon = int(abs((lon_max-lon_min)/(lonres/110)))
target_lat = int(abs((lat_max-lat_min)/(latres/110)))

# Output directories, for coordinates and figures
fig_dir = ut.DefineAndCreateDirectory(os.path.join(basepath + r'\04_output\plume_figures'))

# Decide upon the window sizes and step sizes to be used
windowsizes = list(np.linspace(10, 200, 20).astype(int))
stepsizes = list(np.linspace(20, 100, 9).astype(int))

# Initiate 2D arrays to append the results to
yt = np.zeros((len(windowsizes), len(stepsizes))) # Predicted values (total TROPOMI)
y = np.zeros((len(windowsizes), len(stepsizes))) # Resulted values (overlap GFED/EDGAR and TROPOMI)
#n = np.zeros((len(windowsizes), len(stepsizes))) # Total amount of values (TROPOMI)

# Create a list with all files to apply the analysis on
input_files_directory = os.path.join(basepath + r'00_daily_csv\\')
files = ut.ListFilesInDirectory(input_files_directory, maxfiles=None)
#del files[0:3]

# Setting the time of starting the script
start = datetime.now()

# Reading daily csv files for specified area and day as np.arrays
daily_data = {}
for i, file in enumerate(files):    
    # Reading daily csv's as input array
    daily_data[i] = inpt.CSVtoArray(file, boundaries, target_lon, target_lat)
    
    # Filter measurements taken above the oceans (higher uncertainty)
    if apply_land_sea_mask == True:
        daily_data[i]['CO_ppb'] = mask.land_sea_mask(daily_data[i]['CO_ppb'], boundaries)
        daily_data[i]['count_t'] = mask.land_sea_mask(daily_data[i]['count_t'], boundaries)

print('Total time elapsed reading data: {}'.format(datetime.now()-start))



#%%
# ========================================================
# ALGORITHM
# ========================================================
start_main = datetime.now()

for wi, w in enumerate(windowsizes): # Loop over all window sizes
    for si, s in enumerate(stepsizes): # Loop over all step sizes
        """ Check for plumes based on three criteria: """
        #pred_total_gefdedgar = []
        obs_total_explained = []
        obs_total_tropomi = []
        if s < w: # Of course, the steps of the moving frame can never be larger than the window size!
            for day in daily_data:
                """ 1: At least 2 standard deviations above average of moving window """
                # Apply moving Window operation, with copy of CO_ppb
                arr = copy.deepcopy(daily_data[day]['CO_ppb'])
                plumes, co_average = raster.MovingWindow(arr, mask.identify_enhancements, \
                                                         window = (w,w), step = s) 
    
                # Check if plume was detected in at least 95% of windows
                plumes[plumes >= 0.95] = 1 # If true, plume (1)
                plumes[plumes < 0.95] = 0 # If false, no plume (0)
                
                # Removing nuisances by making sure there is at least one neighbour
                neighbors = raster.CountNeighbors(plumes)  # Identify neighbors of each grid cell
                plumes[neighbors <= 1] = 0 # If there are 1 or fewer neighbors, undo identification as plume
                
            
                #for day in daily_data:
                """ 3: Check if plumes overlap with modelled GFED and EDGAR data, by drawing a buffer around TROPOMI data """
                
                # Read GFED and EDGAR data
                gfed_array = GFED.OpenGFED(GFED_path, boundaries, daily_data[day]['day'], \
                                           daily_data[day]['month'], daily_data[day]['year'], lonres, latres)
                edgar_file = 'v432_CO_2012.0.1x0.1.nc' # EDGAR indsutry filename 
                edgar_array = EDGAR.OpenEDGAR(os.path.join(EDGAR_path + edgar_file), boundaries, lonres, latres)
                
# =============================================================================
#                 # Total cells of gfed & edgar
#                 gfededgar = gfed_array + edgar_array
#                 gfededgar = len(gfededgar[gfededgar != 0])
# =============================================================================
                
                gfed_array[gfed_array > 0] = 10     # GFED plumes = 10
                edgar_array[edgar_array > 0] = 100  # EDGAR plumes = 100
                
                # Draw Buffer around TROPOMI plumes
                plumes_buffered = raster.DrawCircularBuffer(plumes, radius = 7) # Draw circular buffers around TROPOMI plumes
            
                # Delete GFED and EDGAR values outside of plume buffer
                gfed_tropomi = (plumes_buffered * gfed_array) # GFED within TROPOMI
                gfed_tropomi[gfed_tropomi > 0] = 1
                edgar_tropomi = (plumes_buffered * edgar_array) # EDGAR within TROPOMI
                edgar_tropomi[edgar_tropomi > 0] = 1
                
                # Adding all arrays with different plume origins:
                plumes = (plumes + gfed_array + edgar_array + gfed_tropomi + edgar_tropomi).astype(int)
                
                # Now: (0: no plume, 1: TROPOMI plume, 10: GFED plume (within bufferzone of TROPOMI plume, \
                   # 11: TROPOMI + GFED identified plume))
                daily_data[day]['plume_mask'] = plumes
                
                copdict = copy.deepcopy(daily_data[day])
                stats = txt.GetStats(copdict)
                
                # Append data from specific day to list
                #pred_total_gefdedgar.append(stats[0])#gfededgar)
                obs_total_explained.append(((100-stats[3])*stats[0])/100)
                obs_total_tropomi.append(stats[0])
                
                # Generate output figures
                figs.PlotFigures(daily_data[day], os.path.join(basepath + rf'\\w{w}\\s{s}'), subplots=True)
            
            # Add the mean of listed values from timescope to earlier defined arrays
            yt[wi, si] = np.nanmean(np.array(obs_total_tropomi))
            y[wi, si] = np.nanmean(obs_total_explained)
            #n[wi, si] = np.nanmean(np.array(obs_total_tropomi))
                
        else:
            pass
          
    print(f'all steps for window window ({w}x{w}) completed, going to sleep now')
    
    time.sleep(300) # To cool processor, rest the script for 5 minutes...
    
    print(f'Restarting ({datetime.now()})')
    
print('Total time elapsed executing plume masking algorithm: {}'.format(datetime.now()-start_main))

#%%

# Calculate Root mean Square Error and initialize X, Y coorindates
X, Y = np.meshgrid(stepsizes, windowsizes)
RMSE = np.sqrt(((yt-y)**2)/GetTotalLandCells())

# Plot a pcolormesh (2D) and 3D image of results, saved in scripts workingdir
Plot2D(X, Y, ma.array(RMSE, mask=(RMSE == 0)))
RMSE[RMSE == 0] = np.nanmean(RMSE) # Set all RMSE values that are equal to zero to the max, as nans are not supported and 0 will distort the image
Plot3D(X, Y, RMSE)

print('Total time elapsed: {}'.format(datetime.now()-start))

