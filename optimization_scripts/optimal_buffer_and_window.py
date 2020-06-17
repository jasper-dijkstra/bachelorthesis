# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 15:38:49 2020

@author: Jasper Dijkstra

"""

import os, sys
from datetime import datetime
import itertools
import time
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Local imports
sys.path.append("..") # get one directory up
import utilities as ut
import masking_functions as mask
import main as init


# ========================================================
# FUNCTIONS
# ========================================================
def GetTotalLandCells(boundaries, target_lat, target_lon):
    
    # Initialize array from data
    arr = np.ones((target_lat, target_lon))
    
    # Count amount of arrays that are land
    land = mask.land_sea_mask(arr, boundaries)
    total_land_cells = int(len(land[land != 0]))
    
    return total_land_cells


def GetStats(daily_data_dict):
    # Initialize lists to get the average of all data
    total_tropomi = []
    total_gfed_edgar = []
    ext_in_buffer = [] # external databases in buffers
    explained = []
    
    for day in daily_data_dict:
        plumes = np.copy(daily_data_dict[day]['plume_mask'].flatten())
        plumes = plumes[plumes > 0]
        frequency = np.bincount(plumes)
    
        # Make sure frequency has got enough indices to complete this function:
        append_values = 114 - len(frequency)
        frequency = np.lib.pad(frequency, ((0,append_values)), 'constant', constant_values=(0))

        # append total TROPOMI cells to corresponding list
        total = frequency[1] + frequency[12] + frequency[112] + frequency[113]
        total_tropomi.append(total)
        
        # append total GFED and EDGAR cells to corresponding list
        total_gfed_edgar.append(sum(frequency[2:]))
        
        # append total GFED and EDGAR cells within buffer to corresponding list
        total_ext_in_buffer = frequency[11] + frequency[12] + frequency[101] \
            + frequency[102] + frequency[111] + frequency[112]
        ext_in_buffer.append(total_ext_in_buffer)
        
        # Reclass to get other stats as well
        plumes[plumes == 12] = 11
        plumes[plumes == 102] = 101
        plumes[plumes == 112] = 111

        total_tropomi_in_buffer = frequency[1] + frequency[11] + frequency[111]
        total_tropomi_in_gfed = frequency[11]
        total_tropomi_in_edgar = frequency[101]
        total_tropomi_in_edgar_gfed = frequency[111]
        
        # Append total TROPOMI cells that can be explained with the GFED and EDGAR databases
        exp_by_gfed = round(100 * ((total_tropomi_in_gfed + total_tropomi_in_edgar_gfed) / total_tropomi_in_buffer), 1)
        exp_by_edgar = round(100 * ((total_tropomi_in_edgar + total_tropomi_in_edgar_gfed) / total_tropomi_in_buffer), 1)
        explained.append(round((((exp_by_gfed + exp_by_edgar)*total)/100),1))
    
    stats = [np.nanmean(np.array(total_tropomi).astype(int)), \
             np.nanmean(np.array(explained).astype(int)), \
                 np.nanmean(np.array(total_gfed_edgar).astype(int)), \
                     np.nanmean(np.array(ext_in_buffer).astype(int))]

    return stats


def Plot2D(X, Y, score, xlabel='', ylabel='', cblabel='', colormap='rainbow'):
    
    # Do not open figure
    plt.ioff()
    
    # Initialize Figure
    fig = plt.figure()
    ax = fig.gca()
    
    # Plot the surface.
    cs = ax.pcolormesh(X, Y, score, norm=colors.PowerNorm(gamma=0.5), cmap=colormap, linewidth=0, antialiased=False)

    # Customize axes
    # Set Labels, axticks and Distances
    #plt.xticks(X[0])
    ax.set_xlabel(xlabel)
    ax.xaxis.labelpad=10
    plt.yticks(Y[:,0])
    ax.set_ylabel(ylabel)
    ax.yaxis.labelpad=10
    
    # Set colorbar
    cb = plt.colorbar(cs, orientation = 'vertical', fraction=0.046, pad=0.04)
    cb.set_label(cblabel)
      
    # Save the image
    plt.savefig(os.path.join(os.getcwd() + rf'//2d_pcolormesh_{cblabel}.png'), bbox_inches=None, dpi=1200) #'tight'
    
    plt.cla()
    plt.clf()
    
    return

# ========================================================
# USER DEFINED PARAMETERS
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

# == Directories ==
# Main Directory
basepath = ut.DefineAndCreateDirectory(r'C:\Users\jaspd\Documents\Python\00_bachelorthesis\bachelorthesis\THESIS_WORKINGDIR\\')

# Directory where GFED files are stored
GFED_path = os.path.join(basepath, '01_GFED_hdf5_files' + os.path.sep) # Path to gfed .hdf5 files
EDGAR_path = os.path.join(basepath, '02_EDGAR_files' + os.path.sep) # Path to EDGAR .nc files

# Set modelparams
windowsizes = list(np.linspace(30, 200, 18).astype(int))
#stepsizes = list(np.linspace(20, 100, 9).astype(int))
buffersizes = list(np.linspace(2, 15, 14).astype(int))
#stdevs = list(np.round(np.linspace(0.2, 3, 15), 1))

# All loops that will have to be executed
iterations = list(itertools.product(buffersizes, windowsizes))

# Set variables
estimator = np.zeros((len(buffersizes), len(windowsizes))) # y-hat (total TROPOMI identifications)
estimate = np.zeros((len(buffersizes), len(windowsizes))) # y (explained/correct TROPOMI identifications)


# ========================================================
# RUN ALGORITHM
# ========================================================
# Notify the algorithm has started
starttime = datetime.now()
print(f'Algorithm started at: {starttime}')

# Model Behaviour settings
apply_land_sea_mask = True
apply_overlap_filter = False
use_wind_rotations = False

# Initialize other required parameters
boundaries, target_lon, target_lat, files = init.InitializingParameters(lat_min, lat_max, lon_min, lon_max, lonres, latres, basepath)

# Designated breaks (for progress updates)
breaks = start = int(len(iterations)/20)
breaklist = []
for i in range(20):
    breaklist.append(breaks)
    breaks += start

# Read/open the data
daily_data = init.CollectingData(boundaries, target_lon, target_lat, files, basepath, apply_land_sea_mask, use_wind_rotations)

# Iterate over all model parameter combinations
for count, i in enumerate(iterations):
    # params = [
        # buffersize (radius of buffer around TROPOMI plumes),
        # st.devs (minimum amount of st.devs within identification frames),
        # windowsize (size of moving window frame in grid cells),
        # stepsize (steps between each moving window frame in grid cells)
        # ]
    params = [i[0], 1, i[1], 20]
    
    # Detection algorithm
    daily_data = init.Detection(params, daily_data, boundaries, GFED_path, EDGAR_path, lonres, latres, apply_overlap_filter, use_wind_rotations)
    
    # Retrieve statistics from plume detection result
    stats = GetStats(daily_data)
    
    # Identify locations indices to place data in array
    x = np.where(buffersizes == iterations[count][0])[0]
    y = np.where(windowsizes == iterations[count][1])[0]
    
    # Append these statistics to array
    estimator[x,y] = stats[0]
    estimate[x,y] = stats[1]
    
    # Print progress update to console (steps of 5%)
    if count in breaklist:
        index = np.where(np.array(breaklist) == count)[0][0] + 1
        print(f'{5*index}% complete')
        #print('going to sleep for 5 minutes...')
        #time.sleep(300)

print('100% complete')
print(f'Total time elapsed: {datetime.now()-starttime}')

# Calculate Root mean Square Deviation and initialize X, Y coorindates
X, Y = np.meshgrid(windowsizes, buffersizes)

# Where is difference TROPOMI and EXPLAINED the smallest -> Root mean square deviation?
difference = estimator-estimate # Should negative values be filtered?
RMSD = np.sqrt(((difference)**2)/GetTotalLandCells(boundaries, target_lat, target_lon))
RMSD = RMSD / Y

# Tips from Sander (socre function)
y_tot = stats[2] # total sources according to edgar + gfed
y_pos = estimator # total identifications
y_neg = estimator-estimate # total non-explained identifications
#y_neg[y_neg <= 0] = 0 # Remove negative combinations -> more cells explained than identified?

score = (y_pos - y_neg) / Y 

# ============= CREATE OUTPUTS ==============
# Store data as csv's
np.savetxt(os.path.join(os.getcwd() + rf"\RMSD.csv"), RMSD, delimiter=",")
np.savetxt(os.path.join(os.getcwd() + rf"\score.csv"), score, delimiter=",")

# Plot a pcolormesh (2D) image of results, saved in scripts workingdir
xlabel = 'window length and width (x10 km)'
ylabel = 'radius buffer (x10 km)'
cblabel = 'RMSD'

Plot2D(X, Y, ma.array(RMSD, mask=(score == 0)), xlabel=xlabel, ylabel=ylabel, cblabel=cblabel)

# Now also plot score function
cblabel = 'score'
Plot2D(X, Y, ma.array(score, mask=(score == 0)), xlabel=xlabel, ylabel=ylabel, cblabel=cblabel, colormap='RdYlGn')
#Plot2D(X, Y, ma.array(score, mask=(y_neg == 0)), xlabel=xlabel, ylabel=ylabel, cblabel=cblabel, colormap='RdYlGn')
