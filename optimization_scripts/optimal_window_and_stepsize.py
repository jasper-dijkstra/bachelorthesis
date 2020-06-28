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
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Local imports
sys.path.append("..") # get one directory up
import utilities as ut
#import masking_functions as mask
import main as init


# ========================================================
# FUNCTIONS
# ========================================================



def GetStats(daily_data, labeled_plumes, num_features):
    
    # Make a copy of the identified plumes
    plumes = np.copy(daily_data['plumes_explained']) # Get the plumes
    
    # Initiate list with plume types
    plumetypes = []
    
    for i in range(1,num_features+1): # Loop over all labelled plumes
        plume = np.copy(labeled_plumes)
        plume[plume != i] = 0 # Keep only the labelled plume
        
        # Retrieve plume value
        plume_y, plume_x = np.where(plume != 0)
        plumevalue = plumes[plume_y[0], plume_x[0]]
        
        plumetypes.append(int(plumevalue))

    # Get statistics from this information
    daily_data['explained plumes'] = len(np.array(plumetypes)[np.array(plumetypes) != 1])
    daily_data['explained by gfed'] = len(np.array(plumetypes)[np.array(plumetypes) == 11])
    daily_data['explained by edgar'] = len(np.array(plumetypes)[np.array(plumetypes) == 101])
    
    return daily_data


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
CAMS_path = os.path.join(basepath, '02_Store_CAMS_data' + os.path.sep) # Path to CAMS .nc files


# Set modelparams
windowsizes = list(np.linspace(100, 200, 11).astype(int))
stepsizes = list(np.linspace(20, 90, 8).astype(int))
#buffersizes = list(np.linspace(2, 15, 14).astype(int))
#stdevs = list(np.round(np.linspace(0.2, 3, 15), 1))

# All loops that will have to be executed
iterations = list(itertools.product(stepsizes, windowsizes))

# Remove invalid combinations
invalid_list = []
for i in range(len(iterations)):
    if iterations[i][0] >= iterations[i][1]:
        invalid_list.append(i)
iterations = [j for i, j in enumerate(iterations) if i not in invalid_list]


# Set variables
estimator = np.zeros((len(stepsizes), len(windowsizes))) # y-hat (total TROPOMI identifications)
estimate = np.zeros((len(stepsizes), len(windowsizes))) # y (explained/correct TROPOMI identifications)


# ========================================================
# RUN ALGORITHM
# ========================================================
# Notify the algorithm has started
starttime = datetime.now()
print(f'Algorithm started at: {starttime}')

# Redefine model behaviour settings
apply_land_sea_mask = True
incorporate_cams = True
compare_gfed_edgar = True
apply_overlap_filter = False # TODO!
use_wind_rotations = False # TODO!
gen_fig_wind_vector = False
if use_wind_rotations:
    gen_fig_wind_vector = True
    
# Initialize other required parameters
boundaries, target_lon, target_lat, files = init.InitializingParameters(lat_min, lat_max, lon_min, lon_max, lonres, latres, basepath)

# Designated breaks (for progress updates)
breaks = start = int(len(iterations)/20)
breaklist = []
for i in range(20):
    breaklist.append(breaks)
    breaks += start

# Read/open the data
daily_data = init.CollectingData(boundaries, target_lon, target_lat, files, basepath, CAMS_path, apply_land_sea_mask, use_wind_rotations, incorporate_cams)

# Iterate over all model parameter combinations
for count, i in enumerate(iterations):
    # params = [
        # buffersize (radius of buffer around TROPOMI plumes),
        # st.devs (minimum amount of st.devs within identification frames),
        # windowsize (size of moving window frame in grid cells),
        # stepsize (steps between each moving window frame in grid cells)
        # ]
    params = [5, 1.6, i[1], i[0]]
    
    # Detection algorithm
    daily_data = init.Detection(params, daily_data, boundaries, GFED_path, \
                                EDGAR_path, lonres, latres, apply_overlap_filter, \
                                    use_wind_rotations, incorporate_cams, \
                                        compare_gfed_edgar)
    
    # Label the plumes
    labels, nlabels = ndimage.label(daily_data[0]['plumes_explained'])
    
    # Retrieve statistics from plume detection result
    stats = GetStats(daily_data[0], labels, nlabels)
    
    # Identify locations indices to place data in array
    x = np.where(stepsizes == iterations[count][0])[0]
    y = np.where(windowsizes == iterations[count][1])[0]
    
    # Append these statistics to array
    estimator[x,y] = daily_data[0]['total_plumes']
    estimate[x,y] = daily_data[0]['explained plumes']
    
    # Print progress update to console (steps of 5%)
    if count in breaklist:
        index = np.where(np.array(breaklist) == count)[0][0] + 1
        print(f'{5*index}% complete')
        print('going to sleep for 5 minutes...')
        time.sleep(300)

print('100% complete')
print(f'Total time elapsed: {datetime.now()-starttime}')

# Calculate Root mean Square Deviation and initialize X, Y coorindates
X, Y = np.meshgrid(windowsizes, stepsizes)

# Where is difference TROPOMI and EXPLAINED the smallest -> Root mean square deviation?
difference = estimator-estimate
RMSD = np.sqrt(np.divide(((difference)**2),estimator, where=estimator!=0))

# Tips from Sander (score function)
y_pos = estimate # positive identifications
y_neg = estimator-estimate # negative identifications
bufferarea = np.pi * (Y**2)

score = (y_pos - y_neg) / bufferarea

# ============= CREATE OUTPUTS ==============
# Store data as csv's
np.savetxt(os.path.join(os.getcwd() + rf"\RMSD.csv"), RMSD, delimiter=",")
np.savetxt(os.path.join(os.getcwd() + rf"\score.csv"), score, delimiter=",")

# Plot a pcolormesh (2D) image of results, saved in scripts workingdir
xlabel = 'window length and width (x10 km)'
ylabel = 'window step size (x10 km)'
cblabel = 'RMSD'

Plot2D(X, Y, ma.array(RMSD, mask=(RMSD == 0)), xlabel=xlabel, ylabel=ylabel, cblabel=cblabel)

# Now also plot score function
cblabel = 'score'
Plot2D(X, Y, ma.array(score, mask=(score == 0)), xlabel=xlabel, ylabel=ylabel, cblabel=cblabel, colormap='RdYlGn')
