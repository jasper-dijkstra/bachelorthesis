# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 11:12:14 2020

@author: jaspd
"""

import os


def DefineAndCreateDirectory(targetDirectory):
    # Check if directory exists, and else create it
    if not os.path.isdir(targetDirectory):
        os.makedirs(targetDirectory)
        
        # Add notification to log file
    
    # Make sure path ends with separator (//)
    if not targetDirectory.endswith(os.path.sep):
        targetDirectory += os.path.sep
    
    return targetDirectory


def ListCSVFilesInDirectory(inputDirectory, maxfiles=None):
    # Check if directory exists, and else create it
    if not os.path.isdir(inputDirectory):
        print('This directory does not exist!')
        # Add notification to logfile
        return

    files = []
    for file in os.listdir(inputDirectory):
        if maxfiles != None:
            if len(files) == maxfiles: break # Limit the amount of input days
        if file.endswith(".csv"): files.append(os.path.join(inputDirectory, file))
        else:
            continue
    
    return files
    