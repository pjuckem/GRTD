# coding: utf-8

'''
This notebook is used to read-in a list of 2D ascii arrays and convert them to a compressed 3D numpy array file.
There are pros/cons to writing NP arrays as ascii or npz files, so this provides flexibility to feed the desired
format into various other codes.  This code also allows you to add layers to the 3D array and apply a default value
to those layers.
'''

__author__ = 'Juckem'
import os
import numpy as np
import pandas as pd
import gdal
import gen_mod_functions as gm  # from Starn's general model notebooks


# Input
# List the files in the order you want them converted into the 3D array (top of model first)
#ascii_file_list = ['../Models/FWP5L_hK/por_Layer_1_AWHb8_1.ref', #'../Models/FWP5L_hK/por_Layer_2_AWHb8_1.ref',
#                   '../Models/FWP5L_hK/por_Layer_3_AWHb8_1.ref']
ascii_file_list = ['../Models/FWP5L_hK/por_Layer_1_SM96_v7_quad1D.ref', 
'../Models/FWP5L_hK/por_Layer_2_SM96_v7_quad1D.ref',
'../Models/FWP5L_hK/por_Layer_3_SM96_v7_quad1D.ref']
zone_array_src = '../Models/FWP5L_hK/calibration_area_noBR.npz'  # optional. If omitted, stats run for entire model; if included should be 0 & 1s as per Ibound

ref_raster = 'D:/PFJDATA/Projects/NAWQA/Modeling/general-models/FWP_ref_raster.tif'  # New raster produced with proper GT parameters!

# Settings
nrow = 930
ncol = 650
total_layers = 5
value4added_layers = 0.2

# Output
#npz_file = '../Models/FWP5L_hK/por_5LhK_AWHb8_1'
npz_file = '../Models/FWP5L_hK/por_5LhK_SM96_v7_quad1D'
#stats_file = '../Models/FWP5L_hK/por_5LhK_AWHb8_1_stats.dat'
stats_file = '../Models/FWP5L_hK/por_5LhK_SM96_v7_quad1D_stats.dat'
#rasterfiles = ['../Models/FWP5L_hK/por_5LhK_layer1_AWHb8_1.tiff',
#               '../Models/FWP5L_hK/por_5LhK_layer2_AWHb8_1.tiff',
#               '../Models/FWP5L_hK/por_5LhK_layer3_AWHb8_1.tiff']
rasterfiles = ['../Models/FWP5L_hK/por_5LhK_layer1_SM96_v7_quad1D.tiff',
               '../Models/FWP5L_hK/por_5LhK_layer2_SM96_v7_quad1D.tiff',
               '../Models/FWP5L_hK/por_5LhK_layer3_SM96_v7_quad1D.tiff']


# Main

ras = gdal.Open(ref_raster)
shapeproj = ras.GetProjection()
gt = ras.GetGeoTransform()


array = np.ones((total_layers, nrow, ncol))
for i, file in enumerate(ascii_file_list):
    arr = np.loadtxt(file)
    if arr.shape[0] == nrow and arr.shape[1] == ncol:
        array[i] = arr
        gm.make_raster(rasterfiles[i], arr, ncol, nrow, gt, shapeproj, np.nan)

missing_layers = total_layers - 1 - i # -1 to convert from 1-based to zero-based
for l in range(missing_layers):
    array[i+l+1] = array[i+l+1] * value4added_layers

np.save(npz_file, array)


# now perform some stats for the area that will get particles, and print-out the results.
if 'None' not in zone_array_src:
    ext = os.path.splitext(zone_array_src)[1].lower()
    if ext == '.csv':
        zones = pd.read_csv(zone_array_src, header=0)
        zones = zones.unstack().values.reshape(nlay, nrow, ncol)
        print('Zones read from csv file')
    elif (ext == '.zon') | (ext == '.zone'):
        # option to read zones from a MODFLOW zone package file not implemented
        print('Zones read from MODFLOW zone array')
        pass
    elif ext == '.npz':
        print('Zones read from compressed numpy array object')
        # np.save()  # Jeff, why was this added on June 21?  It bombs.
        d = np.load(zone_array_src)
        if len(d.items()) != 1:
            print('There should only be one item type in the npz file and ')
            print('it should be an array dimensioned as (nlay, nrow, ncol)')
        else:
            zones = d.items()[0][1]
else:
    print('No zone information read')
    zones = np.ones((nlay * nrow * ncol))

array = array * zones  # Treating zones as Ibound, so only include 0s and 1s!!!!!
array = np.where(array==0, np.nan, array)
meanarr = np.nanmean(array, axis=(1,2))
medarr = np.nanmedian(array, axis=(1,2))
stdarr = np.nanstd(array, axis=(1,2))
data = np.stack((meanarr, medarr, stdarr), axis=0)
columns = ['mean', 'median', 'standard deviation']
rows = [r+1 for r in range(5)]
df = pd.DataFrame(data=data.T, index=rows, columns=columns)
df.to_csv(stats_file, index_label='Layer')

print('Average porosity of all layers/cells:  {}'.format(np.nanmean(meanarr)))