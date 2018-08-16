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


# Input
# List the files in the order you want them converted into the 3D array (top of model first)
ascii_file_list = ['../Models/FWP5L_hK/por_Layer_1.ref', '../Models/FWP5L_hK/por_Layer_2.ref',
                   '../Models/FWP5L_hK/por_Layer_3.ref']

# Settings
nrow = 930
ncol = 650
total_layers = 5
value4added_layers = 0.2

# Output
npz_file = '../Models/FWP5L_hK/por_5L_hK.npz'

# Main
array = np.ones((total_layers, nrow, ncol))
for i, file in enumerate(ascii_file_list):
    arr = np.loadtxt(file)
    if arr.shape[0] == nrow and arr.shape[1] == ncol:
        array[i] = arr

missing_layers = total_layers - 1 - i # -1 to convert from 1-based to zero-based

for l in range(missing_layers):
    array[i+l+1] = array[i+l+1] * value4added_layers

np.save(npz_file, array)


