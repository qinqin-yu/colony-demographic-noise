#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 19:26:42 2021

@author: qinqinyu
"""

import os
from skimage import io, img_as_uint
from skimage.segmentation import clear_border
from skimage.filters import threshold_otsu
from skimage.filters import threshold_local
from skimage.measure import label
from skimage.morphology import remove_small_objects

def binarize_colony_images(path_folder, experiment):
    openpath = path_folder + 'images/' + experiment + '/'
    savepath = path_folder + 'binary/' + experiment + '/'

    # check if the folder already exists
    try:
        os.makedirs(savepath)
    except OSError:
        if not os.path.isdir(savepath):
            raise

    for filename in os.listdir(openpath):
        if filename[0]!='.':
            image = io.imread(openpath + filename)
            thresh = threshold_local(image, 601)#151)
            binary = image < thresh
            binary2 = clear_border(binary)
            colony = remove_small_objects(binary2, 30000)
            colony = img_as_uint(colony)
            io.imsave(savepath + filename, colony, check_contrast = False)
            
            
folders_to_test = []

path = '../../data/whole_colony_images/'
experiment = 'DE5d'
binarize_colony_images(path, experiment)