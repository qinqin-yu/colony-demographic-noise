#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 18:57:59 2021

@author: qinqinyu
Calculates the relative MSD values by normalizing the absolute MSD value by the
average of the WT MSDs on that plate. Saves the relative MSDs into 
'data/relative_msd_50um/*.csv'.
"""
import csv
import numpy
import matplotlib.pyplot as plt
import re
import os
from glob import glob
plt.rcParams.update({'font.size': 18})

def find_substrings(test_string, substring):
    return [m.start() for m in re.finditer(substring, test_string)]

def import_data_msd_filtered(path, msd_file):
    #Import avg MSD as list
    msd_filename = path + 'extrapolated_msd_50um/' + msd_file
    excluded_filename = path + 'excluded_colonies/' + msd_file
    
    csvfile = open(msd_filename, 'r')
    csvreader = csv.reader(csvfile, delimiter = ",")
    x = list(csvreader)
    msd_all = numpy.array(x).astype('float')

    csvfile = open(excluded_filename, 'r')
    csvreader = csv.reader(csvfile, delimiter = ",")
    toexclude_list = list(csvreader)
    toexclude_array = numpy.array(toexclude_list[1:])
    toexclude_position = toexclude_array[:,1].astype('int')

    for i in range(len(msd_all)):
        if i in toexclude_position:
            msd_all[i] = 'nan'
    
    msd = msd_all[:,0]
    msd_err = msd_all[:,1]
    
    return msd, msd_err
    
def get_relative_drift_value(msd, msd_err, wt_ind):
    # Calculates the drift value relative to the wild type weighted average
    # Get values of wild type drift values, and remove nans
    msd_wt = msd[wt_ind]
    msd_err_wt = msd_err[wt_ind]

    msd_wt = msd_wt[~numpy.isnan(msd_wt)]
    msd_err_wt = msd_err_wt[~numpy.isnan(msd_err_wt)]

    # Calculate weighted mean of wild type values and standard error of the mean
    msd_wmean_wt = numpy.average(msd_wt, weights = 1/(msd_err_wt)**2)
    msd_wmean_ste_wt = numpy.sqrt(1/(numpy.sum(1/(msd_err_wt)**2)))

    # Calculate the relative drift value (scaled to the wild type mean) and error from error propagation
    msd_rel = msd/msd_wmean_wt
    msd_err_rel = msd_rel*numpy.sqrt(((msd_err**2)/msd**2) + ((msd_wmean_ste_wt**2)/(msd_wmean_wt**2)))
    
    return msd_rel, msd_err_rel

data_folder = '../../data/'

strain_list = data_folder +'main_experiments/strain_list.csv'
csvfile = open(strain_list, 'r')
csvreader = csv.reader(csvfile, delimiter = ",")
x = list(csvreader)
strains_all = numpy.array(x)

# Get all non-hidden folders
global_path = data_folder + 'extrapolated_msd_50um/'
files = glob(global_path+'DE*.csv')

for file in files:
    filename = os.path.basename(file)
    msd, msd_err = import_data_msd_filtered(data_folder, filename)
    plate_label = filename.split('.')[0][:-1]
    plate_strains = strains_all[strains_all[:,0]==plate_label,:]
    wt_ind = numpy.where(plate_strains[:,2] == 'WT')[0]
    msd_rel, msd_err_rel = get_relative_drift_value(msd, msd_err, wt_ind)

    msd_rel_combined = numpy.transpose([msd_rel, msd_err_rel])
    numpy.savetxt(data_folder + 'relative_msd_50um/'+filename, msd_rel_combined, delimiter=",")