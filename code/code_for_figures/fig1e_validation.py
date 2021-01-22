#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 09:47:09 2021

@author: qinqinyu

Calculates the fraction of diversity preserved (FDP) from the number of 
sectros and inoculum area. Plots the average FDP vs MSD for each genotype, 
fits inverse square root relationship from Hallatschek et al, Evolution, 2010. Saves mean FDP and MSD
for each genotype to a csv.
"""

import csv
import numpy
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import statistics
from scipy.stats import pearsonr, spearmanr
import math
import glob
import os
import re
import pandas as pd
from scipy.optimize import curve_fit
plt.rcParams.update({'font.size': 25})

data_folder = '../../data/validation/'
figure_folder = '../../figures/'

def func(x,a):
    return a/numpy.sqrt(x)

genotypes = ['BW25113',
            'MG1655',
            'DH5$\\alpha$',
            '$\Delta gpmI$',
             '$\Delta recB$',
            '$\Delta pgm$',
            '$\Delta tolQ$',
            '$\Delta ychJ$',
            '$\Delta lpcA$',
            '$\Delta dsbA$',
            '$\Delta rfaF$',
            '$\Delta tatB$']

genotypes_df = ['BW25113',
            'MG1655',
            'DH5alpha',
            'gpmI',
             'recB',
            'pgm',
            'tolQ',
            'ychJ',
            'lpcA',
            'dsbA',
            'rfaF',
            'tatB']

# Comparison of bead MSD with the fraction of diversity preserved
# Filtering out measurements where msd_err>0.5*msd

file1 = data_folder + 'msd_50um.csv'
csvfile = open(file1, 'r')
csvreader = csv.reader(csvfile, delimiter = ",")
temp = list(csvreader)
msd = numpy.array(temp).astype('float')[0:96,0]
msd_err = numpy.array(temp).astype('float')[0:96,1]
msd = numpy.transpose(msd)
msd_err = numpy.transpose(msd_err)

file2 = data_folder + 'number_sectors.csv'
csvfile = open(file2, 'r')
csvreader = csv.reader(csvfile, delimiter = ",")
sectors = numpy.zeros((96,))
sectors_err = numpy.zeros((96,))
rownum = 0
for row in csvreader:
    if rownum>0:
        sectors[rownum-1] = row[1]
        sectors_err[rownum-1] = row[2]
    rownum+=1
    
df_areas = pd.read_csv(data_folder + 'inoculum_area.csv')
areas = df_areas['Area'].values
radii = numpy.sqrt(areas/math.pi)
conv_8x = 0.1240 #0.1240 pixels/um for 8x;
radii_um = radii/conv_8x
cell_size = 1.7
num_cells = 2*math.pi*radii_um/cell_size
fdp = 2*sectors/num_cells  #factor of 2 comes from the probaiblity that two neighboring cells in the inoculum have the same color
fdp_err = 2*sectors_err/num_cells

#Filter out high error MSD measurements
for j in range(len(msd)):
    if msd_err[j]>msd[j]:
        msd[j] = 'nan'
        msd_err[j] = 'nan'

numgenotypes = int(len(sectors)/8)

xavg = numpy.zeros((numgenotypes, 1))
yavg = numpy.zeros((numgenotypes, 1))
xstd = numpy.zeros((numgenotypes, 1))
xste = numpy.zeros((numgenotypes, 1))
ystd = numpy.zeros((numgenotypes, 1))
yste = numpy.zeros((numgenotypes, 1))

#Calculate weighted averages
for i in range(numgenotypes):
    ind = list(range(i,96,12))
    xvals = msd[ind]
    yvals = fdp[ind]
    xerr = msd_err[ind]
    yerr = fdp_err[ind]
  
    xvals = xvals[~numpy.isnan(xvals)]
    yvals = yvals[(~numpy.isnan(yvals))]
    
    xerr = xerr[~numpy.isnan(xerr)]
    xweight = numpy.divide(1, numpy.square(xerr))
    
    if len(xvals)>6:
        xavg[i] = numpy.sum(numpy.multiply(xvals, xweight))/numpy.sum(xweight)
        xste[i] = numpy.sqrt(numpy.divide(1, numpy.sum(xweight)))#/len(xweight)**(1/2)
    else:
        xavg[i] = numpy.nan
        xste[i] = numpy.nan
    
    yerr = yerr[~numpy.isnan(yerr)]
    yweight = numpy.divide(1,numpy.square(yerr))
    
    yavg[i] = numpy.sum(numpy.multiply(yvals, yweight))/numpy.sum(yweight)
    yste[i] = numpy.sqrt(numpy.divide(1, numpy.sum(yweight)))#/len(yweight)**(1/2)

ind = numpy.logical_and(~numpy.isnan(xavg), ~numpy.isnan(yavg))

xavg_nonan = xavg[ind]
xste_nonan = xste[ind]
yavg_nonan = yavg[ind]
yste_nonan = yste[ind]
genotypes_array = numpy.array(genotypes_df)
genotypes_nonan = genotypes_array[ind.flatten()]

popt,pcov = curve_fit(func, xavg_nonan, yavg_nonan)

perr = numpy.sqrt(numpy.diag(pcov))

chisq = sum((1/(yste_nonan)**2)*(yavg_nonan - func(xavg_nonan,popt))**2)

reduced_chisq = chisq/(len(xavg_nonan)-1)

xfit = numpy.linspace(0,3,100)

fig,ax = plt.subplots(1,1,figsize = (10,7))
ax.errorbar(xavg_nonan, yavg_nonan*10**3, yste_nonan*10**3, xste_nonan, c = 'k', fmt = 'o')
ax.errorbar(xavg[~ind], yavg[~ind]*10**3, yste[~ind]*10**3, xste[~ind], c = 'gray', fmt = 'o')
ax.errorbar(xavg[7], yavg[7]*10**3, yste[7]*10**3, xste[7], c = 'tab:blue', fmt = 'o')
ax.errorbar(xavg[0], yavg[0]*10**3, yste[0]*10**3, xste[0], c = 'tab:orange', fmt = 'o')

spearman = spearmanr(xavg_nonan, yavg_nonan)
ax.plot(xfit, func(xfit, popt)*10**3, color = 'k', label = 'fit: $y = \\frac{a}{\sqrt{x}}$')#' \n a = ' + str("{:.1e}".format(popt[0])) + '+/-' + str("{:.1e}".format(perr[0])) + '\n$\chi^2_\\nu = $' + str(round(reduced_chisq, 2)))
ax.legend(loc = 'upper right')

for i in range(len(genotypes)):
    ax.text(xavg[i]+0.02, yavg[i]*10**3+0.5, genotypes[i])
ax.text(0.25, 1, '$\\rho$ = ' + str(round(spearman[0], 2)) + ' (p = ' + str(round(spearman[1], 3)) + ')')
ax.set_xlim([0,3])
ax.set_ylim([0,24])
ax.set_xlabel('Mean squared displacement, MSD \n (L = 50$\mu m$) [$\mu m^2$]')
ax.set_ylabel('Fraction of diveristy preserved \n ($\\times 10^{-3}$)')
ax.set_title('a = ' + str("{:.1e}".format(popt[0])) + '+/-' + str("{:.1e}".format(perr[0])) + ', $\chi^2_\\nu = $' + str(round(reduced_chisq, 2)))
plt.savefig(figure_folder + 'fig1e_validation.pdf')
plt.show()

df_coeff = pd.DataFrame(numpy.array([[popt[0], perr[0]]]), columns = ['a', 'aerr'])
# df_coeff.to_csv('/Volumes/PORTKeio3/20200224_0900_validation_a/summary_data/msd2fdp_conv.csv', index = False)

df_msd_fdp = pd.DataFrame(data = numpy.transpose(numpy.array([genotypes_nonan, xavg_nonan, xste_nonan, yavg_nonan, yste_nonan])), index = numpy.array(range(len(xavg_nonan))), columns = ['genotype', 'MSD (L = 50um) (microns)', 'MSD standard error (L = 50um) (microns)', 'FDP', 'FDP standard error'])
df_msd_fdp.to_csv(data_folder + 'validation.csv', index = False)
