#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 17:29:59 2021

@author: qinqinyu

Fits the establishment probability as a function of fitness coefficient. 
Uses fit to get the establishment probability at fixed fitness coefficients
across all genotypes. Plots Figure 4b: establishment probability vs the MSD
for each genotype.
"""

# Plotting only Keio collection strains

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="whitegrid",font_scale=1.2)

data_folder = '../../data/establishment_probability/'
figure_folder = '../../figures/'


genotypes_filename = data_folder + '/genotypes_order.csv'
genotypes = pd.read_csv(genotypes_filename, names = ['Genotype'])

pest_all = pd.read_csv(data_folder + 'establishment_probability_raw.csv', index_col = 0)
fitness_all = pd.read_csv(data_folder + 'fitness.csv', index_col = 0)

filename = data_folder + '../validation/validation.csv'
df_validation = pd.read_csv(filename)
fits = pd.DataFrame(columns = ['genotype', 's0.05', 's0.1', 's0.2', 's0.05_err', 's0.1_err', 's0.2_err'])

s = [0.05, 0.1, 0.15, 0.2]

fig, ax = plt.subplots(3, 3, figsize = (10, 8))

i = 0
j = 0

for genotype in genotypes['Genotype'].values:
    fitness = fitness_all[fitness_all['Genotype']==genotype]
    fitness_tiled = np.tile(fitness['Fitness coefficient'].values, 3)
    pest = pest_all[pest_all['Genotype']==genotype]
    pest_p0_002 = pest[pest['Expected initial mutant fraction'] == 0.002]['Establishment probability'].values
    pest_p0_005 = pest[pest['Expected initial mutant fraction'] == 0.005]['Establishment probability'].values
    pest_p0_020 = pest[pest['Expected initial mutant fraction'] == 0.020]['Establishment probability'].values
    pest_ordered = np.concatenate((pest_p0_002, pest_p0_005, pest_p0_020))
    
    idx = np.logical_and(np.logical_and(fitness_tiled<0.5, fitness_tiled>-0.1), ~np.isnan(pest_ordered))

    x = fitness_tiled[idx]
    y = pest_ordered[idx]

    if len(x)>0:
        z, cov = np.polyfit(x, y, 1, full = False, cov = True)
        zerr = np.sqrt(np.diag(cov))
        p = np.poly1d(z)
        ppp = np.poly1d(z + zerr)
        pmm = np.poly1d(z - zerr)
        ppm = np.poly1d(z + np.array([1, -1])*zerr)
        pmp = np.poly1d(z + np.array([-1, 1])*zerr)
        
        dpp = ppp(s)
        dmm = pmm(s)
        dpm = ppm(s)
        dmp = pmp(s)

        d = np.array([dpp, dmm, dpm, dmp])
        dmax = np.max(d, axis = 0)
        dmin = np.min(d, axis = 0)
        
        err = (dmax-dmin)/2
        
        xfit = np.linspace(-0.1, 0.5, 10)
        fit = pd.DataFrame([np.append(np.append(genotype, p(s)),err)], columns = ['genotype', 's0.05', 's0.1', 's0.15', 's0.2', 's0.05_err', 's0.1_err', 's0.15_err', 's0.2_err'])
        fits = fits.append(fit)
        
        if ~np.isnan(z[0]):
            ax[i,j].scatter(fitness['Fitness coefficient'], pest_p0_002, label = '$p_0$ = 0.002')# + str(round(pest[pest['p0_expected'] == 0.002]['p0_actual'].values[0], 4)))
            ax[i,j].scatter(fitness['Fitness coefficient'], pest_p0_005, label = '$p_0$ = 0.005')# + str(round(pest[pest['p0_expected'] == 0.005]['p0_actual'].values[0], 4)))
            ax[i,j].scatter(fitness['Fitness coefficient'], pest_p0_020, label = '$p_0$ = 0.02')# + str(round(pest[pest['p0_expected'] == 0.020]['p0_actual'].values[0], 4)))

            ax[i,j].plot(xfit, p(xfit), 'k-', label = 'linear fit')

            ax[i,j].fill_between([-0.5, -0.1], 0, 1, color = 'gray', alpha = 0.2)
            ax[i,j].fill_between([0.5, 1.5], 0, 1, color = 'gray', alpha = 0.2)
            ax[i,j].set_xlim([-0.5, 1.5])
            ax[i,j].set_ylim([0, 1])
            ax[i,j].set_title(genotype)
    
            if (j+1)%3 == 0:
                i += 1
                j = 0
            else:
                j += 1
fits.reset_index(inplace = True, drop = True)

ax[2,1].set_xlabel('fitness')
ax[1,0].set_ylabel('establishment probability')

ax[1,2].legend(loc = (1.05, 0.15))
plt.tight_layout()
plt.savefig(figure_folder + 'si_establishment_probability_fit.pdf')
plt.show()

merged = df_validation.merge(fits, on = 'genotype')

merged_nowt = merged[(merged['genotype'] != 'BW25113') & (merged['genotype'] != 'DH5alpha') & (merged['genotype'] != 'MG1655')]
merged_wt = merged[merged['genotype'] == 'BW25113']

plt.errorbar(merged_nowt['MSD (L = 50um) (microns)'].values.astype('float'), merged_nowt['s0.05'].values.astype('float'), yerr = merged_nowt['s0.05_err'].values.astype('float'), xerr = merged_nowt['MSD standard error (L = 50um) (microns)'].values.astype('float'), linestyle = '', marker = 'o', label = 's = 0.05', alpha = 0.8)
plt.errorbar(merged_nowt['MSD (L = 50um) (microns)'].values.astype('float'), merged_nowt['s0.1'].values.astype('float'), yerr = merged_nowt['s0.1_err'].values.astype('float'), xerr = merged_nowt['MSD standard error (L = 50um) (microns)'].values.astype('float'), linestyle = '', marker = 'o', label = 's = 0.1', alpha = 0.8)
plt.errorbar(merged_nowt['MSD (L = 50um) (microns)'].values.astype('float'), merged_nowt['s0.15'].values.astype('float'), yerr = merged_nowt['s0.15_err'].values.astype('float'), xerr = merged_nowt['MSD standard error (L = 50um) (microns)'].values.astype('float'), linestyle = '', marker = 'o', label = 's = 0.15', alpha = 0.8)

plt.errorbar(merged_wt['MSD (L = 50um) (microns)'].values.astype('float'), merged_wt['s0.05'].values.astype('float'), yerr = merged_wt['s0.05_err'].values.astype('float'), xerr = merged_wt['MSD standard error (L = 50um) (microns)'].values.astype('float'), linestyle = '', marker = 'o', mec = 'tab:blue', mfc = 'white', color = 'tab:blue', alpha = 0.8)
plt.errorbar(merged_wt['MSD (L = 50um) (microns)'].values.astype('float'), merged_wt['s0.1'].values.astype('float'), yerr = merged_wt['s0.1_err'].values.astype('float'), xerr = merged_wt['MSD standard error (L = 50um) (microns)'].values.astype('float'), linestyle = '', marker = 'o', mec = 'tab:orange',  mfc = 'white', color = 'tab:orange', alpha = 0.8)
plt.errorbar(merged_wt['MSD (L = 50um) (microns)'].values.astype('float'), merged_wt['s0.15'].values.astype('float'), yerr = merged_wt['s0.15_err'].values.astype('float'), xerr = merged_wt['MSD standard error (L = 50um) (microns)'].values.astype('float'), linestyle = '', marker = 'o', mec = 'tab:green',  mfc = 'white', color = 'tab:green', alpha = 0.8)

D = np.logspace(-3, 1, 100)
phi0 = 2*10**-2
L = 50 # microns
colors = ['tab:blue', 'tab:orange', 'tab:green']

plt.errorbar([], [], yerr = [], xerr = [], marker = 'o', color = 'k', mfc = 'white', mec = 'k', linestyle = '', label = 'WT')
plt.errorbar([], [], yerr = [], xerr = [], marker = 'o', color = 'k', linestyle = '', label = 'single gene deletion')

plt.legend()
plt.xlim([0, 3])
plt.ylim([0, 0.4])
plt.xlabel('Mean squared displacement, MSD \n (L = 50$\mu m$) [$\mu m^2$]')
plt.ylabel('Establishment probability, $p_{est}$')
# plt.tight_layout()
plt.savefig(figure_folder + 'fig4b_establishment_prob.pdf')
plt.show()