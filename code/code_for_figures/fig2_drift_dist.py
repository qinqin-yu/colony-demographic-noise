#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 16:49:07 2021

Plots the distribution of MSDs for WT and KOs, calculates the probability that
the distributions are different and the medians are different. 

Note that the probability that the medians are different changes slightly every
time because it is calculated by drawing a set of random samples.
@author: qinqinyu
"""

import numpy
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import ks_2samp, norm
plt.rcParams.update({'font.size': 25})

data_folder = '../../data/main_experiments/'
figure_folder = '../../figures/'

msd_wmean_file = data_folder + 'relative_msd_phenotypic_traits.csv'
df = pd.read_csv(msd_wmean_file)

df = df.dropna(subset = ['Weighted mean of relative MSD'])
df = df.reset_index(drop = True)

experiment = df['Plate'].values
genotype = df['Genotype'].values
wmean = df['Weighted mean of relative MSD'].values
wmean_ste = df['Standard error of relative MSD'].values

rand_de_wt_ind = [i for i in range(len(df)) if 'DE' in experiment[i] and experiment[i]!='DE5' and genotype[i]=='WT']
rand_de_ko_ind = [i for i in range(len(df)) if 'DE' in experiment[i] and experiment[i]!='DE5' and genotype[i]!='WT']

n_wt = len(rand_de_wt_ind)
n_ko = len(rand_de_ko_ind)

x = wmean[rand_de_wt_ind]
y = wmean[rand_de_ko_ind]

(D, p) = ks_2samp(x, y)

print('Median WT = ' + str(numpy.median(x)))
print('Median KO = ' + str(numpy.median(y)))

print('Mean WT = ' + str(numpy.mean(x)))
print('Mean KO = ' + str(numpy.mean(y)))

print('Var WT = ' + str(numpy.var(x)))
print('Var KO = ' + str(numpy.var(y)))

print('Number of KOs lower than lowest WT = ' + str(numpy.count_nonzero(y<numpy.min(x))))
quantiles_input = [0, 25, 50, 75, 100]
quantiles = numpy.percentile(y, quantiles_input)
print('Quantiles = ' + str(quantiles))

# Calculate how significant the medians are different from one another
# Jackknife
n_subsets = 100000
subsample = numpy.empty((n_wt, n_subsets))
subsample.fill(numpy.nan)
for q in range(n_subsets):
    subsample[:, q] = numpy.random.choice(y, size = n_wt, replace = False)
subsample_median = numpy.median(subsample, axis = 0)
subsample_median_sorted = numpy.sort(subsample_median)
lower_ci = subsample_median_sorted[int(n_subsets*0.025)]
upper_ci = subsample_median_sorted[int(n_subsets*0.975)]
print('Confidence interval for KO median = ' + '[' + str(lower_ci) + ', ' + str(upper_ci) + ']')
p_medians = numpy.sum(subsample_median>numpy.median(x))/n_subsets

# Fit Gaussian to WT dist
mean, std = norm.fit(x)
mean_ko, std_ko = norm.fit(y)
xnorm = numpy.linspace(0, 2, num = 1000)
ynorm = norm.pdf(xnorm, mean, std)
ynorm_cdf = norm.cdf(xnorm, mean, std)

plt.figure(figsize = (13,10))
plt.hist(x, bins = numpy.linspace(0,2,30), density = True, cumulative = 1, histtype = 'step', label = 'wild type, n = ' + str(n_wt), linewidth = 2)
plt.hist(y, bins = numpy.linspace(0,2,30), density = True, cumulative = 1, histtype = 'step', label = 'single gene knockouts \n n = ' + str(n_ko), linewidth = 2)
plt.plot(xnorm, ynorm_cdf, color = 'tab:blue', linestyle = ':')
plt.plot([numpy.median(x), numpy.median(x)], [0,1], color = 'tab:blue')
plt.plot([numpy.median(y), numpy.median(y)], [0,1], color = 'tab:orange')

plt.xlabel('Relative mean squared displacement, MSD/$\\langle$MSD$_{WT}\\rangle$ \n (L = 50$\mu m$)')
plt.ylabel('Cumulative distribution function')
plt.xlim([0,1.5])
plt.ylim([0,1])
plt.text(0.1, 0.1, '$p_{dist}$ = ' + str(numpy.format_float_scientific(p, precision = 1)))
plt.text(1.05, 0.1, '$p_{median}$ = ' + "{:.1e}".format(p_medians))
plt.legend(loc = 'upper left')

plt.tight_layout()
plt.savefig(figure_folder + 'fig2_drift_dist.pdf')
    