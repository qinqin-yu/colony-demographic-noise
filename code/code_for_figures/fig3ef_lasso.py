#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 17:12:05 2021

@author: qinqinyu

Calculates a Lasso regression model using the phenotypic traits to predict MSD.
Plots Figure 3e (linear model coefficients).
Calculates the predicted MSD from this model and plots Figure 3f.
"""
import numpy
from sklearn import linear_model
from sklearn.linear_model import LassoCV
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.utils import shuffle
import matplotlib.pyplot as plt
import time
import pandas as pd
from matplotlib.ticker import MultipleLocator
plt.rcParams.update({'font.size': 25})

data_folder = '../../data/main_experiments/'
figure_folder = '../../figures/'

filename = data_folder + 'relative_msd_phenotypic_traits.csv'
df = pd.read_csv(filename)
labels = ['Relative colony area mean', 'Relative front roughness mean', 'Relative growth layer depth mean', 'Relative colony thickness mean', 'Cell major axis', 'Cell minor axis', 'Cell aspect ratio', 'Cell volume', 'Cell surface area']
data = df[labels].values.astype('float')
target = df['Weighted mean of relative MSD'].values.astype('float')

experiment = df['Plate']
position = df['Position']

# Only looking at strains from DE experiments
data = data[[i for i in range(len(experiment)) if experiment[i][0:2] == 'DE' and position[i]!='WT']]
target = target[[i for i in range(len(experiment)) if experiment[i][0:2] == 'DE' and position[i]!='WT']]

# remove rows with nan values
data_isnan = numpy.isnan(data)
data_entry_isnan = data_isnan.any(axis = 1)

target_isnan = numpy.isnan(target)
joint = ~numpy.logical_or(data_entry_isnan, target_isnan)
joint_ind = numpy.where(joint)[0]

data_nonan = data[joint_ind,:]
target_nonan = target[joint_ind]

# Standardize
scaler = StandardScaler()
scaler.fit(data)
data_transformed = scaler.transform(data_nonan)

# Shuffle data
data_transformed_shuffled, target_nonan_shuffled = shuffle(data_transformed, target_nonan)

# This is to avoid division by zero while doing numpy.log10
EPSILON = 1e-8#1e-4

# LassoCV: coordinate descent

# Compute paths
print("Computing regularization path using the coordinate descent lasso...")
t1 = time.time()
model = LassoCV(cv=10).fit(data_transformed_shuffled, target_nonan_shuffled)
t_lasso_cv = time.time() - t1

# Display results
m_log_alphas = -numpy.log10(model.alphas_ + EPSILON)

plt.figure(figsize = (8, 5))
ymin, ymax = 0, 0.06
plt.plot(m_log_alphas, model.mse_path_, ':')
plt.plot(m_log_alphas, model.mse_path_.mean(axis=-1), 'k',
         label='Average across the folds', linewidth=2)
plt.axvline(-numpy.log10(model.alpha_ + EPSILON), linestyle='--', color='k',
            label='alpha: CV estimate')

plt.legend()

plt.xlabel('-log(alpha)')
plt.ylabel('Mean square error')
plt.title('Mean square error on each fold: \n coordinate descent '
          '(train time: %.2fs)' % t_lasso_cv)
plt.axis('tight')
plt.ylim(ymin, ymax)

clf = linear_model.Lasso(alpha = model.alpha_)
clf.fit(data_transformed, target_nonan)
predicted = clf.predict(data_transformed)
r2 = r2_score(target_nonan, predicted)

xlim = [0, numpy.max(target_nonan)+0.1]
ylim = xlim.copy()

plt.rcParams.update({'font.size': 25})
plt.figure(figsize = (7,7))
plt.scatter(predicted, target_nonan, color = 'k', edgecolors = None, alpha = 0.5)
plt.plot(xlim, ylim, 'k-')
plt.xlabel('Predicted MSD/MSD$_{WT}$')
plt.ylabel('True MSD/MSD$_{WT}$')
plt.title('$\\alpha$ = ' + str(numpy.format_float_scientific(model.alpha_, precision = 1)))
plt.text(numpy.max(target_nonan)-0.7, 0.1, 'R$^2$ = ' + str(round(r2,3)))
plt.xlim(xlim)
plt.ylim(ylim)

plt.savefig(figure_folder + 'fig3f_lasso_prediction.pdf')
plt.show()

labels = numpy.array(['Colony area', 
          'Front roughness', 
          'Growth layer depth',
          'Colony thickness',
          'Major axis length', 
          'Minor axis length', 
          'Aspect ratio',
          'Volume',
          'Surface area'])

contribution = abs(clf.coef_)/sum(abs(clf.coef_))
sort_ind = numpy.argsort(contribution)

y_pos = numpy.arange(len(clf.coef_[sort_ind]))

plt.rcParams.update({'font.size': 20})
fig, ax = plt.subplots(figsize = (5,5))
ax.barh(y_pos, clf.coef_[sort_ind], align='center', color = 'k', alpha = 0.5)
ax.plot([0,0],[max(y_pos)+1, min(y_pos)-1], 'k-')

ax.xaxis.set_major_locator(MultipleLocator(0.1))
ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.set_yticks(y_pos)
ax.set_yticklabels(labels[sort_ind])
ax.set_xlabel('Linear regression coefficient')
ax.set_ylim([min(y_pos)-1, max(y_pos)+1])
ax.set_xlim([-0.15, 0.15])

plt.grid(which = 'both', axis = 'x')

plt.savefig(figure_folder + 'fig3e_lasso_coeffs.pdf')

plt.show()