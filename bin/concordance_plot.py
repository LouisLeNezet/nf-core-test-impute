#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import gzip

x_maf = []
y_aggr = []
y_avg = []

#plt.rc('xtick', labelsize=20)
#plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)

with gzip.open('input.rsquare.grp.txt.gz','rt') as csvfile:
    plots = csv.reader(csvfile, delimiter=' ')
    for row in plots:
        y_aggr.append(float(row[4]))
        x_maf.append(float(row[2]))
	
fig, axs = plt.subplots(1)

axs.semilogx(x_maf,y_aggr, '--',  marker='o', lw=1, label = 'Aggregate r2', markersize=10, alpha=0.8)
   
mybins = [0.0002, 0.0005, 0.001, 0.005, 0.02, 0.05, 0.10, 0.20, 0.50]
labels = ['0.02', '0.05', '0.1', '0.5', '2', '5', '10', '20', '50']
axs.set_xlabel('Minor allele frequency (%)')
axs.set_ylabel('$r^2$ imputed vs true genotypes')
axs.grid(linestyle='--')
axs.legend(loc="lower right",prop={'size': 22})
axs.set_ylim([0.1,1.05])
axs.set_xlim([0.0005,0.6])
axs.set_xticks(mybins)
axs.set_xticklabels(labels)
axs.minorticks_off()
#axs.set_xticklabels(r_labels)
axs.minorticks_off()

axs.set_title('Imputation accuracy',fontsize=20)
fig.set_size_inches(16, 12)
fig.savefig('accplot.png',dpi=300)
