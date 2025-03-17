#!/usr/bin/env python

import argparse

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
font = {'family' : 'Arial',
        'size'   : 12}
matplotlib.rc('font', **font)

def get_args():
    parser = argparse.ArgumentParser(description='Plot result of MCL clustering')
    parser.add_argument("summary_filename", help="Filename of summary file for region")
    parser.add_argument("png_filename", help="Figure png filename")
    parser.add_argument("pdf_filename", help="Figure pdf filename")
    parser.add_argument("region", help="Region name")
    return parser.parse_args()

args = get_args()

cluster_info = pd.read_csv(args.summary_filename)

fig, ax = plt.subplots(3, 2, sharex = True, figsize = (8,8))
ax[0,0].plot(cluster_info['inflation'], cluster_info['num_clusters'], '.-')
ax[0,0].set_title('Number of clusters')

ax[0,1].plot(cluster_info['inflation'], cluster_info['efficiency'], '.-')
ax[0,1].set_title('Efficiency')
ax[0,1].set_ylim([0,1])

ax[1,0].plot(cluster_info['inflation'], cluster_info['mass_fraction'], '.-')
ax[1,0].set_title('Mass fraction')
ax[1,0].set_ylim([0,1])

ax[1,1].plot(cluster_info['inflation'], cluster_info['area_fraction'], '.-')
ax[1,1].set_title('Area fraction')
ax[1,1].set_ylim([0,1])

ax[2,0].plot(cluster_info['inflation'], cluster_info['percent_difference'], '.-')
ax[2,0].plot(cluster_info['inflation'], [1]*len(cluster_info['inflation']), 'k--', label = '1%')
ax[2,0].set_title('% difference in clustering\nwith next smallest infl. par.')
ax[2,0].set_xlabel('Inflation')
ax[2,0].legend()
ax[2,0].set_yscale('log')

ax[2,1].set_xlabel('Inflation')

fig.suptitle(args.region)

plt.tight_layout()
plt.savefig(args.png_filename, dpi = 300)
plt.savefig(args.pdf_filename)
plt.show()