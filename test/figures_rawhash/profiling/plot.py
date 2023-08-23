import sys
import matplotlib.pyplot as plt 
import matplotlib 
import numpy as np 
import seaborn as sns
import pandas as pd
import csv
from collections import OrderedDict
import pathlib

matplotlib.rc("figure", facecolor="white")
width = 0.5
palette=sns.color_palette('tab10')
colors = ["#6799ce", "#aa4499", "#db8043", "#f0c041", "#4b71bb", "#a3a3a3", "#5FA137"]
f="profiling.csv"

with open(f, 'r') as csvfile:
	csvf = pd.read_csv(csvfile)
	idx=csvf['idx']
	profiling=csvf['Profiling (sec)']
	tools=csvf['Step']
	data=csvf['Data']

fig, ax = plt.subplots(1, 1, figsize=(6,5))

####################################################################################################
##Overlapping Indexing Memory
####################################################################################################

ax.bar(idx[0], profiling[0], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[1], profiling[1], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[2], profiling[2], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[3], profiling[3], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[4], profiling[4], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
# ax.set_yscale("log")

left,right = ax.get_ylim()
ax.set_ylim(70, 100)
ax.margins(x=0.01)
ax.grid(axis='y', linestyle='dotted', linewidth=0.1) 
ax.grid(axis='x', linestyle='--', linewidth=2) 

ax.set_xticks([2,4,6,8])
ax.tick_params(axis="x", which="both", pad=10, direction="in", left=True, labelleft=True) 
ax.tick_params(axis="y", which="both", pad=10, direction="in", rotation=0, top=True) 
ax.set_xticklabels([])


plt.subplots_adjust(top=0.99, bottom=0.01, left=0.005, right=0.995, hspace=0.05, wspace=0.2)
fig.savefig("profiling.pdf")

plt.show()
