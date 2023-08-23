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
f="mean_time.csv"

with open(f, 'r') as csvfile:
	csvf = pd.read_csv(csvfile)
	idx=csvf['idx']
	mean_time=csvf['Mean time spent per read (ms)']
	tools=csvf['Tool']
	data=csvf['Data']
	
	#calculating the improvements by X
	mean_timeX=[]
	mean_timeXVal=[]
	toolmean_time = {"UNCALLED": [], "Sigmap": [], "RawHash": [], "RawHash2-Minimizer": []}

	for i in range(len(data)):
		if tools[i] == 'RawHash2':
			rawhashmean_time=mean_time[i]
			mean_timeXVal.append(1)
			mean_timeX.append(' ')
		else:
			mean_timeXVal.append(mean_time[i]/rawhashmean_time)
			mean_timeX.append(r'$%.1f{\times}$' % mean_timeXVal[i])
			toolmean_time[tools[i]].append(mean_time[i]/rawhashmean_time)

# UNCALLED CPU and memory results for latex
print("\\newcommand\\avgmtU{$%.1f\\times$\\xspace}\n\\newcommand\\maxmtU{$%.1f\\times$\\xspace}\n\\newcommand\\minmtU{$%.1f\\times$\\xspace}\n" % (sum(toolmean_time["UNCALLED"])/len(toolmean_time["UNCALLED"]), max(toolmean_time["UNCALLED"]), min(toolmean_time["UNCALLED"])))

# Sigmap CPU and memory results for latex
print("\\newcommand\\avgmtS{$%.1f\\times$\\xspace}\n\\newcommand\\maxmtS{$%.1f\\times$\\xspace}\n\\newcommand\\minmtS{$%.1f\\times$\\xspace}\n" % (sum(toolmean_time["Sigmap"])/len(toolmean_time["Sigmap"]), max(toolmean_time["Sigmap"]), min(toolmean_time["Sigmap"])))

fig, ax = plt.subplots(1, 1, figsize=(9,3))

####################################################################################################
##Overlapping Indexing Memory
####################################################################################################

ax.bar(idx[0:5], mean_time[0:5], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[5:10], mean_time[5:10], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[10:15], mean_time[10:15], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[15:20], mean_time[15:20], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[20:25], mean_time[20:25], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[25:30], mean_time[25:30], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[30:35], mean_time[30:35], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.set_yscale("log")

cidx=0
bidx=0
for c in ax.containers:
	ax.bar_label(c, labels=mean_timeX[bidx:bidx+len(c.datavalues)], fmt='%gs', color='white', fontsize=12, weight='bold', rotation='vertical', label_type='center', padding=-24)
	bidx += len(c.datavalues)

left,right = ax.get_ylim()
ax.set_ylim(0.01, right+2*right)
ax.margins(x=0.01)
ax.grid(axis='y', linestyle='dotted', linewidth=0.1) 
ax.grid(axis='x', linestyle='--', linewidth=2) 

ax.set_xticks([6,12,18,24,30,36])
ax.tick_params(axis="x", which="both", pad=10, direction="in", left=True, labelleft=True) 
ax.tick_params(axis="y", which="both", pad=10, direction="in", rotation=0, top=True) 
ax.set_xticklabels([])


plt.subplots_adjust(top=0.99, bottom=0.01, left=0.005, right=0.995, hspace=0.05, wspace=0.2)
fig.savefig("mean_time.pdf")

plt.show()
