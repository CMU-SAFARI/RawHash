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
colors = ["#5FA137", "#f0c041", "#6799ce", "#aa4499", "#db8043", "#a3a3a3", "#4b71bb"]
f="throughput.csv"

with open(f, 'r') as csvfile:
	csvf = pd.read_csv(csvfile)
	idx=csvf['idx']
	throughput=csvf['Throughput (bp/sec)']
	tools=csvf['Tool']
	data=csvf['Data']
	
	#calculating the improvements by X
	throughputX=[]
	throughputXVal=[]
	toolThroughput = {"UNCALLED": [], "Sigmap": [], "RawHash": [], "RawHash2-Minimizer": []}

	for i in range(len(data)):
		if tools[i] == 'RawHash2':
			rawhashThroughput=throughput[i]
			# throughputXVal.append(1)
			# throughputX.append(' ')
			throughputXVal.append(rawhashThroughput/450)
			throughputX.append(r'$%.1f{\times}$' % throughputXVal[i])
		else:
			# throughputXVal.append(rawhashThroughput/throughput[i])
			throughputXVal.append(throughput[i]/450)
			throughputX.append(r'$%.1f{\times}$' % throughputXVal[i])
			toolThroughput[tools[i]].append(rawhashThroughput/throughput[i])

# UNCALLED CPU and memory results for latex
print("\\newcommand\\avgthrU{$%.1f\\times$\\xspace}\n\\newcommand\\maxthrU{$%.1f\\times$\\xspace}\n\\newcommand\\minthrU{$%.1f\\times$\\xspace}\n" % (sum(toolThroughput["UNCALLED"])/len(toolThroughput["UNCALLED"]), max(toolThroughput["UNCALLED"]), min(toolThroughput["UNCALLED"])))

# Sigmap CPU and memory results for latex
print("\\newcommand\\avgthrS{$%.1f\\times$\\xspace}\n\\newcommand\\maxthrS{$%.1f\\times$\\xspace}\n\\newcommand\\minthrS{$%.1f\\times$\\xspace}\n" % (sum(toolThroughput["Sigmap"])/len(toolThroughput["Sigmap"]), max(toolThroughput["Sigmap"]), min(toolThroughput["Sigmap"])))

fig, ax = plt.subplots(1, 1, figsize=(9,3))

####################################################################################################
##Overlapping Indexing Memory
####################################################################################################

ax.bar(1, 450, color=colors[5], width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[0:5]+2, throughput[0:5], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[5:10]+2, throughput[5:10], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[10:15]+2, throughput[10:15], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[15:20]+2, throughput[15:20], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[20:25]+2, throughput[20:25], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[25:30]+2, throughput[25:30], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[30:35]+2, throughput[30:35], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.set_yscale("log")

cidx=0
bidx=0
for c in ax.containers:
	if len(c.datavalues) > 1:
		ax.bar_label(c, labels=throughputX[bidx:bidx+len(c.datavalues)], fmt='%gs', color='white', fontsize=12, weight='bold', rotation='vertical', label_type='center', padding=-28)
		bidx += len(c.datavalues)
	else:
		ax.bar_label(c, labels=['450 bp/sec'], fmt='%gs', color='white', fontsize=12, rotation='vertical', label_type='center', padding=-33)

left,right = ax.get_ylim()
ax.set_ylim(1, right+2*right)
ax.margins(x=0.01)
ax.grid(axis='y', linestyle='dotted', linewidth=0.1) 
ax.grid(axis='x', linestyle='--', linewidth=2) 

ax.set_xticks([2,8,14,20,26,32,38])
ax.tick_params(axis="x", which="both", pad=10, direction="in", left=True, labelleft=True) 
ax.tick_params(axis="y", which="both", pad=10, direction="in", rotation=0, top=True) 
ax.set_xticklabels([])


plt.subplots_adjust(top=0.99, bottom=0.01, left=0.005, right=0.995, hspace=0.05, wspace=0.2)
fig.savefig("throughput.pdf")

plt.show()
