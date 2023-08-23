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
colors = ["#5FA137", "#f0c041", "#6799ce", "#db8043", "#aa4499", "#4b71bb", "#a3a3a3"]
f="sequenced_chunks.csv"

with open(f, 'r') as csvfile:
	csvf = pd.read_csv(csvfile)
	idx=csvf['idx']
	sequenced_chunks=csvf['Sequenced chunks (#)']
	tools=csvf['Tool']
	data=csvf['Data']
	
	#calculating the improvements by X
	mean_seqbaseX=[]
	mean_seqbaseXVal=[]
	toolmean_seqbase = {"Sigmap": [], "RawHash": [], "RawHash2-Minimizer": []}

	for i in range(len(data)):
		if tools[i] == 'RawHash2':
			rawhashmean_seqbase=sequenced_chunks[i]
			mean_seqbaseXVal.append(1)
			mean_seqbaseX.append(' ')
		else:
			mean_seqbaseXVal.append(sequenced_chunks[i]/rawhashmean_seqbase)
			mean_seqbaseX.append(r'$%.1f{\times}$' % mean_seqbaseXVal[i])
			toolmean_seqbase[tools[i]].append(sequenced_chunks[i]/rawhashmean_seqbase)

# Sigmap CPU and memory results for latex
print("\\newcommand\\avgmtU{$%.1f\\times$\\xspace}\n\\newcommand\\maxmtU{$%.1f\\times$\\xspace}\n\\newcommand\\minmtU{$%.1f\\times$\\xspace}\n" % (sum(toolmean_seqbase["Sigmap"])/len(toolmean_seqbase["Sigmap"]), max(toolmean_seqbase["Sigmap"]), min(toolmean_seqbase["Sigmap"])))

fig, ax = plt.subplots(1, 1, figsize=(9,3))

####################################################################################################
##Overlapping Indexing Memory
####################################################################################################

ax.bar(idx[0:4], sequenced_chunks[0:4], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[4:8], sequenced_chunks[4:8], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[8:12], sequenced_chunks[8:12], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[12:16], sequenced_chunks[12:16], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[16:20], sequenced_chunks[16:20], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
# ax.set_yscale("log")

cidx=0
bidx=0
for c in ax.containers:
	ax.bar_label(c, labels=mean_seqbaseX[bidx:bidx+len(c.datavalues)], fmt='%gs', color='black', fontsize=18, weight='bold', rotation='vertical', label_type='edge', padding=0)
	bidx += len(c.datavalues)

left,right = ax.get_ylim()
ax.set_ylim(0.01, right+0.2*right)
ax.margins(x=0.01)
ax.grid(axis='y', linestyle='dotted', linewidth=0.1) 
ax.grid(axis='x', linestyle='--', linewidth=2) 

ax.set_xticks([5,10,15,20])
ax.tick_params(axis="x", which="both", pad=10, direction="in", left=True, labelleft=True) 
ax.tick_params(axis="y", which="both", pad=10, direction="in", rotation=0, top=True) 
ax.set_xticklabels([])


plt.subplots_adjust(top=0.99, bottom=0.01, left=0.005, right=0.995, hspace=0.05, wspace=0.2)
fig.savefig("sequenced_chunks.pdf")

plt.show()
