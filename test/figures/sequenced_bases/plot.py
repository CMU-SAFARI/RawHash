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
f="sequenced_bases.csv"

with open(f, 'r') as csvfile:
	csvf = pd.read_csv(csvfile)
	idx=csvf['idx']
	sequenced_bases=csvf['Sequenced bases (#)']
	tools=csvf['Tool']
	data=csvf['Data']
	
	#calculating the improvements by X
	mean_seqbaseX=[]
	mean_seqbaseXVal=[]
	toolmean_seqbase = {"UNCALLED": []}

	for i in range(len(data)):
		if tools[i] == 'RawHash':
			rawhashmean_seqbase=sequenced_bases[i]
			mean_seqbaseXVal.append(1)
			mean_seqbaseX.append(' ')
		else:
			mean_seqbaseXVal.append(sequenced_bases[i]/rawhashmean_seqbase)
			mean_seqbaseX.append(r'$%.1f{\times}$' % mean_seqbaseXVal[i])
			toolmean_seqbase[tools[i]].append(sequenced_bases[i]/rawhashmean_seqbase)

# UNCALLED CPU and memory results for latex
print("\\newcommand\\avgmtU{$%.1f\\times$\\xspace}\n\\newcommand\\maxmtU{$%.1f\\times$\\xspace}\n\\newcommand\\minmtU{$%.1f\\times$\\xspace}\n" % (sum(toolmean_seqbase["UNCALLED"])/len(toolmean_seqbase["UNCALLED"]), max(toolmean_seqbase["UNCALLED"]), min(toolmean_seqbase["UNCALLED"])))

fig, ax = plt.subplots(1, 1, figsize=(9,3))

####################################################################################################
##Overlapping Indexing Memory
####################################################################################################

ax.bar(idx[0:2], sequenced_bases[0:2], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[2:4], sequenced_bases[2:4], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[4:6], sequenced_bases[4:6], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[6:8], sequenced_bases[6:8], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
ax.bar(idx[8:10], sequenced_bases[8:10], color=colors, width=0.99, align="center", edgecolor="black", linewidth=1.2)
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

ax.set_xticks([3,6,9,12])
ax.tick_params(axis="x", which="both", pad=10, direction="in", left=True, labelleft=True) 
ax.tick_params(axis="y", which="both", pad=10, direction="in", rotation=0, top=True) 
ax.set_xticklabels([])


plt.subplots_adjust(top=0.99, bottom=0.01, left=0.005, right=0.995, hspace=0.05, wspace=0.2)
fig.savefig("sequenced_bases.pdf")

plt.show()
