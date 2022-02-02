import pandas as pd 
import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec
from matplotlib import patches
import glob
import os
from itertools import combinations
import seaborn as sns

############################################################
############################################################
#Generate zoomed manhattan and replication matrix
############################################################
############################################################
ancestry_color = {'african':['#FF7F00','#FFBE7D'],'european':['#377EB8','#9AD1FF'],'south_asian':['#E41A1C','#FC9191'],'east_asian':['#FF56E3','#FFB1F2'],'oceanian':['#4DAF4A','#A3F7A0'],'native_american':['#984EA3','#D996FF'],'hispanic':['#FFE800','#FFF1AA'],'shared':['#000000']}
thresh_archive = {'hispanic': 5.829801736604661e-09, 'native_american': 5.8789486840678846e-09, 'oceanian': 7.177547254100102e-09, 'south_asian': 5.217164471109952e-08, 'east_asian': 8.91472177153351e-09, 'african': 4.073609803256053e-09, 'european': 2.5865244147212555e-08}

merged = pd.read_csv('CRP.merged.snps.txt',sep = '\t')
merged = merged.drop(['native_american'],axis = 1)
dim = len(merged.columns)-3
fig,ax = plt.subplots(nrows = dim, ncols = 1,sharex=True)
fig.set_size_inches((8,dim*2))
plt.xlim([0,merged.index.tolist()[-1]])
plt.subplots_adjust(left = 0.1,right = 0.9,bottom=0.1,top = 0.9,wspace = 0.2,hspace = 0.4)

merged = merged.drop('SNP',axis = 1)
merged = merged.astype(float)
merged = merged.sort_values(by = ['#CHROM','POS'])
merged = merged.reset_index(drop = True)
merged = merged.reset_index()
for i,j in enumerate(merged.columns[3:]):
	merged['log'] = np.negative(np.log10(merged[j]))
	grouped = merged.groupby(('#CHROM'))

	for num, (name, group) in enumerate(grouped):
		ax[i].scatter(group['index'],group['log'],color=ancestry_color[j][num % 2],s = 5)
		ax[i].spines['right'].set_visible(False)
		ax[i].spines['top'].set_visible(False)
		ax[i].yaxis.set_ticks_position('left')
		ax[i].xaxis.set_ticks_position('bottom')
		ax[i].axhline(y=np.negative(np.log10(thresh_archive[j])),color = 'grey',linestyle = '--')
		# ax[i].set_ylim(ymin=3)
xlabels = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,19,21]
xpos = []
for i in xlabels:
	temp = merged[merged['#CHROM'] == i]
	max_pos = np.max(temp['index'])
	min_pos = np.min(temp['index'])
	xpos.append((max_pos + min_pos)/2)
plt.sca(ax[-1])
plt.xticks(xpos,xlabels)
# plt.suptitle(trait)
merged = merged.drop(['log'],axis = 1)
plt.tight_layout()
plt.savefig('CRP.test.png')
plt.clf()
fig,ax = plt.subplots(nrows = dim, ncols = 1,sharex=True)

new = merged[(merged['POS'] >= 159661585) & (merged['POS'] <= 159703379) & (merged['#CHROM'] == 1)]
# plt.xlim([0,new.index.tolist()[-1]])
fig.set_size_inches((3,dim*2))
plt.subplots_adjust(left = 0.1,right = 0.9,bottom=0.1,top = 0.9,wspace = 0.2,hspace = 0.4)
new = new.reset_index(drop = True)

for i,j in enumerate(new.columns[3:]):
	new['ancestry'] = new[j] < thresh_archive[j]
	others = ['african','european','south_asian','east_asian','hispanic','oceanian']
	others.remove(j)
	new['replicated'] = new[others].isnull().astype(int).sum(axis = 1)
	new['test_replication'] = new['replicated'] <= 5
	new['label'] = new[['ancestry','test_replication']].all(axis =1).astype(int)
	color_map = {1:ancestry_color[j][0],0:ancestry_color[j][0]}
	shape_map = {1:'*',0:'o'}
	new['color'] = new['label'].map(color_map)
	new['shape'] = new['label'].map(shape_map)
	new['log'] = np.negative(np.log10(new[j]))
	grouped = new.groupby(('#CHROM'))
	for num, (name, group) in enumerate(grouped):
		# shared = group[group['shape'] == '*']
		ax[i].axhline(y=np.negative(np.log10(thresh_archive[j])),color = 'grey',linestyle = '--')
		shared = group[group['POS'] == 159684665.0]
		ax[i].scatter(group['POS'],group['log'],color=group['color'],marker = 'o',s = 5)
		ax[i].scatter(shared['POS'],shared['log'],color='black',marker = '*',s = 100)
		ax[i].spines['right'].set_visible(False)
		ax[i].spines['top'].set_visible(False)
		ax[i].yaxis.set_ticks_position('left')
		ax[i].xaxis.set_ticks_position('bottom')
		ax[i].set_xlim(group.iloc[0]['POS'],group.iloc[-1]['POS'])
ax[0].axvline(x=159682078,color = 'black')
ax[0].axvline(x=159684379)
		# ax[i].set_ylim(ymin=3)
print(new.iloc[0]['POS'])
print(new.iloc[-1]['POS'])
plt.xticks([159682078],[19])
plt.sca(ax[0])
plt.yticks([0,20,40,60,80])
plt.sca(ax[1])
plt.yticks([0,100,200])
plt.sca(ax[2])
plt.yticks([4,6,8])


# merged = merged.drop(['#CHROM','POS','index','log'],axis = 1)
plt.tight_layout()
plt.savefig('CRP.zoom.pdf')
plt.clf()

#############Heatmaps with imputed data
# merged = merged[(merged['POS'] >= 159682078) & (merged['POS'] >= 159684379) & (merged['#CHROM'] == 1)]
heat = pd.DataFrame(np.zeros((6,6)),index = ['african','european','south_asian','east_asian','hispanic','oceanian'],columns = ['african','european','south_asian','east_asian','hispanic','oceanian'])
total_replicates = set()
for i in combinations(['african','european','south_asian','east_asian','hispanic','oceanian'],2):
	temp = merged[(merged[i[0]] <= thresh_archive[i[0]]) & (merged[i[1]] <= thresh_archive[i[1]])]
	replicates = temp.shape[0]
	# heat.loc[i[0],i[1]] = replicates
	heat.loc[i[1],i[0]] = replicates
	for snp in temp['SNP'].tolist():
		total_replicates.add(snp)
print(len(total_replicates))
for i in ['african','european','south_asian','east_asian','hispanic','oceanian']:
	temp = merged[merged[i] <= thresh_archive[i]]
	replicates = temp.shape[0]
	heat.loc[i,i] = replicates
print(heat)
plt.figure()
mask = np.triu(np.ones_like(heat, dtype=np.bool),k=1)
sns.heatmap(heat,cmap = ['white'],cbar = False,annot = True,mask = mask)
plt.savefig('CRP.heatmap.numbers.imputed.pdf')
plt.clf()

for i in combinations(['african','european','south_asian','east_asian','hispanic','oceanian'],2):
	denom = np.min([heat.loc[i[0],i[0]],heat.loc[i[1],i[1]]])
	heat.loc[i[1],i[0]] = heat.loc[i[1],i[0]]/denom

for i in ['african','european','south_asian','east_asian','hispanic','oceanian']:
	heat.loc[i,i] = 0
mask = np.triu(np.ones_like(heat, dtype=np.bool))
plt.figure()
sns.heatmap(heat,cmap = 'Reds',cbar = False,annot = False,mask = mask)
plt.savefig('CRP.heatmap.color.imputed.pdf')
plt.clf()

#########heatmaps without imputed data
# merged = merged[(merged['POS'] >= 159682078) & (merged['POS'] >= 159684379) & (merged['#CHROM'] == 1)]
first_genotyped = pd.read_csv('CRP.african.genotyped.txt',header = None)
first_genotyped = first_genotyped[first_genotyped[0] != '.']
temp = merged[['SNP','african']].dropna()
hold = temp[temp['SNP'].isin(first_genotyped[0].tolist())]
total_replicates = set()
for ancestry in ['european','south_asian','east_asian','hispanic','oceanian']:
	first_genotyped = pd.read_csv('CRP.' + ancestry + '.genotyped.txt',header = None)
	first_genotyped = first_genotyped[first_genotyped[0] != '.']
	temp = merged[['SNP',ancestry]].dropna()
	temp = temp[temp['SNP'].isin(first_genotyped[0].tolist())]
	hold = hold.merge(temp,on = 'SNP',how = 'outer')


heat = pd.DataFrame(np.zeros((6,6)),index = ['african','european','south_asian','east_asian','hispanic','oceanian'],columns = ['african','european','south_asian','east_asian','hispanic','oceanian'])
for i in combinations(['african','european','south_asian','east_asian','hispanic','oceanian'],2):
	temp = hold[(hold[i[0]] <= thresh_archive[i[0]]) & (hold[i[1]] <= thresh_archive[i[1]])]
	replicates = temp.shape[0]
	# heat.loc[i[0],i[1]] = replicates
	heat.loc[i[1],i[0]] = replicates
	for snp in temp['SNP'].tolist():
		total_replicates.add(snp)
print(len(total_replicates))

for i in ['african','european','south_asian','east_asian','hispanic','oceanian']:
	temp = hold[hold[i] <= thresh_archive[i]]
	replicates = temp.shape[0]
	heat.loc[i,i] = replicates

plt.figure()
mask = np.triu(np.ones_like(heat, dtype=np.bool),k=1)
sns.heatmap(heat,cmap = ['white'],cbar = False,annot = True,mask = mask)
plt.savefig('CRP.heatmap.numbers.genotyped.pdf')
plt.clf()

for i in combinations(['african','european','south_asian','east_asian','hispanic','oceanian'],2):
	denom = np.min([heat.loc[i[0],i[0]],heat.loc[i[1],i[1]]])
	heat.loc[i[1],i[0]] = heat.loc[i[1],i[0]]/denom

for i in ['african','european','south_asian','east_asian','hispanic','oceanian']:
	heat.loc[i,i] = 0
mask = np.triu(np.ones_like(heat, dtype=np.bool))
plt.figure()
sns.heatmap(heat,cmap = 'Reds',cbar = True,annot = False,mask = mask)
plt.savefig('CRP.heatmap.color.genotyped.pdf')
plt.clf()




