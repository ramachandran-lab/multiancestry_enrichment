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

############################################################
############################################################
#Read in and combine the different vartiant level outputs
############################################################
############################################################
ancestry_color = {'african':['#FF7F00','#FFBE7D'],'european':['#377EB8','#9AD1FF'],'south_asian':['#E41A1C','#FC9191'],'east_asian':['#FF56E3','#FFB1F2'],'oceanian':['#4DAF4A','#A3F7A0'],'native_american':['#984EA3','#D996FF'],'hispanic':['#FFE800','#FFF1AA'],'shared':['#000000']}
traits = ['CRP']
for trait in traits:
	thresh_archive = {}
	if os.path.exists(trait + '.merged.snps.txt'):
		merged = pd.read_csv(trait + '.merged.snps.txt',sep = '\t')
		temp = [i.strip().split('\t') for i in open(trait + '.thresholds.txt','r')]
		thresh_archive = {}
		for i in temp:
			thresh_archive[i[0]] = float(i[1])
	else:
		files = glob.glob('*/' + trait + '*.snp.stats.txt')
		if 'african/' + trait + '.page.snp.stats.txt' in files:
			files.remove('african/' + trait + '.snp.stats.txt')
			page_data = True
		else:
			page_data = False

		if page_data:
			african = pd.read_csv('african/' + trait + '.page.snp.stats.txt',sep = '\t')[['CHROM','POS','RS_ID','PVALUE']]
			african.columns = ['#CHROM','POS','SNP','african']
			thresh = 0.05/african.shape[0]
			african = african[african['african'] <= 0.001]
			african['id'] = african['#CHROM'].astype(str) + '_' + african['POS'].astype(str)
		else:
			african = pd.read_csv('african/' + trait + '.snp.stats.txt',sep = '\t')[['#CHROM','POS','ID','P']]
			african.columns = ['#CHROM','POS','SNP','african']
			thresh = 0.05/african.shape[0]

			african = african[african['african'] <= 0.001]
			african['id'] = african['#CHROM'].astype(str) + '_' + african['POS'].astype(str)
		thresh_archive['african'] = thresh

		european = pd.read_csv('european/' + trait + '.snp.stats.txt',sep = '\t')[['#CHROM','POS','ID','P']]
		european.columns = ['#CHROM','POS','SNP','european']
		thresh = 0.05/european.shape[0]
		thresh_archive['european'] = thresh

		european = european[european['european'] <= 0.001]
		european['id'] = european['#CHROM'].astype(str) + '_' + european['POS'].astype(str)
		merged = african.merge(european,on = 'SNP',how = 'outer')
		merged.loc[merged["#CHROM_x"].isnull(),'#CHROM_x'] = merged['#CHROM_y']
		merged.loc[merged['POS_x'].isnull(),'POS_x'] = merged['POS_y']
		merged = merged[['#CHROM_x','POS_x','SNP','african','european']]
		merged.columns = ['#CHROM','POS','SNP','african','european']

		south_asian = pd.read_csv('south_asian/' + trait + '.snp.stats.txt',sep = '\t')[['#CHROM','POS','ID','P']]
		south_asian.columns = ['#CHROM','POS','SNP','south_asian']
		thresh = 0.05/south_asian.shape[0]
		thresh_archive['south_asian'] = thresh
		south_asian = south_asian[south_asian['south_asian'] <= 0.001]
		south_asian['id'] = south_asian['#CHROM'].astype(str) + '_' + south_asian['POS'].astype(str)
		merged = merged.merge(south_asian,on = 'SNP',how = 'outer')
		merged.loc[merged["#CHROM_x"].isnull(),'#CHROM_x'] = merged['#CHROM_y']
		merged.loc[merged['POS_x'].isnull(),'POS_x'] = merged['POS_y']
		merged = merged[['#CHROM_x','POS_x','SNP','african','european','south_asian']]
		merged.columns = ['#CHROM','POS','SNP','african','european','south_asian']


		east_asian = pd.read_csv('east_asian/' + trait + '.snp.stats.txt',sep = '\t')[['#CHROM','POS','SNP','P']]
		east_asian.columns = ['#CHROM','POS','SNP','east_asian']
		thresh = 0.05/east_asian.shape[0]
		thresh_archive['east_asian'] = thresh
		east_asian = east_asian[east_asian['east_asian'] <= 0.001]
		east_asian['id'] = east_asian['#CHROM'].astype(str) + '_' + east_asian['POS'].astype(str)
		merged = merged.merge(east_asian,on = 'SNP',how = 'outer')
		merged.loc[merged["#CHROM_x"].isnull(),'#CHROM_x'] = merged['#CHROM_y']
		merged.loc[merged['POS_x'].isnull(),'POS_x'] = merged['POS_y']
		merged = merged[['#CHROM_x','POS_x','SNP','african','european','south_asian','east_asian']]
		merged.columns = ['#CHROM','POS','SNP','african','european','south_asian','east_asian']


		if 'oceanian/' + trait + '.page.snp.stats.txt' in files:
			oceanian = pd.read_csv('oceanian/' + trait + '.page.snp.stats.txt',sep = '\t')[['CHROM','POS','RS_ID','PVALUE']]
			oceanian.columns = ['#CHROM','POS','SNP','oceanian']
			thresh = 0.05/oceanian.shape[0]
			thresh_archive['oceanian'] = thresh
			oceanian = oceanian[oceanian['oceanian'] <= 0.001]
			oceanian['id'] = oceanian['#CHROM'].astype(str) + '_' + oceanian['POS'].astype(str)
			stat_oceanian = True
			merged = merged.merge(oceanian,on = 'SNP',how = 'outer')
			merged.loc[merged["#CHROM_x"].isnull(),'#CHROM_x'] = merged['#CHROM_y']
			merged.loc[merged['POS_x'].isnull(),'POS_x'] = merged['POS_y']
			merged = merged[['#CHROM_x','POS_x','SNP','african','european','south_asian','east_asian','oceanian']]
			merged.columns = ['#CHROM','POS','SNP','african','european','south_asian','east_asian','oceanian']
		else:
			stat_oceanian = False

		if 'native_american/' + trait + '.page.snp.stats.txt' in files:
			native_american = pd.read_csv('native_american/' + trait + '.page.snp.stats.txt',sep = '\t')[['CHROM','POS','RS_ID','PVALUE']]
			native_american.columns = ['#CHROM','POS','SNP','native_american']
			thresh = 0.05/native_american.shape[0]
			thresh_archive['native_american'] = thresh
			native_american = native_american[native_american['native_american'] <= 0.001]
			native_american['id'] = native_american['#CHROM'].astype(str) + '_' + native_american['POS'].astype(str)
			stat_native_american = True
			if stat_oceanian == True:
				merged = merged.merge(native_american,on = 'SNP',how = 'outer')
				merged.loc[merged["#CHROM_x"].isnull(),'#CHROM_x'] = merged['#CHROM_y']
				merged.loc[merged['POS_x'].isnull(),'POS_x'] = merged['POS_y']
				merged = merged[['#CHROM_x','POS_x','SNP','african','european','south_asian','east_asian','oceanian','native_american']]
				merged.columns = ['#CHROM','POS','SNP','african','european','south_asian','east_asian','oceanian','native_american']

			else:
				merged = merged.merge(native_american,on = 'SNP',how = 'outer')
				merged.loc[merged["#CHROM_x"].isnull(),'#CHROM_x'] = merged['#CHROM_y']
				merged.loc[merged['POS_x'].isnull(),'POS_x'] = merged['POS_y']
				merged = merged[['#CHROM_x','POS_x','SNP','african','european','south_asian','east_asian','native_american']]
				merged.columns = ['#CHROM','POS','SNP','african','european','south_asian','east_asian','native_american']


		else:
			stat_native_american = False
		if 'hispanic/' + trait + '.page.snp.stats.txt' in files:
			hispanic = pd.read_csv('hispanic/' + trait + '.page.snp.stats.txt',sep = '\t')[['CHROM','POS','RS_ID','PVALUE']]
			hispanic.columns = ['#CHROM','POS','SNP','hispanic']
			thresh = 0.05/hispanic.shape[0]
			thresh_archive['hispanic'] = thresh
			hispanic = hispanic[hispanic['hispanic'] <= 0.001]
			hispanic['id'] = hispanic['#CHROM'].astype(str) + '_' + hispanic['POS'].astype(str)
			
			if stat_oceanian and stat_native_american:
				merged = merged.merge(hispanic,on = 'SNP',how = 'outer')
				merged.loc[merged["#CHROM_x"].isnull(),'#CHROM_x'] = merged['#CHROM_y']
				merged.loc[merged['POS_x'].isnull(),'POS_x'] = merged['POS_y']
				merged = merged[['#CHROM_x','POS_x','SNP','african','european','south_asian','east_asian','oceanian','native_american','hispanic']]
				merged.columns = ['#CHROM','POS','SNP','african','european','south_asian','east_asian','oceanian','native_american','hispanic']

			elif stat_oceanian and not stat_native_american:
				merged = merged.merge(hispanic,on = 'SNP',how = 'outer')
				merged.loc[merged["#CHROM_x"].isnull(),'#CHROM_x'] = merged['#CHROM_y']
				merged.loc[merged['POS_x'].isnull(),'POS_x'] = merged['POS_y']
				merged = merged[['#CHROM_x','POS_x','SNP','african','european','south_asian','east_asian','oceanian','hispanic']]
				merged.columns = ['#CHROM','POS','SNP','african','european','south_asian','east_asian','oceanian','hispanic']
				
			elif not stat_oceanian and stat_native_american:
				merged = merged.merge(hispanic,on = 'SNP',how = 'outer')
				merged.loc[merged["#CHROM_x"].isnull(),'#CHROM_x'] = merged['#CHROM_y']
				merged.loc[merged['POS_x'].isnull(),'POS_x'] = merged['POS_y']
				merged = merged[['#CHROM_x','POS_x','SNP','african','european','south_asian','east_asian','native_american','hispanic']]
				merged.columns = ['#CHROM','POS','SNP','african','european','south_asian','east_asian','native_american','hispanic']

			elif not stat_oceanian and not stat_native_american:
				merged = merged.merge(hispanic,on = 'SNP',how = 'outer')
				merged.loc[merged["#CHROM_x"].isnull(),'#CHROM_x'] = merged['#CHROM_y']
				merged.loc[merged['POS_x'].isnull(),'POS_x'] = merged['POS_y']
				merged = merged[['#CHROM_x','POS_x','SNP','african','european','south_asian','east_asian','hispanic']]
				merged.columns = ['#CHROM','POS','SNP','african','european','south_asian','east_asian','hispanic']

		merged.to_csv(trait + '.merged.snps.txt',sep = '\t',index = False)
		newfile = open(trait + '.thresholds.txt','w')
		for key,value in thresh_archive.items():
			newfile.write(str(key) + '\t' + str(value) + '\n')
	dim = len(merged.columns)-3
	fig,ax = plt.subplots(nrows = dim, ncols = 1,sharex=True)
	fig.set_size_inches((8,dim*2))
	plt.xlim([0,merged.index.tolist()[-1]])
	plt.subplots_adjust(left = 0.1,right = 0.9,bottom=0.1,top = 0.9,wspace = 0.2,hspace = 0.4)
	
	upset_save = merged
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
			ax[i].axhline(y=np.negative(np.log10(thresh_archive[j])),color = 'black',linestyle = '--')
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
	merged = merged.drop(['#CHROM','POS','index','log'],axis = 1)
	plt.tight_layout()
	plt.savefig(trait +'.png')
	plt.clf()

	merged = pd.read_csv(trait + '.merged.snps.txt',sep = '\t')
	# merged = merged.drop(['native_american'],axis = 1)
	first = merged[merged['african'] <= thresh_archive['african']]

	first = first[['SNP','#CHROM','POS','african']]
	# print(first)
	keep_columns = ['SNP','#CHROM','POS','african']
	for i in merged.columns[4:]:
		temp = merged[merged[i] <= thresh_archive[i]]
		temp = temp[['SNP','#CHROM','POS',i]]
		first = first.merge(temp,on = 'SNP', how = 'outer')
		keep_columns.append(j)
		first.loc[first['POS_x'].isnull(),'POS_x'] = first['POS_y']
		first.loc[first["#CHROM_x"].isnull(),'#CHROM_x'] = first['#CHROM_y']
		first = first.rename(columns = {'#CHROM_x':'#CHROM','POS_x':'POS'})
		first = first.drop(['#CHROM_y','POS_y'],axis = 1)


	
	first = first.set_index('SNP')
	for j in first.columns[2:]:
		first[j] = first[j] > 0
		first[j] = first[j].astype(int)
	# print(first)
	first.to_csv(trait + '.upset.input.csv')
	# first['occurences'] = first.sum(axis=1)


