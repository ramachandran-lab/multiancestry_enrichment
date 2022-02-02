import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import scipy.stats
from matplotlib.legend import _get_legend_handles_labels

linetype = [':','-.','--','-']
f = [0.01,0.1,0.25,0.5]

# b_alt = 0.1
# counter = 0
# for f in f:
# 	sigma = np.sqrt(1 - 2*f*(1-f)*b_alt**2) #error sd after SNP effect is accounted for (see next part for explanation)
# 	ns = np.arange(500, 400000, 10) #candidate values for n
# 	ses = sigma/np.sqrt(ns*2*f*(1-f)) #SEs corresponding to each candidate n
# 	q_thresh = scipy.stats.chi2.isf(5e-8, 1) #chi-sqr threshold corresp alpha=5e-8
# 	pwr = scipy.stats.ncx2.sf(q_thresh,1,(b_alt/ses)**2) #power at alpha=5e-8 for VECTOR of SE values
# 	ninety = ns[np.min(np.where(pwr >= 0.90))]
# 	print(ninety)
# 	plt.plot(ns,pwr,linestyle = linetype[counter],color = 'grey')
# 	counter+=1

# plt.tight_layout()
# plt.savefig('power.pdf')
# plt.clf()

ndict = {'hispanic':18377,'african':10032,'south_asian':5716,'oceanian':1915,'aian':604}
colordict = {'african':['#FF7F00','#FFBE7D'],'south_asian':['#E41A1C','#FC9191'],'oceanian':['#4DAF4A','#A3F7A0'],'aian':['#984EA3','#D996FF'],'hispanic':['#FFE800','#FFF1AA']}

fig,ax = plt.subplots(nrows = 3, ncols = 1,sharex = True)
linetype = [':','-.','--','-']

counter = 0
for b_alt in [0.1,0.5,1]:
	ow = 0
	for f in [0.01,0.1,0.25,0.5]:
		sigma = np.sqrt(1 - 2*f*(1-f)*b_alt**2) #error sd after SNP effect is accounted for (see next part for explanation)
		ns = np.arange(100, 300000, 10) #candidate values for n
		ses = sigma/np.sqrt(ns*2*f*(1-f)) #SEs corresponding to each candidate n
		q_thresh = scipy.stats.chi2.isf(5e-8, 1) #chi-sqr threshold corresp alpha=5e-8
		pwr = scipy.stats.ncx2.sf(q_thresh,1,(b_alt/ses)**2) #power at alpha=5e-8 for VECTOR of SE values
		ninety = ns[np.min(np.where(pwr >= 0.90))]
		ax[counter].plot(ns,pwr,linestyle = linetype[ow],color = 'grey',label = f)
		ow +=1
	

	for key,value in ndict.items():
		ax[counter].axvline(value,color = colordict[key][0])

	counter+=1
handles, labels = ax[1].get_legend_handles_labels()
plt.legend(handles = handles, labels = labels,loc = 'center right')

plt.xlim([0,30000])
# plt.tight_layout()
plt.savefig('power.zoom.pdf')