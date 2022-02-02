import pandas as pd 
import numpy as np  
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
from collections import Counter 



color_dict = {'Anthropometric':'#1b9e77','Blood pressure':'#d95f02','Hematological':'#e7298a','Metabolic':'#7570b3','Kidney':'#66a61e','Other biochemical':'#e6ab02'}
marker_dict = {'Anthropometric':'*','Blood pressure':'p','Hematological':'^','Metabolic':'s','Kidney':'d','Other biochemical':'8'}
category_dict = pd.read_csv('trait_categories.txt',sep = '\t')
category_dict = {r['Trait']:r['Category'] for i,r in category_dict.iterrows()}


trait_order = ['Basophil','Eosinophil','Neutrophil','Monocyte','Lymphocyte','WBC','RBC','MCV','Hemoglobin','Hematocrit','MCH','MCHC','PLC','DBP','SBP','CRP','EGFR','Urate','HBA1C','LDL','HDL','Triglyceride','Cholesterol','BMI','Height']

###########Once the scale files are gnereated you only have to run this to plot
#1a

proportions = pd.read_csv('scale_proportions.txt',sep = '\t',index_col = 'Traits')
fig,ax = plt.subplots(figsize = (4,4), nrows = 1, ncols = 1)
proportions = proportions.reset_index()
proportions['color'] = proportions['Traits'].map(category_dict)
proportions['color'] = proportions['color'].map(color_dict)
proportions['shape'] = proportions['Traits'].map(category_dict)
proportions['shape'] = proportions['shape'].map(marker_dict)

ax.plot([0,1],[0,1], transform=ax.transAxes, color = 'black', linestyle = '--')
proportions = proportions.set_index('Traits')


for x,y in proportions.iterrows():
	ax.scatter(y['Variants'],y['Genes'],color = y['color'],marker=y['shape'])

plt.plot(np.mean(y['Variants']),0,marker = '*', color = 'black', markersize = 14)
plt.plot(0,np.mean(y['Genes']),marker = '*', color = 'black', markersize = 14)

ax.set_ylim([0,0.2])
ax.set_xlim([0,0.2])
ax.set_xticks([0,0.05,0.10,0.15,0.20])
ax.set_yticks([0,0.05,0.10,0.15,0.20])

plt.savefig('Figure1a.pdf')

#1b
proportions = pd.read_csv('scale_proportions_alt.txt',sep = '\t',index_col = 'Traits')
fig,ax = plt.subplots(figsize = (4,4), nrows = 1, ncols = 1)
proportions = proportions.reset_index()
proportions['color'] = proportions['Traits'].map(category_dict)
proportions['color'] = proportions['color'].map(color_dict)
proportions['shape'] = proportions['Traits'].map(category_dict)
proportions['shape'] = proportions['shape'].map(marker_dict)

ax.plot([0,1],[0,1], transform=ax.transAxes, color = 'black', linestyle = '--')


for x,y in proportions.iterrows():
	ax.scatter(y['Variants'],y['Genes'],color = y['color'],marker=y['shape'])

plt.plot(np.mean(y['Variants']),0,marker = '*', color = 'black', markersize = 14)
plt.plot(0,np.mean(y['Genes']),marker = '*', color = 'black', markersize = 14)

ax.set_ylim([0,0.2])
ax.set_xlim([0,0.2])
ax.set_xticks([0,0.05,0.10,0.15,0.20])
ax.set_yticks([0,0.05,0.10,0.15,0.20])

plt.savefig('Figure1b.pdf')

#1c
proportions = pd.read_csv('scale_proportions_alt_nothresh.txt',sep = '\t',index_col = 'Traits')
fig,ax = plt.subplots(figsize = (4,4), nrows = 1, ncols = 1)
proportions = proportions.reset_index()
proportions['color'] = proportions['Traits'].map(category_dict)
proportions['color'] = proportions['color'].map(color_dict)
proportions['shape'] = proportions['Traits'].map(category_dict)
proportions['shape'] = proportions['shape'].map(marker_dict)

ax.plot([0,1],[0,1], transform=ax.transAxes, color = 'black', linestyle = '--')


for x,y in proportions.iterrows():
	ax.scatter(y['Variants'],y['Genes'],color = y['color'],marker=y['shape'])

plt.plot(np.mean(y['Variants']),0,marker = '*', color = 'black', markersize = 14)
plt.plot(0,np.mean(y['Genes']),marker = '*', color = 'black', markersize = 14)

ax.set_ylim([0,0.2])
ax.set_xlim([0,0.2])
ax.set_xticks([0,0.05,0.10,0.15,0.20])
ax.set_yticks([0,0.05,0.10,0.15,0.20])

plt.savefig('Figure1c.pdf')

#1d
proportions = pd.read_csv('scale_proportions_nominal.txt',sep = '\t',index_col = 'Traits')
fig,ax = plt.subplots(figsize = (4,4), nrows = 1, ncols = 1)
proportions = proportions.reset_index()
proportions['color'] = proportions['Traits'].map(category_dict)
proportions['color'] = proportions['color'].map(color_dict)
proportions['shape'] = proportions['Traits'].map(category_dict)
proportions['shape'] = proportions['shape'].map(marker_dict)

proportions[['Traits','Variants','Genes']].to_csv('scale_proportions_nominal.txt',sep = '\t',index = False)
ax.plot([0,1],[0,1], transform=ax.transAxes, color = 'black', linestyle = '--')


for x,y in proportions.iterrows():
	ax.scatter(y['Variants'],y['Genes'],color = y['color'],marker=y['shape'])

plt.plot(np.mean(y['Variants']),0,marker = '*', color = 'black', markersize = 14)
plt.plot(0,np.mean(y['Genes']),marker = '*', color = 'black', markersize = 14)

ax.set_ylim([0,0.1])
ax.set_xlim([0,1])
# ax.set_xticks([0,0.05,0.10,0.15,0.20])
# ax.set_yticks([0,0.05,0.10,0.15,0.20])

plt.savefig('Figure1d.pdf')

