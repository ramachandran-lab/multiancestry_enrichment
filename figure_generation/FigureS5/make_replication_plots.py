import pandas as pd 
import numpy as np  
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
from collections import Counter 


fig,ax = plt.subplots(figsize = (9,9), nrows = 2, ncols = 3)

color_dict = {'Anthropometric':'#1b9e77','Blood pressure':'#d95f02','Hematological':'#e7298a','Metabolic':'#7570b3','Kidney':'#66a61e','Other biochemical':'#e6ab02'}

category_dict = pd.read_csv('trait_categories.txt',sep = '\t')
category_dict = {r['Trait']:r['Category'] for i,r in category_dict.iterrows()}
patches = []

trait_order = ['Basophil','Eosinophil','Neutrophil','Monocyte','Lymphocyte','WBC','RBC','MCV','Hemoglobin','Hematocrit','MCH','MCHC','PLC','DBP','SBP','CRP','EGFR','Urate','HBA1C','LDL','HDL','Triglyceride','Cholesterol','BMI','Height']

###########Once the scale files are gnereated you only have to run this to plot
counts = pd.read_csv('scale_counts.txt',sep = '\t',index_col = 'Trait')
proportions = pd.read_csv('scale_proportions.txt',sep = '\t',index_col = 'Trait')
print(counts['Genes'].tolist())

ax[0,0].barh(ind,counts['Variants'].tolist(),color = colors,alpha = 1.00,edgecolor = 'black')
ax[0,1].barh(ind,counts['Clumps'].tolist(),color = colors,alpha = 1.00,edgecolor = 'black')
ax[0,2].barh(ind,counts['Genes'].tolist(),color = colors,alpha = 1.00,edgecolor = 'black')

plt.sca(ax[0,0])
plt.yticks(ind,trait_order)
plt.sca(ax[0,1])
plt.yticks(ind,trait_order)
plt.sca(ax[0,2])
plt.yticks(ind,trait_order)

ax[1,0].barh(ind,proportions['Variants'].tolist(),color = colors,alpha = 1,edgecolor = 'black')
plt.sca(ax[1,0])
plt.yticks(ind,trait_order)

ax[1,1].barh(ind,proportions['Clumps'].tolist(),color = colors,alpha = 1,edgecolor = 'black')
plt.sca(ax[1,1])
plt.yticks(ind,trait_order)

ax[1,2].barh(ind,proportions['Genes'].tolist(),color = colors,edgecolor = 'black')
plt.sca(ax[1,2])
plt.yticks(ind,trait_order)

ax[1,0].set_xlim([0,0.2])
ax[1,1].set_xlim([0,0.2])
ax[1,2].set_xlim([0,0.2])

plt.yticks(ind,trait_order)
plt.tight_layout()


plt.savefig('Figure1.pdf')

