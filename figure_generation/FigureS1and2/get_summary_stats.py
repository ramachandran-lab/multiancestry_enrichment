import pandas as pd 
import numpy as np 
import sys
import matplotlib.pyplot as plt 

new = pd.DataFrame(np.zeros((len([j for j in range(2012,2021)]),4)), index = [k for k in range(2012,2021)], columns = ['on_both','PMID_no_UKB','UKB_no_pubmed','Total_pubs'])

for i in range(2012,2021):
	df = pd.read_csv('UK_Biobank_PMIDs_' + str(i) + '.txt',header = None,sep = '\t')
	ids = df.iloc[6:]
	new.loc[i,'PMID_no_UKB'] = df.iloc[3][0].split(' ')[-1]
	new.loc[i,'UKB_no_pubmed'] = int(df.iloc[1][0].split(' ')[-1]) - (int(df.iloc[4][0].split(' ')[-1]) - int(df.iloc[3][0].split(' ')[-1]))
	new.loc[i,'on_both'] = int(df.iloc[4][0].split(' ')[-1]) - int(df.iloc[3][0].split(' ')[-1]) - new.loc[i,'UKB_no_pubmed']

	new.loc[i,'Total_pubs'] = df.iloc[4][0].split(' ')[-1]
print(new)
new.to_csv('summary.txt',sep = '\t')

new = pd.read_csv('summary.txt',sep = '\t')
new.columns = ['Year','on_both','PMID_no_UKB','UKB_no_pubmed','Total_pubs']
print(new)
labels = new['Year']
both = new['on_both']
PMID_no_UKB = new['PMID_no_UKB']
UKB_no_pubmed = new['UKB_no_pubmed']
width = 0.35

fig, ax = plt.subplots()
ind = np.arange(2012,2021)
ax.bar(labels, both, width, label='PMID and UKB Website',color = 'purple',alpha = 0.75,align = 'center',edgecolor = 'black',hatch = 'x')
ax.bar(labels, PMID_no_UKB, width, bottom=both, label='PMID and not on UKB Website',color = 'red',alpha = 0.75,align = 'center',edgecolor = 'black',hatch = '/')
ax.bar(labels, UKB_no_pubmed,width, bottom = both+PMID_no_UKB,label = 'No PMID and on UKB Website',color = 'blue',alpha = 0.75,align = 'center',edgecolor = 'black',hatch = 'o')
ax.set_ylabel('Number of studies')
ax.set_title('UK Biobank Studies from 2012 to 2020')
plt.xticks(ind,['2012','2013','2014','2015','2016','2017','2018','2019','2020'])
ax.legend(loc = 'upper left')
plt.tight_layout()
plt.savefig('PMID_UKB.pdf')



new_array = []
ukb = pd.read_csv('UKBiobank_entries.csv').set_index('PUBMEDID')
pmids = set()
for k in ukb.index.tolist():
	pmids.add(int(k))

for i in range(2012,2021):
	df = pd.read_csv('UK_Biobank_PMIDs_' + str(i) + '.txt',header = None,sep = '\t')
	print(df)
	ids = df.iloc[5:][0]
	ids = ids.tolist()
	newids = set()
	for j in ids:
		newids.add(int(j))

	europe_only = 0
	multi = 0
	x = pmids.intersection(newids)
	for study in x:
		temp = ukb.loc[study]
		if temp.shape == (14,):
			ancestries = [temp['ANCESTRY']]
			if ancestries == 'European_ancestry':
				europe_only +=1
 		else:
			ancestries = temp['ANCESTRY'].unique()
			if len(ancestries) == 1:
				if ancestries == 'European_ancestry':
					europe_only += 1
				else:
					pass
			elif len(ancestries) > 1:
				multi +=1
	new_array.append([i,europe_only,multi])

new_array = np.array(new_array)
print(new_array)
width = 0.35
labels = new_array[:,0]
fig, ax = plt.subplots()
ind = np.arange(2012,2021)
ax.bar(labels, new_array[:,1], width, label='European only',color = 'blue',alpha = 0.75,align = 'center',edgecolor = 'black',hatch = 'x')
ax.bar(labels, new_array[:,2], width, label='Multi-ethnic',color = 'red',alpha = 0.75,align = 'center',bottom = new_array[:,1],edgecolor = 'black',hatch = '.')
ax.set_ylabel('Number of studies')
ax.set_title('UK Biobank Studies logged in the GWAS Catalog 2012 to 2020')
plt.xticks(ind,['2012','2013','2014','2015','2016','2017','2018','2019','2020'])
ax.legend(loc = 'upper left')
plt.tight_layout()
plt.savefig('PMID_UKB_GWACAT.pdf')


