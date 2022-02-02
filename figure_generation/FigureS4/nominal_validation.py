import pandas as pd 
import numpy as np 

traits = ['MCV','PLC','Height', 'BMI', 'DBP', 'SBP', 'WBC', 'RBC', 'Hemoglobin', 'Hematocrit', 'MCH', 'MCHC', 'Lymphocyte', 'Monocyte', 'Neutrophil', 'Eosinophil', 'Basophil','Urate','Triglyceride','Cholesterol','LDL','HDL','HBA1C','EGFR','CRP']


thresholds = pd.DataFrame(np.zeros((25,4)), columns = ['thresh','tested','replicated','missed'],index = traits)
for trait in traits:
	significance = pd.read_csv('../variant_manhattan_plot/' + trait + '.upset.input.csv')
	otr_sig = significance[significance['european'] == 0]
	significance = significance[significance['european'] == 1]
	pfile = pd.read_csv('../variant_manhattan_plot/' + trait + '.merged.snps.txt',sep = '\t').set_index('SNP')
	pfile = pfile.loc[significance['SNP'].tolist()]
	pfile = pfile.drop(['european','#CHROM','POS'], axis = 1)

	pfile['nans'] = pfile.isnull().sum(axis = 1)
	pfile = pfile[pfile['nans'] != len(pfile.columns)-1]
	thresholds.loc[trait,'tested'] = pfile.shape[0]
	thresh = 0.05/pfile.shape[0]
	pfile = pfile.drop(['nans'], axis = 1)
	siggy = pfile < thresh
	siggy = siggy*1
	siggy.to_csv('nominal.' + trait + '.upset.input.csv',index = True)
	siggy['replicates'] = siggy.sum(axis = 1)
	replicated = siggy[siggy['replicates'] >= 1]
	thresholds.loc[trait,'replicated'] = replicated.shape[0]
	thresholds.loc[trait,'missed'] = otr_sig.shape[0]

thresholds['thresh'] = 0.05/thresholds['tested']
thresholds['proportion'] = thresholds['replicated']/thresholds['tested']
thresholds[['thresh','tested','replicated','proportion','missed']].to_csv('nominal.metadata.csv')
thresholds[['thresh','tested','replicated','proportion','missed']].to_csv('nominal.metadata.latex.csv',sep = '&')


thresholds = pd.DataFrame(np.zeros((25,4)), columns = ['thresh','tested','replicated','missed'],index = traits)
for trait in traits:
	significance = pd.read_csv('../gene_manhattan_plot/' + trait + '.upset.input.csv')
	otr_sig = significance[significance['european'] == 0]
	significance = significance[significance['european'] == 1]
	pfile = pd.read_csv('../gene_manhattan_plot/' + trait + '.merged.genes.txt',sep = '\t').set_index('gene')
	pfile = pfile.loc[significance['gene'].tolist()]
	pfile = pfile.drop(['european','chr'], axis = 1)

	pfile['nans'] = pfile.isnull().sum(axis = 1)
	pfile = pfile[pfile['nans'] != len(pfile.columns)-1]
	thresholds.loc[trait,'tested'] = pfile.shape[0]
	thresh = 0.05/pfile.shape[0]
	pfile = pfile.drop(['nans'], axis = 1)
	siggy = pfile < thresh
	siggy = siggy*1
	siggy.to_csv('nominal.' + trait + '.gene.upset.input.csv',index = True)
	siggy['replicates'] = siggy.sum(axis = 1)
	replicated = siggy[siggy['replicates'] >= 1]
	thresholds.loc[trait,'replicated'] = replicated.shape[0]
	thresholds.loc[trait,'missed'] = otr_sig.shape[0]

thresholds['thresh'] = 0.05/thresholds['tested']
thresholds['proportion'] = thresholds['replicated']/thresholds['tested']
thresholds.to_csv('nominal.metadata.gene.csv')
thresholds.to_csv('nominal.metadata.gene.latex.csv',sep = '&')

