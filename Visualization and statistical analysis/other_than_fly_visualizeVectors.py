# Sarah M. Alghamdi
#------------------------------------
# args
# sim_file disease_file genes_file outputs_prefix
#------------------------------------
import json
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import pandas as pd
from sklearn.manifold import TSNE

#------------------------------------
# separate the vectors
file_name = "vectors.text"
phenotpes = 'Data/All_clean_filtered.txt'
phenotpes_annotaion = "FYPO" #MP"
gene_annotation = "SP" #"MGI" #
init = 'y_' #f_' # m_, f_ , 
mapping_file = '../yeast_human_mapping.json' #"../human_1_n_fish_mapping.json" #"../human_1_n_mouse_mapping.json" # 



with open(file_name) as f:
	content = f.readlines()
vectors = {}

F_human_p,F_OMIM,F_human,p_genes = {},set([]),{},set([])

for line in content:
	l = line.strip().split()
	vectors[l[0]] = l[1].split(',')

with open(mapping_file,'r') as f:
	F_human = json.load(f)
with open("../human_positives.tsv") as f:
	content = f.readlines()	
for p in content:
	hg = p.split()[1]
	if hg in F_human.keys():
		for g in F_human[hg]:
			for o in p.split()[0].split('|'):
				om=o.strip()
				if(om in F_human_p):
					F_human_p[om].append(g)
					F_OMIM.add(om)
					p_genes.add(g)
				else:
					F_human_p[om]=[g]
					F_OMIM.add(om)
					p_genes.add(g)

p_pheno = set([])
n_pheno = set([])



with open(phenotpes) as f:
	content = f.readlines()	
for p in content:
	ls = p.strip().split()
	if ls[0] in p_genes:
		p_pheno.add(ls[1])
	else:
		if gene_annotation in ls[0]:
			n_pheno.add(ls[1])

n_genes = set([])
for key in vectors:
	if gene_annotation in key and key not in p_genes:
		n_genes.add(key)



pos = list(p_pheno)
frequencies_pos = [0]*len(pos)
neg = list(n_pheno)
frequencies_neg = [0]*len(neg)

gene2phenotype = {}
with open(phenotpes) as f:
	content = f.readlines()	
for p in content:
	ls = p.strip().split()
	if(gene_annotation in ls[0]):
		if(ls[0] not in gene2phenotype):
			gene2phenotype[ls[0]] = set([])
		gene2phenotype[ls[0]].add(ls[1])

for p in range(len(pos)):
	for gene in p_genes:
		if gene in gene2phenotype:
			if pos[p] in gene2phenotype[gene]:
				frequencies_pos[p]+=1
for n in range(len(neg)):
	for gene in n_genes:
		if gene in gene2phenotype:
			if neg[n] in gene2phenotype[gene]:
				frequencies_neg[n]+=1

with open(init+'pos', 'w') as fb:
        json.dump(pos,fb)
with open(init+'neg', 'w') as fb:
        json.dump(neg,fb)

with open(init+'frequencies_neg', 'w') as fb:
        json.dump(frequencies_neg,fb)
with open(init+'frequencies_pos', 'w') as fb:
        json.dump(frequencies_pos,fb)


print("statictics")
print("disease-associated genes", len(p_genes))
print("non-disease-associated genes", len(n_genes))
print("disease-associated genes phenotypes", len(p_pheno))
print("non-disease-associated genes phenotypes", len(n_pheno))


##viguilizing words in model

P_Gs = { key:value for key,value in vectors.items() if key in p_genes}
N_Gs = { key:value for key,value in vectors.items() if key in n_genes}

Ds = { key:value for key,value in vectors.items() if key in F_OMIM}

PN_P_U = { key:value for key,value in vectors.items() if (key in p_pheno or key in n_pheno)}

phenotypes= list(PN_P_U.keys())
phenotypes_vectors = []
phenotypes_p,phenotypes_n,phenotypes_i = [],[],[]
for pheno in phenotypes:
	phenotypes_vectors.append(PN_P_U[pheno])
	if(pheno in p_pheno and pheno in n_pheno):
		phenotypes_p.append('purple')
		phenotypes_n.append('maroon')
		phenotypes_i.append('darkcyan')
	elif(pheno in p_pheno and pheno not in n_pheno):
		phenotypes_p.append('purple')
		phenotypes_n.append('tan')
		phenotypes_i.append('tan')
	elif(pheno in n_pheno and pheno not in p_pheno):
		phenotypes_p.append('tan')
		phenotypes_n.append('maroon')
		phenotypes_i.append('tan')
	else:
		phenotypes_p.append('tan')
		phenotypes_n.append('tan')
		phenotypes_i.append('tan')


#figure 1 
X = list(N_Gs.values())+list(Ds.values())+list(P_Gs.values())+phenotypes_vectors
print("number of entities to rpresent",len(X))
color1 = ['maroon']*len(N_Gs.keys())+['darkcyan']*len(Ds.keys())+['purple']*len(P_Gs.keys())+['tan']*len(phenotypes)
color2 = ['tan']*len(Ds.keys())+['tan']*len(P_Gs.keys())+['tan']*len(N_Gs.keys())+phenotypes_p
color3 = ['tan']*len(Ds.keys())+['tan']*len(P_Gs.keys())+['tan']*len(N_Gs.keys())+phenotypes_n
color4 = ['tan']*len(Ds.keys())+['tan']*len(P_Gs.keys())+['tan']*len(N_Gs.keys())+phenotypes_i
color5 = ['tan']*len(N_Gs.keys())+['darkcyan']*len(Ds.keys())+['tan']*len(P_Gs.keys())+['tan']*len(phenotypes)
color6 = ['maroon']*len(N_Gs.keys())+['tan']*len(Ds.keys())+['tan']*len(P_Gs.keys())+['tan']*len(phenotypes)
color7 = ['tan']*len(N_Gs.keys())+['tan']*len(Ds.keys())+['purple']*len(P_Gs.keys())+['tan']*len(phenotypes)

tsne = TSNE(n_components=2)
X_tsne = tsne.fit_transform(X)
np.save('tsne_fly.npy',X_tsne)
print("done TSNE")
df = pd.concat([pd.DataFrame(X_tsne),pd.Series(color1)],axis=1)
df.columns = ['x', 'y', 'color']
#plt.rcParams.update({'font.size': 0.8})
#colors=['coral','beige','lime','khaki']
plt.rcParams.update({'lines.markersize': 1})
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.scatter(df['x'], df['y'],c=df['color'],alpha=0.5)
#for i, txt in enumerate(df['word']):
#    ax.annotate(txt, (df['x'].iloc[i], df['y'].iloc[i]))
plt.show()



#figure 2 
df2 = pd.concat([pd.DataFrame(X_tsne),pd.Series(color2)],axis=1)
df2.columns = ['x', 'y', 'color']
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.scatter(df2['x'], df2['y'],c=df2['color'],alpha=0.5)
plt.show()

#figure 3 
df3 = pd.concat([pd.DataFrame(X_tsne),pd.Series(color3)],axis=1)
df3.columns = ['x', 'y', 'color']
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.scatter(df3['x'], df3['y'],c=df3['color'],alpha=0.5)
plt.show()

#figure 4 
df4 = pd.concat([pd.DataFrame(X_tsne),pd.Series(color4)],axis=1)
df4.columns = ['x', 'y', 'color']
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.scatter(df4['x'], df4['y'],c=df4['color'],alpha=0.5)
plt.show()


#figure 5 disease 
df5 = pd.concat([pd.DataFrame(X_tsne),pd.Series(color5)],axis=1)
df5.columns = ['x', 'y', 'color']
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.scatter(df5['x'], df5['y'],c=df5['color'],alpha=0.5)
plt.show()

#figure 6 p genes 
df6 = pd.concat([pd.DataFrame(X_tsne),pd.Series(color6)],axis=1)
df6.columns = ['x', 'y', 'color']
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.scatter(df6['x'], df6['y'],c=df6['color'],alpha=0.5)
plt.show()

#figure 7 n genes 
df7 = pd.concat([pd.DataFrame(X_tsne),pd.Series(color7)],axis=1)
df7.columns = ['x', 'y', 'color']
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.scatter(df7['x'], df7['y'],c=df7['color'],alpha=0.5)
plt.show()
