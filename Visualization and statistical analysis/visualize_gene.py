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
# name = sys.argv[1]

# with open('inputs/vectors_fixed_ite50.json','r') as f:
# 	all_vectors = json.load(f)

# MGs = { key:value for key,value in all_vectors.items() if "MGI" in key and "http" not in key}
# HDs = { key:value for key,value in all_vectors.items() if (("OMIM" in key or "ORPHA" in key or "DECIPHER" in key) and "http" not in key)}
# #YGs = { key:value for key,value in all_vectors.items() if "SP" in key and "http" not in key}
# #FBGs = { key:value for key,value in all_vectors.items() if "FBgn" in key and "http" not in key}
# #ZFGs = { key:value for key,value in all_vectors.items() if "ZDB-GENE" in key and "http" not in key}


# print(len(MGs.keys()))
# print(len(HDs.keys()))
# #print(len(YGs.keys()))
# #print(len(FBGs.keys()))
# #print(len(ZFGs.keys()))



# ##viguilizing words in model

# X = list(MGs.values())+list(HDs.values())
# print("number of entities to rpresent",len(X))
# color = ['coral']*len(MGs.keys())+['beige']*len(HDs.keys())

# tsne = TSNE(n_components=2)
# X_tsne = tsne.fit_transform(X)
# print("done TSNE")
# df = pd.concat([pd.DataFrame(X_tsne),pd.Series(color)],axis=1)
# df.columns = ['x', 'y', 'color']
# #plt.rcParams.update({'font.size': 0.8})
# #colors=['coral','beige','lime','khaki']
# plt.rcParams.update({'lines.markersize': 0.005})
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1)

# ax.scatter(df['x'], df['y'],c=df['color'])

# #for i, txt in enumerate(df['word']):
# #    ax.annotate(txt, (df['x'].iloc[i], df['y'].iloc[i]))
# plt.show()

#-----------------------------------------------------
#other
organism = "human"


file1 = "Data/"+organism+"_gene_pheno_tab.txt"
file2 = "Data/"+organism+"normalized_information_content.tsv"
file3 = "Data/"+organism+"ranking_dic_genes_DL2Vec_withoutORPHA_fix_negatives_5_.json"
file4 = "Data/"+organism+"disease_gene_associations.tsv"

# visualize avarage IC and avarage rank

gene2annotations = {}
with open(file1) as f:
	content = f.readlines()
for mapping in content:
	if (len(mapping.split())>1):
		ls = mapping.split()
		gene2annotations[ls[0]] = ls[1:]

annotaion2IC = {}
with open(file2) as f:
	content = f.readlines()
for mapping in content:
	if (len(mapping.split())==2):
		ls = mapping.split()
		annotaion2IC[ls[0]] = float(ls[1])

'''
gene2rank = {}
with open(file3,'r') as f:
	gene2rank= json.load(f)
'''

id2name = {}
with open("class_id2name.json",'r') as f:
	id2name= json.load(f)

included_annotaions={}
for gene in gene2annotations:
	ls = [gene]
	ls2 = [gene]
	for annot in gene2annotations[gene]:
		key = annot.replace("http://purl.obolibrary.org/obo/",'').replace('_',":")
		if key in id2name:
			name = id2name[key]
			if not ((name.startswith('viable')) or (name.startswith('Viable')) or (name.startswith('normal')) or (name.startswith('Normal'))):
				ls.append(annot)
				ls2.append(name)
			#else:
				#print (name)
		#else:
			#print(key)

	if(len(ls)>1):
		print("\t".join(ls2))
		included_annotaions[gene]=ls[1:]

'''
#filter input for PLR training

print(sys.argv[1])
with open(sys.argv[1]) as f:
	content = f.readlines()
for mapping in content:
	if (len(mapping.split())>1):
		ls = mapping.split()
		if (ls[1] in included_annotaions):
			print(mapping.strip())
'''

'''
#filter gene sets
genes = []
with open("input/yeast/yeastGs_keys.json",'r') as f:
	genes= json.load(f)


filtered_genes=[]
for gene in genes:
	if gene in included_annotaions:
		filtered_genes.append(gene)

with open("yeast2Gs_keys.json","w") as f:
	json.dump(filtered_genes,f)
'''

#convert the annotations to list





#plot 1

'''
entity2annotation = {}
with open("annotations_without_OrphaandDecipher.txt") as f:
	content = f.readlines()
for mapping in content:
	if (len(mapping.split())==2):
		ls = mapping.split()
		#print(ls)
		entity = ls[0]
		annotation = ls[1]
		if entity not in entity2annotation:
			entity2annotation[entity]=set([])
		entity2annotation[entity].add(annotation)
		if("OMIM" in entity):
			print(entity+'\t'+annotation)
		elif(entity in included_annotaions):
			if(annotation.replace("<","").replace(">","") in included_annotaions[entity]):
				print(entity+'\t'+annotation)

'''

'''
Genes_avarage_IC = []
Genes_avarage_Rank = []

for gene in gene2rank:
	annotations_IC = []
	for annotation in gene2annotations[gene]:
		annotations_IC.append(annotaion2IC[annotation])

	Genes_avarage_IC.append(sum(annotations_IC)/len(annotations_IC))
	Genes_avarage_Rank.append(sum(gene2rank[gene])/len(gene2rank[gene]))

plt.scatter(Genes_avarage_IC, Genes_avarage_Rank, c="g", alpha=0.5)
plt.xlabel("Genes Average IC")
plt.ylabel("Genes Average Rank")
plt.gca().invert_yaxis()
plt.show()
'''

#plot 2
'''
class_IC = []
class_average_Rank = []

annotaion2ranks = {}
for gene in gene2rank:
	for annotation in gene2annotations[gene]:
		if annotation not in annotaion2ranks:
			annotaion2ranks[annotation] = []
		annotaion2ranks[annotation] = annotaion2ranks[annotation]+gene2rank[gene]


for annotation in annotaion2ranks:
	class_IC.append(annotaion2IC[annotation])
	class_average_Rank.append(sum(annotaion2ranks[annotation])/len(annotaion2ranks[annotation]))

plt.scatter(class_IC, class_average_Rank, c="g", alpha=0.5)
plt.xlabel("Class IC")
plt.ylabel("Class Average Rank")
plt.gca().invert_yaxis()
plt.show()
'''


#plot 3
'''
class_IC = []
class_Rank = []

annotaion2ranks = {}
for gene in gene2rank:
	for annotation in gene2annotations[gene]:
		for rank in gene2rank[gene]:
			class_IC.append(annotaion2IC[annotation])
			class_Rank.append(rank)

plt.scatter(class_IC, class_Rank, c="g", alpha=0.5)
plt.xlabel("Class IC")
plt.ylabel("Class Rank")
plt.gca().invert_yaxis()
plt.show()
'''

#plot 4
'''
associated_genes = set([])
with open("disease_gene_associations.tsv") as f:
	content = f.readlines()
for mapping in content:
	if (len(mapping.split())==2):
		ls = mapping.split()
		associated_genes.add(ls[1])


Genes_avarage_IC = []
Genes_avarage_Rank = []
colors = []

for gene in gene2rank:
	annotations_IC = []
	for annotation in gene2annotations[gene]:
		annotations_IC.append(annotaion2IC[annotation])

	Genes_avarage_IC.append(sum(annotations_IC)/len(annotations_IC))
	Genes_avarage_Rank.append(sum(gene2rank[gene])/len(gene2rank[gene]))
	if gene in associated_genes:
		colors.append('g')
	else:
		colors.append('r')

plt.scatter(Genes_avarage_IC, Genes_avarage_Rank, c=colors, alpha=0.3)
plt.xlabel("Genes Average IC")
plt.ylabel("Genes Average Rank")
plt.gca().invert_yaxis()
plt.show()
'''

#plot 5
'''
associated_genes = set([])
with open("disease_gene_associations.tsv") as f:
	content = f.readlines()
for mapping in content:
	ls = mapping.split()
	associated_genes.add(ls[1])

print(associated_genes)

Genes_sum_IC = []
Genes_avarage_Rank = []
colors = []

for gene in gene2rank:
	if gene in gene2annotations:
		annotations_IC = []
		for annotation in gene2annotations[gene]:
			annotations_IC.append(annotaion2IC[annotation])

		Genes_sum_IC.append(sum(annotations_IC))
		Genes_avarage_Rank.append(sum(gene2rank[gene])/len(gene2rank[gene]))
		if gene in associated_genes:
			colors.append('g')
			#print(gene,sum(annotations_IC))
		else:
			colors.append('r')

plt.scatter(Genes_sum_IC, Genes_avarage_Rank, c=colors, alpha=0.3)
plt.xlabel("Genes Sum IC")
plt.ylabel("Genes Average Rank")
plt.gca().invert_yaxis()
plt.show()
'''

# plot 6 histogram

n_bins = 50
associated_genes = set([])
with open(file4) as f:
	content = f.readlines()
for mapping in content:
	ls = mapping.split()
	associated_genes.add(ls[1])
Genes_sum_IC_associated = []
Genes_sum_IC_non_associated = []


for gene in gene2annotations:
	annotations_IC = []
	for annotation in gene2annotations[gene]:
		annotations_IC.append(annotaion2IC[annotation])

	if gene in associated_genes:
		Genes_sum_IC_associated.append(sum(annotations_IC))
	else:
		Genes_sum_IC_non_associated.append(sum(annotations_IC))



plt.hist(Genes_sum_IC_associated, bins=n_bins , density=1, alpha=0.5, label='associated genes')
plt.hist(Genes_sum_IC_non_associated, bins=n_bins, density=1 ,alpha=0.5, label='non-associated genes')
plt.legend(loc='upper right')
plt.show()

import numpy as np
from scipy import stats

Genes_sum_IC_associated.sort()
Genes_sum_IC_non_associated.sort()

t2, p2 = stats.ttest_ind(a=Genes_sum_IC_associated,b=Genes_sum_IC_non_associated, alternative="greater")
print("t = " + str(t2))
print("p = " + str(p2))

