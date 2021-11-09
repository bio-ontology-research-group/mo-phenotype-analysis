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

#-----------------------------------------------------
#other
organism = sys.argv[1]
file1 = "Data/"+organism+"_resnik_input.tsv"
file2 = "Data/"+organism+"_2021_IC.tsv"
#file3 = "input/"+organism+"/ranking_dic_genes_DL2Vec_withoutORPHA_fix_negatives_5_.json"
file4 = "Data/"+organism+"/disease_gene_associations.tsv"
file5 = "Data/"+organism+"Gs_keys.json"


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
with open( "Data/class_id2name.json",'r') as f:
	id2name= json.load(f)

included_annotaions={}
for gene in gene2annotations:
	ls = [gene]
	for annot in gene2annotations[gene]:
		key = annot.replace("http://purl.obolibrary.org/obo/",'').replace('_',":")
		if key in id2name:
			name = id2name[key]
			if not ((name.startswith('viable')) or (name.startswith('Viable')) or (name.startswith('normal')) or (name.startswith('Normal'))):
				ls.append(annot)
			#else:
				#print (name)
		#else:
			#print(key)

	if(len(ls)>1):
		#print("\t".join(ls))
		included_annotaions[gene]=ls[1:]

#matrex generation

diseases = set([])
Genes_frequencies = {}

with open(file4) as f:
	content = f.readlines()
for mapping in content:
	ls = mapping.split()
	diseases.add(ls[0])
	if ls[1] not in Genes_frequencies:
		Genes_frequencies[ls[1]]=1
	else:
		Genes_frequencies[ls[1]]+=1


genelist = []
with open(file5,'r') as f:
	genelist= json.load(f)



Genes_sum_IC = []

for gene in genelist:
	annotations_IC = []
	for annotation in gene2annotations[gene]:
		annotations_IC.append(annotaion2IC[annotation])
	Genes_sum_IC.append(sum(annotations_IC))


diseases_dictionary ={}
for disease in diseases:
	diseases_dictionary[disease]=Genes_sum_IC

with open(organism+"_naive_IC.json","w") as f:
	json.dump(diseases_dictionary,f)


