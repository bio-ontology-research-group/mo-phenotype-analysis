# Sarah M. Alghamdi
#------------------------------------
# This code is used to make the statistical analysis of MIS_HP (the most informative super class from HP)
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
organism = "fly"


file1 = "Data/"+organism+"_gene_pheno_tab.txt"
file2 = "Data/"+organism+"_unnormalized_IC.txt"
file3 = "Data/"+"Pheno_e_all_parents.json"
file4 = "Data/"+organism+"/disease_gene_associations.tsv"


avg_test = True
sum_test =False 
abs_diff = True
MIS_HP = False


# visualize avarage IC and avarage rank

gene2annotations = {}
with open(file1) as f:
	content = f.readlines()
for mapping in content:
	if (len(mapping.split())>1):
		ls = mapping.split()
		gene2annotations[ls[0]] = list(set(ls[1:]))

annotaion2IC = {}
with open(file2) as f:
	content = f.readlines()
for mapping in content:
	if (len(mapping.split())==2):
		ls = mapping.split()
		annotaion2IC[ls[0].replace("http://aber-owl.net/phenotype.owl#",'').replace('http://purl.obolibrary.org/obo/','').replace("_",":")] = float(ls[1])

id2name = {}
with open("Data/class_id2name.json",'r') as f:
	id2name = json.load(f)

superclasses = {}
with open(file3,'r') as f:
	superclasses = json.load(f)	

included_annotaions={}
for gene in gene2annotations:
	ls = [gene]
	ls2 = [gene]
	for annot in gene2annotations[gene]:
		key = annot.replace("http://aber-owl.net/phenotype.owl#",'').replace('http://purl.obolibrary.org/obo/','').replace("_",":")
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


entity2annotation = {}
with open("Data/annotations_without_OrphaandDecipher.txt") as f:
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
		#if("OMIM" in entity):
		#	print(entity+'\t'+annotation)
		#elif(entity in included_annotaions):
		#	if(annotation.replace("<","").replace(">","") in included_annotaions[entity]):
		#		print(entity+'\t'+annotation)



# plot 6 histogram

n_bins = 50
associated_genes = set([])
with open(file4) as f:
	content = f.readlines()
for mapping in content:
	ls = mapping.split()
	associated_genes.add(ls[1])
Genes_avg_IC_associated = []
Genes_avg_IC_non_associated = []


for gene in gene2annotations:
	annotations_IC = []
	for annotation in gene2annotations[gene]:
		annot = annotation.replace("http://aber-owl.net/phenotype.owl#",'').replace('http://purl.obolibrary.org/obo/','').replace("_",":")
		MIHP = ''
		IC = -1
		if annot in superclasses:
			for cl in superclasses[annot]:
				if("HP" in cl and IC <= annotaion2IC[cl] ):
					MIHP = cl
					IC = annotaion2IC[cl]
			if (MIHP is not ''):
				print(annot,MIHP,abs(annotaion2IC[annot]-IC))
				if(abs_diff):
					annotations_IC.append(abs(annotaion2IC[annot]-IC))
				if(MIS_HP):
					annotations_IC.append(IC)
		#else:
			#print("annotation not found to have superclasses", annot)
	if (len(annotations_IC)>0):
		if gene in associated_genes:
			if(avg_test):
				Genes_avg_IC_associated.append(sum(annotations_IC)*1.0/len(annotations_IC))
			if(sum_test):
				Genes_avg_IC_associated.append(sum(annotations_IC))
		else:
			if(avg_test):
				Genes_avg_IC_non_associated.append(sum(annotations_IC)*1.0/len(annotations_IC))
			if(sum_test):
				Genes_avg_IC_non_associated.append(sum(annotations_IC))
	#else:
		#print("gene not have annotation with superclasses from HP", gene)
		

print("number of annotated associated genes",len(Genes_avg_IC_associated))
print("number of annotated non associated genes",len(Genes_avg_IC_non_associated))


plt.hist(Genes_avg_IC_associated, bins=n_bins , density=1, alpha=0.5, label='associated genes')
plt.hist(Genes_avg_IC_non_associated, bins=n_bins, density=1 ,alpha=0.5, label='non-associated genes')
plt.legend(loc='upper right')
plt.show()

import numpy as np
from scipy import stats

Genes_avg_IC_associated.sort()
Genes_avg_IC_non_associated.sort()

#t2, p2 = stats.ttest_ind(a=Genes_avg_IC_associated,b=Genes_avg_IC_non_associated)#, alternative="greater")
#print("t = " + str(t2))
#print("p = " + str(p2))

print("max of IC for associated genes" , max(Genes_avg_IC_associated))
print("max of IC for non associated genes" , max(Genes_avg_IC_non_associated))

print("min of IC for associated genes" , min(Genes_avg_IC_associated))
print("min of IC for non associated genes" , min(Genes_avg_IC_non_associated))

print("mean of IC for associated genes" , sum(Genes_avg_IC_associated)/len(Genes_avg_IC_associated))
print("mean of IC for non associated genes" , sum(Genes_avg_IC_non_associated)/len(Genes_avg_IC_non_associated))


def meanof50(ls):
	ls.sort()
	start = int(len(ls)/4)
	end = int(3*len(ls)/4)
	print("median = "+str(ls[int(len(ls)/5)]))
	mean50 = sum(ls[start:end])/len(ls[start:end])
	return mean50



print("mean of IC from 25 percintile tell the 75 percintile for associated genes" , meanof50(Genes_avg_IC_associated))
print("mean of IC from 25 percintile tell the 75 percintile for non associated genes" , meanof50(Genes_avg_IC_non_associated))

