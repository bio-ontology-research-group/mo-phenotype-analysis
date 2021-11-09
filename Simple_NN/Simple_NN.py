# Sarah M. Alghamdi
#------------------------------------
# this is a simple implemantation of MLP takes vectors of genes and diseases and positives dictionary file
# it does the training on 10-fold-cross validation and calculates the rank roc auc 
# 
#
# inputs :
# genes_vectors_filename : json dictionary {"gene_id":vector of real numbers as a list}
# diseases_vectors_filename : json dictionary {"disease_id":vector of real numbers as a list}
# positives_filename : json dictionary {"disease_id": list of gene ids}
#------------------------------------------------------------------------

from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
import json
import sys
import numpy as np 
import random
import math

genes_vectors_filename = sys.argv[1]
diseases_vectors_filename = sys.argv[2]
positives_filename = sys.argv[3]
method = sys.argv[4]
hard = False
ratio = 5

Vector_size = 200


# This method is used to generate negative samples 
# input:  
# # genes_keys: 
# # diseases_keys:
# # positives:
# # hard: 
# output:
# # negatives:
# # new_positives:
# # pos_count:
# # neg_count:
# 
# When data are generated the negative genes are selected in 2 ways: hard choise will select the negative genes from the disease associated genes only,
# not hard when the selection of the genes are from associated and non associated genes. 
def generate_negatives(genes_keys, diseases_keys, positives, hard):
	negatives = {}
	new_positives = {}
	pos_count = 0
	neg_count = 0
	disease_associated_genes = set([])
	for disease in positives:
		if (disease in diseases_keys):
			for gene in positives[disease]:
				if(gene in genes_keys):
					if(disease not in new_positives):
						new_positives[disease]=set([])
					pos_count+=1
					disease_associated_genes.add(gene)
					new_positives[disease].add(gene)
	non_disease_associated_genes = set([])
	for gene in genes_keys:
		if gene not in disease_associated_genes:
			non_disease_associated_genes.add(gene)

	#genes can be associated or non associated genes
	if not hard: 
		for disease in diseases_keys:
			if disease in positives:
				negatives[disease] = set([])
				for gene in genes_keys:
					neg_count+=1
					negatives[disease].add(gene)

	#genes are only the associated genes
	if hard:
		for disease in diseases_keys:
			if disease in positives:
				negatives[disease] = set([])
				for gene in genes_keys:
					if (gene not in positives[disease]) and gene not in non_disease_associated_genes:
						neg_count+=1
						negatives[disease].add(gene)
						break
	return negatives,new_positives, pos_count, neg_count

def get_input_analysis(genes_vectors_filename, diseases_vectors_filename, positives_filename):
	genes_vectors = {}
	with open(genes_vectors_filename,'r') as f:
		genes_vectors = json.load(f)

	diseases_vectors = {}
	with open(diseases_vectors_filename,'r') as f:
		diseases_vectors = json.load(f)

	positives = {}
	with open(positives_filename,'r') as f:
		positives = json.load(f)

	diseases_keys = list(diseases_vectors.keys())
	genes_keys = list(genes_vectors.keys())

	new_positives={}
	for disease in positives:
		if (disease in diseases_keys):
			for gene in positives[disease]:
				if(gene in genes_keys):
					if(disease not in new_positives):
						new_positives[disease]=set([])
					new_positives[disease].add(gene)

	new_disease_keys = [x for x in diseases_keys if x in new_positives]

	print(len(new_disease_keys), len(genes_keys) , len(new_positives.keys()))

	return new_disease_keys,genes_keys,new_positives


def get_input(genes_vectors_filename, diseases_vectors_filename ,positives_filename, ratio):
	genes_vectors = {}
	with open(genes_vectors_filename,'r') as f:
		genes_vectors = json.load(f)

	diseases_vectors = {}
	with open(diseases_vectors_filename,'r') as f:
		diseases_vectors = json.load(f)

	positives = {}
	with open(positives_filename,'r') as f:
		positives = json.load(f)

	diseases_keys = list(diseases_vectors.keys())
	genes_keys = list(genes_vectors.keys())

	negatives, new_positives, pos_count, neg_count = generate_negatives(genes_keys, diseases_keys, positives, hard)


	# Defining Feature Matrex
	X= np.empty(((ratio+1)*pos_count,Vector_size*2))
	y= np.empty((ratio+1)*pos_count)

	negative_diseases = list(negatives.keys())
	sample_number=0
	for disease in new_positives:
		for gene in new_positives[disease]:
			x = np.concatenate((diseases_vectors[disease],genes_vectors[gene]),axis=0)
			X[sample_number]=x
			y[sample_number]=1
			sample_number+=1


			for i in range(ratio):
				n = random.randint(0,len(negative_diseases))
				n_disease = negative_diseases[n-1]
				n = random.randint(0,len(negatives[n_disease]))
				n_gene = list(negatives[n_disease])[n-1]
				x = np.concatenate((diseases_vectors[n_disease],genes_vectors[n_gene]),axis=0)
				X[sample_number]=x
				y[sample_number]=0
				sample_number+=1
	return X,y

def get_training_folds(genes_vectors_filename, diseases_vectors_filename ,positives,diseases_keys,genes_keys, ratio, fold):
	genes_vectors = {}
	with open(genes_vectors_filename,'r') as f:
		genes_vectors = json.load(f)

	diseases_vectors = {}
	with open(diseases_vectors_filename,'r') as f:
		diseases_vectors = json.load(f)

	start = int(len(diseases_keys)*fold/10)
	end = int(len(diseases_keys)*(fold+1)/10) - 1


	testing_disease_keys = diseases_keys[start:end]
	training_disease_keys = [x for x in diseases_keys if x not in testing_disease_keys]

	print(start,end,len(testing_disease_keys),len(training_disease_keys))

	negatives, new_positives, pos_count, neg_count = generate_negatives(genes_keys, training_disease_keys, positives, hard)


	# Defining Feature Matrex
	X= np.empty(((ratio+1)*pos_count,Vector_size*2))
	y= np.empty((ratio+1)*pos_count)

	negative_diseases = list(negatives.keys())
	sample_number=0

	for disease in new_positives:
		for gene in new_positives[disease]:
			x = np.concatenate((diseases_vectors[disease],genes_vectors[gene]),axis=0)
			X[sample_number]=x
			y[sample_number]=1
			sample_number+=1


			for i in range(ratio):
				n = random.randint(1,len(negative_diseases))
				n_disease = negative_diseases[n-1]
				n = random.randint(1,len(negatives[n_disease]))
				n_gene = list(negatives[n_disease])[n-1]
				x = np.concatenate((diseases_vectors[n_disease],genes_vectors[n_gene]),axis=0)
				X[sample_number]=x
				y[sample_number]=0
				sample_number+=1

	index = 0
	X_test= np.empty((len(testing_disease_keys)*len(genes_keys),Vector_size*2))
	y_test= np.empty(len(testing_disease_keys)*len(genes_keys))
	test_guide = {}
	for disease in testing_disease_keys:
		test_guide[disease] = {}
		for gene in genes_keys:
			test_guide[disease][gene] = index
			x = np.concatenate((diseases_vectors[disease],genes_vectors[gene]),axis=0)
			X_test[index]=x
			if(disease in new_positives):
				if(gene in new_positives[disease]):
					y_test[index]=1
				else:
					y_test[index]=0
			else:
				y_test[index]=0
			index+=1



	return X,y , X_test, y_test, test_guide



def computeROC_AUC(OGs_HDs_sim,HDs_keys,OGs_keys,h_positives):
    OGs_HDs_sim_t = OGs_HDs_sim
    ranks = np.argsort(OGs_HDs_sim_t, axis=1)
    TPR = [0]
    FPR = [0]
    prev= [0,0]
    P=0
    N=0
    positive_genes = set([])
    includedDiseases =  np.zeros(len(HDs_keys))
    print(OGs_HDs_sim_t.shape,len(HDs_keys),len(OGs_keys))

    positives_matrix = np.zeros([len(HDs_keys),len(OGs_keys)])
    for og in range(0,len(OGs_keys)):
        for hd in range(0,len(HDs_keys)):
            if(HDs_keys[hd] in  h_positives):
                if(OGs_keys[og] in h_positives[HDs_keys[hd]]):
                    positives_matrix[hd][og] = 1
                    includedDiseases[hd]=1
                    positive_genes.add(OGs_keys[og])
    ranking_dic_disease={}
    ranking_dic_gene={}
    positives_ranks = {}

    for hd in range(0,len(HDs_keys)):
        if(includedDiseases[hd]==1):
            p = np.sum(positives_matrix[hd])
            P+= p
            N+= len(OGs_keys)-p
    print("p",P)    
    print("included_Disease", np.sum(includedDiseases))
    for r in range(0,len(OGs_keys)):
        TP = prev[0]
        FP = prev[1]
        for hd in range(0,len(HDs_keys)):
            if(includedDiseases[hd]==1):
                g=len(OGs_keys)-r-1
                if (positives_matrix[hd,ranks[hd][g]] == 1):
                    TP+=1
                    if(HDs_keys[hd] not in ranking_dic_disease):
                        ranking_dic_disease[HDs_keys[hd]] = r+1
                        if (OGs_keys[g] not in ranking_dic_gene):
                            ranking_dic_gene[OGs_keys[g]]=[]
                        ranking_dic_gene[OGs_keys[g]].append(r+1)
                else:
                    FP+=1
        prev = [TP,FP]
        TPR.append(TP/P)
        FPR.append(FP/N)
    return TPR,FPR,np.trapz(TPR,FPR),P,N


'''
X, y = get_input(genes_vectors_filename, diseases_vectors_filename, positives_filename, ratio)
X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, shuffle=True ,random_state=42)
clf = MLPClassifier(hidden_layer_sizes=200, activation= "logistic", solver = "adam", alpha=0.0001, learning_rate= 'constant',learning_rate_init=0.001, random_state=42, max_iter=300).fit(X_train, y_train)
clf.predict_proba(X_test[:1])
clf.score(X_test, y_test)
print("Accuricy score is :")
print(clf.score(X_test, y_test))

'''

disease = 0
genes = 0

HDs_keys,OGs_keys,positives = get_input_analysis(genes_vectors_filename, diseases_vectors_filename, positives_filename)
OGs_HDs_sim = np.empty((len(HDs_keys),len(OGs_keys)))

for fold in range(10):
	print("-------------statring fold--------------")
	print(fold)
	X_train, y_train, X_test, y_test, test_guid = get_training_folds(genes_vectors_filename, diseases_vectors_filename, positives,HDs_keys, OGs_keys, ratio, fold)

	clf = MLPClassifier(hidden_layer_sizes=(Vector_size,), activation= "logistic", solver = "adam", alpha=0.0001, learning_rate= 'constant',learning_rate_init=0.001, random_state=42, max_iter=300, early_stopping=False).fit(X_train, y_train)
	result = clf.predict_proba(X_test)
	#print(clf.classes_)
	#print(result)
	print("filling the results")
	for d in range(0,len(HDs_keys)):
		disease = HDs_keys[d]
		if disease in test_guid:
			for g in range(len(OGs_keys)):
				gene=OGs_keys[g]
				index = test_guid[disease][gene]
				OGs_HDs_sim[d][g] = result[index][1]





TPR,FPR,ROC,P,N = computeROC_AUC(OGs_HDs_sim,HDs_keys,OGs_keys,positives)


#-----------------------------------------------------------
def precisionAt(TPR,P,FPR,N,at):
    return (TPR[at]*P)/((TPR[at]*P)+(FPR[at]*N))

def recallAt(TPR,P,FPR,N,at):
    return TPR[at]

def hitAt(TPR,P,FPR,N,at):
    return TPR[at]*P
#---------------------------------------------------------



print("Human Evaluation")
print("ROC_AUC = ",ROC)
print("CI = ",2*(math. sqrt(ROC*(1-ROC)/(min([P,N])))))
print("Precision@1 @10 @50 @100 = ",precisionAt(TPR,P,FPR,N,1),precisionAt(TPR,P,FPR,N,10),precisionAt(TPR,P,FPR,N,50),precisionAt(TPR,P,FPR,N,100))
print("Recall@1 @10 @50 @100 = ",recallAt(TPR,P,FPR,N,1),recallAt(TPR,P,FPR,N,10),recallAt(TPR,P,FPR,N,50),recallAt(TPR,P,FPR,N,100))
print("Hits@1 @10 @50 @100 = ",hitAt(TPR,P,FPR,N,1),hitAt(TPR,P,FPR,N,10),hitAt(TPR,P,FPR,N,50),hitAt(TPR,P,FPR,N,100))
print("Total positives = ", P)