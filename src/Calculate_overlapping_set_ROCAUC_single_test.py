# Sarah M. Alghamdi
#------------------------------------
# arg: 
# org_simle // sim_file // hunam_disease_file // organism_gene_file // organism-human-mapping // organism-mouse-mapping
#------------------------------------
import json
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
#------------------------------------
test_tag = "Resnik_similarity"
#test_tag = "Resnik_IC"
#test_tag = "OPA2Vec-Cosine"
#test_tag = "OWL2Vec-Cosine"
#test_tag = "DL2Vec-PLR"
#test_tag = "DL2Vec-cosine"

#---------------------------------------------------------
# ROC_AUC
# computeROC_AUC
# inputs:
#   OGs_HDs_sim : a numpy matrex for genes-disease similarity rows for organism genes, columns for human disease
#   HDs_keys    : a list representing the human disease header for OGs_HDs_sim
#   OGs_keys    : a list representing the gene disease header for OGs_HDs_sim
#   p_OMIMs     : list of posiitive OMIMs that are found eaither assosiated with human gene or assosiated with organism gene
#   p_OGs       : list of posiitive organism genes that are found assosated with disease
#   o_positives : positives pairs of (organism gene, disease)
#   h_positives : positive pairs of (human gene, disease)
#   o_or_h      : 'o' if you want the ROC_AUC computed considering the organism genes as positives , 'h' for human genes as positives
# outputs:
#   TPR         : a list of the true positive rate at the consedered ranks thresholds
#   FPR         : a list of the false positive rate at the consedered ranks thresholds
#   ROC_AUC     : the AUC calculated using trapz method

def computeROC_AUC(OGs_HDs_sim,HDs_keys,OGs_keys,o_genes,o_positives,h_positives,o_or_h,organism):
    OGs_HDs_sim_t = OGs_HDs_sim
    if(test_tag == "OPA2Vec-Cosine" or test_tag == "OWL2Vec-Cosine"  or test_tag == "Resnik_similarity" or test_tag == "DL2Vec-cosine"):
        OGs_HDs_sim_t = OGs_HDs_sim.transpose()
    ranks = np.argsort(OGs_HDs_sim_t, axis=1)
    TPR = [0]
    FPR = [0]
    prev= [0,0]
    P=0
    N=0
    includedDiseases =  np.zeros(len(HDs_keys))
    includedGenes = np.zeros(len(OGs_keys))
    positive_genes = set([])
    #check if the gene have human orthog
    for og in range(0,len(OGs_keys)):
        if (OGs_keys[og] in o_genes):
            includedGenes[og]=1

    print(OGs_HDs_sim_t.shape,len(HDs_keys),len(OGs_keys))

    positives_matrix = np.zeros([len(HDs_keys),len(OGs_keys)])
    for og in range(0,len(OGs_keys)):
        if(includedGenes[og]==1):
            for hd in range(0,len(HDs_keys)):
                if o_or_h == 'o':
                    if(HDs_keys[hd] in  o_positives):
                        if(OGs_keys[og] in o_positives[HDs_keys[hd]] and includedGenes[og]):
                            positives_matrix[hd][og] = 1
                            includedDiseases[hd]=1
                            positive_genes.add(OGs_keys[og])

                if o_or_h == 'h':
                    #print(HDs_keys[hd])
                    if(HDs_keys[hd] in  h_positives):
                        if(OGs_keys[og] in h_positives[HDs_keys[hd]] and includedGenes[og]):
                            positives_matrix[hd][og] = 1
                            includedDiseases[hd]=1
                            positive_genes.add(OGs_keys[og])
                            


    print("D = ",np.sum(includedDiseases),"G = ",np.sum(includedGenes))
    for hd in range(0,len(HDs_keys)):
        if(includedDiseases[hd]==1):
            p = np.sum(positives_matrix[hd])
            P+= p
            N+= sum(includedGenes)-p

    ls = [HDs_keys[x] for x in range(0,len(HDs_keys)) if includedDiseases[x]==1]
    with open(sys.argv[1]+".json","w") as fb:
        json.dump(ls,fb)

        
    print("p",P)
    print("positive_genes",len(positive_genes))
    print("n",N)
    # r is the rank considered
    ranking_dic_disease={}

    last_index = len(HDs_keys)*[len(OGs_keys)]
    for r in range(0,int(sum(includedGenes))):
        #print("rank",r)
        TP = prev[0]
        FP = prev[1]

        for hd in range(0,len(HDs_keys)):
            if(includedDiseases[hd]==1):
                #loop throgh diseases
                index = last_index[hd]
                while True:
                    index-=1
                    #print(index)
                    if(includedGenes[ranks[hd][index]] == 1):
                        break

                last_index[hd]=index

                if (positives_matrix[hd,ranks[hd][index]] == 1):
                    TP+=1
                    if(HDs_keys[hd] not in ranking_dic_disease):
                        ranking_dic_disease[HDs_keys[hd]] = r+1
                else:
                    FP+=1
                
        prev = [TP,FP]

        TPR.append(TP/P)
        FPR.append(FP/N)
        if(math.isnan(TP/N) or math.isnan(FP/N)):
            print(TP,FP,N,r,d,og,OGs_keys[og])
            return
    with open('overlapped_disease_ranking_'+test_tag+'_'+organism+'.json', 'w') as fb:
        json.dump(ranking_dic_disease,fb)


    print('genes ',np.sum(includedGenes),' associations ',P)
    return TPR,FPR,np.trapz(TPR,FPR),P,N

#-----------------------------------------------------------

def precisionAt(TPR,P,FPR,N,at):
    return (TPR[at]*P)/((TPR[at]*P)+(FPR[at]*N))

def recallAt(TPR,P,FPR,N,at):
    return TPR[at]

def hitAt(TPR,P,FPR,N,at):
    return TPR[at]*P



#---------------------------------------------------------
m,z,fly,y,h = False,False,False,False,False
if("mouse" in sys.argv[1]):
    m = True
if("fish" in sys.argv[1]):
    z = True
if("fly" in sys.argv[1]):
    fly = True
if("yeast" in sys.argv[1]):
    y = True
#---------------------------------------------------------
with open('Data/annotation_allspecies_human_genes.tsv', 'r') as f:
    content = f.readlines()
evaluation_set = []
genes = {}
for line in content:
    hg = line.split()[0]
    if hg not in genes:
        genes[hg] = []
    phenotype = line.split()[1]
    if "MP_" in phenotype:
        genes[hg].append("m")
    if "FBcv_" in phenotype or "FBbtAB_" in phenotype:
        genes[hg].append("f")
    if "#PHENO:" in phenotype:
        genes[hg].append("z")
    if "FYPO_" in phenotype:
        genes[hg].append("y")

for hg in genes:
    include = True
    if m and "m" not in genes[hg]:
        include = False
    if z and "z" not in genes[hg]:
        include = False
    if fly and "f" not in genes[hg]:
        include = False
    if y and "y" not in genes[hg]:
        include = False
    if (include):
        evaluation_set.append(hg)
print("the Evaluation data set informations")
print(list(set(evaluation_set)))
print("number og genes",len(list(set(evaluation_set))))

#--------------------------------------------------------
ROC_ls = []
Genes_ls = []
CI_ls = []
Associations_ls=[]
H_ls=[[],[],[]]
P_ls=[[],[],[]]
R_ls=[[],[],[]]

m_ROC_ls = []
m_Genes_ls = []
m_CI_ls = []
m_Associations_ls=[]
m_H_ls=[[],[],[]]
m_P_ls=[[],[],[]]
m_R_ls=[[],[],[]]
#--------------------------------------------------------
    # mouse human data - OPA2Vec similarities
humanGenes,OMIMs,MouseGenes = [],[],[]
m_human_mapping,h_positives,m_positives = {},{},{}
# getting human positive set
with open('Data/human_1_n_mouse_mapping.json','r') as f:
        m_human_mapping = json.load(f)

with open("Data/human_positives.tsv") as f:
    content = f.readlines()

count = 0
for p in content:
    hg = p.split()[1]
    if hg in evaluation_set:
        humanGenes.append(hg)
        for om in p.split()[0].split('|'):
            count+=1
            if(om in h_positives):
                h_positives[om].append(hg)
                OMIMs.append(om)
            else:
                h_positives[om]=[hg]
                OMIMs.append(om)
print(len(list(set(humanGenes))))
print(len(list(set(OMIMs))))
print(count)

'''
# getting mice positive set
count=0
OMIMs=[]
with open("new_mouse_positives") as f:
    content = f.readlines()
for p in content:
    mg = p.split()[1]
    MouseGenes.append(mg)
    for om in p.split()[0].split('|'):
        count+=1
        if(om in m_positives):
            m_positives[om].append(mg)
            OMIMs.append(om)
        else:
            m_positives[om]=[mg]
            OMIMs.append(om)
print(len(list(set(MouseGenes))))
print(len(list(set(OMIMs))))
print(count)
'''
HGs_keys = ""
HDs_keys = ""
HGs_HDs_sim = ""

if(test_tag == "DL2Vec-PLR"):#
    with open("Data/"+test_tag+"_human.json",'r') as f:
        sim_dic = json.load(f)
    items =  list(sim_dic.items())
    HDs_keys = [x[0] for x in items]
    HGs_HDs_sim = np.array([x[1] for x in items])
    with open("Data/"+test_tag+"humanGs_keys.json",'r') as f:
        HGs_keys = json.load(f)

if(test_tag == "OPA2Vec-Cosine"):#
    HGs_HDs_sim = np.load("Data/"+test_tag+"_HG_HD_sim.npy")
    with open("Data/"+test_tag+"H_HD_keys.json",'r') as f:
        HDs_keys = json.load(f)
    with open("Data/"+test_tag+"_HG_keys.json",'r') as f:
        HGs_keys = json.load(f)

if(test_tag == "OWL2Vec-Cosine"):#
    HGs_HDs_sim = np.load("Data/"+test_tag+"_HG_HD_sim.npy")
    with open("Data/"+test_tag+"H_HD_keys.json",'r') as f:
        HDs_keys = json.load(f)
    with open("Data/"+test_tag+"_HG_keys.json",'r') as f:
        HGs_keys = json.load(f)

if(test_tag == "DL2Vec-cosine"):
    HGs_HDs_sim = np.load("Data/"+test_tag+"_HG_HD_sim.npy")
    with open("Data/"+test_tag+"H_HD_keys.json",'r') as f:
        HDs_keys = json.load(f)
    with open("Data/"+test_tag+"_HG_keys.json",'r') as f:
        HGs_keys = json.load(f)

if(test_tag == "Resnik_similarity"):#
    HGs_HDs_sim = np.load("Data/"+test_tag+"_HG_HD_sim.npy")
    with open("Data/"+test_tag+"H_HD_keys.json",'r') as f:
        HDs_keys = json.load(f)
    with open("Data/"+test_tag+"_HG_keys.json",'r') as f:
        HGs_keys = json.load(f)

if(test_tag == "Resnik_IC"):#
    with open("Data/"+test_tag+".json",'r') as f:
        sim_dic = json.load(f)
    items =  list(sim_dic.items())
    HDs_keys = [x[0] for x in items]
    HGs_HDs_sim = np.array([x[1] for x in items])
    with open("Data/"+test_tag+"humanGs_keys.json",'r') as f:
        HGs_keys = json.load(f)




TPR,FPR,ROC,P,N = computeROC_AUC(HGs_HDs_sim,HDs_keys,HGs_keys,evaluation_set,m_positives, h_positives, "h","human")


with open('TPR_FPR_human_'+test_tag+"_"+sys.argv[1]+'.json', 'w') as fb:
    json.dump([TPR,FPR],fb)

print("All phenotypes results")
print("ROC_AUC = ",ROC)
print("CI = ",2*(math. sqrt(ROC*(1-ROC)/(min([P,N])))))
print("Precision@1 @10 @50 @100 = ",precisionAt(TPR,P,FPR,N,1),precisionAt(TPR,P,FPR,N,10),precisionAt(TPR,P,FPR,N,50),precisionAt(TPR,P,FPR,N,100))
print("Recall@1 @10 @50 @100 = ",recallAt(TPR,P,FPR,N,1),recallAt(TPR,P,FPR,N,10),recallAt(TPR,P,FPR,N,50),recallAt(TPR,P,FPR,N,100))
print("Hits@1 @10 @50 @100 = ",hitAt(TPR,P,FPR,N,1),hitAt(TPR,P,FPR,N,10),hitAt(TPR,P,FPR,N,50),hitAt(TPR,P,FPR,N,100))
print("Total positives = ", P)
    
ROC_ls.append(ROC)
CI_ls.append(2*(math. sqrt(ROC*(1-ROC)/(min([P,N])))))
Genes_ls.append(len(evaluation_set))
Associations_ls.append(P)
H_ls[0].append(hitAt(TPR,P,FPR,N,1))
H_ls[1].append(hitAt(TPR,P,FPR,N,10))
H_ls[2].append(hitAt(TPR,P,FPR,N,100))
P_ls[0].append(precisionAt(TPR,P,FPR,N,1))
P_ls[1].append(precisionAt(TPR,P,FPR,N,10))
P_ls[2].append(precisionAt(TPR,P,FPR,N,100))
R_ls[0].append(recallAt(TPR,P,FPR,N,1))
R_ls[1].append(recallAt(TPR,P,FPR,N,10))
R_ls[2].append(recallAt(TPR,P,FPR,N,100))


if(m):
    M_human_p = {}
    M_mouse_p = {}
    M_MGenes = []
    MGenes = []
    M_OMIM = []

    for hg in evaluation_set:
        if hg in m_human_mapping.keys():
            for g in m_human_mapping[hg]:
                MGenes.append(g)


    with open("Data/human_positives.tsv") as f:
        content = f.readlines()
    for p in content:
        hg = p.split()[1]
        if hg in evaluation_set:
            if hg in m_human_mapping.keys():
                for g in m_human_mapping[hg]:
                    #print('h_positive addition')
                    for om in p.split()[0].split('|'):
                        if(om in M_human_p):
                            M_human_p[om].append(g)
                            M_OMIM.append(om)
                        else:
                            M_human_p[om]=[g]
                            M_OMIM.append(om)

    MGs_keys = ""
    HDs_keys = ""
    MGs_HDs_sim = ""

    if(test_tag == "DL2Vec-PLR"):#
        with open("Data/"+test_tag+"mouseGs_keys.json",'r') as f:
            MGs_keys = json.load(f)

        with open("Data/"+test_tag+"mouse.json",'r') as f: #HDs_keys.json
            sim_dic = json.load(f)

        items =  list(sim_dic.items())
        HDs_keys = [x[0] for x in items]
        MGs_HDs_sim = np.array([x[1] for x in items])


    if(test_tag == "OPA2Vec-Cosine"):#
        MGs_HDs_sim = np.load("Data/"+test_tag+"_MGs_HDs_sim.npy")
        with open("Data/"+test_tag+"_HDs_keys.json",'r') as f:
            HDs_keys = json.load(f)
        with open("Data/"+test_tag+"_MGs_keys.json",'r') as f:
            MGs_keys = json.load(f)

    if(test_tag == "OWL2Vec-Cosine"):#
        MGs_HDs_sim = np.load("Data/"+test_tag+"_MGs_HDs_sim.npy")
        with open("Data/"+test_tag+"_HDs_keys.json",'r') as f:
            HDs_keys = json.load(f)
        with open("Data/"+test_tag+"_MGs_keys.json",'r') as f:
            MGs_keys = json.load(f)

    if(test_tag == "Resnik_similarity"):#
        MGs_HDs_sim = np.load("Data/"+test_tag+"_MGs_HDs_sim.npy")
        with open("Data/"+test_tag+"_HDs_keys.json",'r') as f:
            HDs_keys = json.load(f)
        with open("Data/"+test_tag+"_MGs_keys.json",'r') as f:
            MGs_keys = json.load(f)

    if(test_tag == "DL2Vec-cosine"):
        MGs_HDs_sim = np.load("Data/"+test_tag+"_MGs_HDs_sim.npy")
        with open("Data/"+test_tag+"_HDs_keys.json",'r') as f:
            HDs_keys = json.load(f)
        with open("Data/"+test_tag+"_MGs_keys.json",'r') as f:
            MGs_keys = json.load(f)

    if(test_tag == "Resnik_IC"):#
        with open("Data/"+test_tag+"mouseGs_keys.json",'r') as f:
            MGs_keys = json.load(f)

        with open("Data/"+test_tag+"mouse.json",'r') as f: #HDs_keys.json
            sim_dic = json.load(f)

        items =  list(sim_dic.items())
        HDs_keys = [x[0] for x in items]
        MGs_HDs_sim = np.array([x[1] for x in items])



    TPR,FPR,ROC,P,N =  computeROC_AUC(MGs_HDs_sim,HDs_keys,MGs_keys,MGenes, M_mouse_p, M_human_p,"h","mouse")
    
    #with open('TPR_FPR_Human_'+test_tag+"_"+sys.argv[1]+'.json', 'w') as fb:
    #    json.dump([TPR,FPR],fb)
    print("Mouse phenotypes results")
    print("Human Evaluation")
    print("ROC_AUC = ",ROC)
    print("CI = ",2*(math. sqrt(ROC*(1-ROC)/(min([P,N])))))
    print("Precision@1 @10 @50 @100 = ",precisionAt(TPR,P,FPR,N,1),precisionAt(TPR,P,FPR,N,10),precisionAt(TPR,P,FPR,N,50),precisionAt(TPR,P,FPR,N,100))
    print("Recall@1 @10 @50 @100 = ",recallAt(TPR,P,FPR,N,1),recallAt(TPR,P,FPR,N,10),recallAt(TPR,P,FPR,N,50),recallAt(TPR,P,FPR,N,100))
    print("Hits@1 @10 @50 @100 = ",hitAt(TPR,P,FPR,N,1),hitAt(TPR,P,FPR,N,10),hitAt(TPR,P,FPR,N,50),hitAt(TPR,P,FPR,N,100))
    print("Total positives = ", P)

    #plt.scatter(FPR,TPR)
    #plt.ylabel('TPR')
    #plt.xlabel('FPR')
    #plt.savefig('ROC_figures/'+sys.argv[4]+'.png')

    ROC_ls.append(ROC)
    CI_ls.append(2*(math. sqrt(ROC*(1-ROC)/(min([P,N])))))
    Genes_ls.append(len(MGenes))
    Associations_ls.append(P)
    H_ls[0].append(hitAt(TPR,P,FPR,N,1))
    H_ls[1].append(hitAt(TPR,P,FPR,N,10))
    H_ls[2].append(hitAt(TPR,P,FPR,N,100))
    P_ls[0].append(precisionAt(TPR,P,FPR,N,1))
    P_ls[1].append(precisionAt(TPR,P,FPR,N,10))
    P_ls[2].append(precisionAt(TPR,P,FPR,N,100))
    R_ls[0].append(recallAt(TPR,P,FPR,N,1))
    R_ls[1].append(recallAt(TPR,P,FPR,N,10))
    R_ls[2].append(recallAt(TPR,P,FPR,N,100))

'''
print("ROC_AUC = ")
print(roc_auc)
plt.plot(tpr,fpr)
plt.ylabel('TPR')
plt.xlabel('FPR')
plt.show()
'''

#----------------------------------------------------------
# zebrafish -> mouse -> human data

if(z):
    ZFGenes = []
    m_ZFGenes = []
    ZF_human = {}
    ZF_human_p = {}
    ZF_mouse = {}
    ZF_mouse_p = {}
    ZF_OMIM = []
   
    with open('Data/human_1_n_fish_mapping.json','r') as f:
        ZF_human = json.load(f)
    #with open('fish_human_mapping.json', 'w') as fp:
    #    json.dump(ZF_human, fp)
        
    

    #print(ZF_human)
    # getting human positive set
    #h_index = []

    for hg in evaluation_set:
        if hg in ZF_human.keys():
            for g in ZF_human[hg]:
                ZFGenes.append(g)


    with open("Data/human_positives.tsv") as f:
        content = f.readlines()
    for p in content:
        hg = p.split()[1]
        #print(hg)
        if hg in evaluation_set:
            if hg in ZF_human.keys():
                for g in ZF_human[hg]:
                    #print('h_positive addition')
                    #ZFGenes.append(g)
                    for om in p.split()[0].split('|'):
                        if(om in ZF_human_p):
                            ZF_human_p[om].append(g)
                            ZF_OMIM.append(om)
                        else:
                            ZF_human_p[om]=[g]
                            ZF_OMIM.append(om)


    ZFGs_keys = ""
    HDs_keys = ""
    ZFGs_HDs_sim = ""

    if(test_tag == "DL2Vec-PLR"):#
        with open("Data/"+test_tag+"fish.json",'r') as f: #HDs_keys.json
            sim_dic = json.load(f)

        items =  list(sim_dic.items())
        HDs_keys = [x[0] for x in items]
        ZFGs_HDs_sim = np.array([x[1] for x in items])
    
        with open("Data/"+test_tag+"fishGs_keys.json",'r') as f: #ZFGs_keys.json
            ZFGs_keys = json.load(f)

    if(test_tag == "OPA2Vec-Cosine"):#
        ZFGs_HDs_sim = np.load("Data/"+test_tag+"_ZFGs_HDs_sim.npy")
        with open("Data/"+test_tag+"_HDs_keys.json",'r') as f:
            HDs_keys = json.load(f)
        with open("Data/"+test_tag+"_ZFGs_keys.json",'r') as f:
            ZFGs_keys = json.load(f)
        
    if(test_tag == "OWL2Vec-Cosine"):#
        ZFGs_HDs_sim = np.load("Data/"+test_tag+"_ZFGs_HDs_sim.npy")
        with open("Data/"+test_tag+"_HDs_keys.json",'r') as f:
            HDs_keys = json.load(f)
        with open("Data/"+test_tag+"_ZFGs_keys.json",'r') as f:
            ZFGs_keys = json.load(f)


    if(test_tag == "Resnik_similarity"):#
        ZFGs_HDs_sim = np.load("Data/"+test_tag+"_ZFGs_HDs_sim.npy")
        with open("Data/"+test_tag+"_HDs_keys.json",'r') as f:
            HDs_keys = json.load(f)
        with open("Data/"+test_tag+"_ZFGs_keys.json",'r') as f:
            ZFGs_keys = json.load(f)

    if(test_tag == "DL2Vec-cosine"):
        ZFGs_HDs_sim = np.load("Data/"+test_tag+"_ZFGs_HDs_sim.npy")
        with open("Data/"+test_tag+"_HDs_keys.json",'r') as f:
            HDs_keys = json.load(f)
        with open("Data/"+test_tag+"_ZFGs_keys.json",'r') as f:
            ZFGs_keys = json.load(f)

    if(test_tag == "Resnik_IC"):#
        with open("Data/"+test_tag+"fish.json",'r') as f: #HDs_keys.json
            sim_dic = json.load(f)

        items =  list(sim_dic.items())
        HDs_keys = [x[0] for x in items]
        ZFGs_HDs_sim = np.array([x[1] for x in items])
    
        with open("Data/"+test_tag+"fishGs_keys.json",'r') as f: #ZFGs_keys.json
            ZFGs_keys = json.load(f)

   
    
    TPR,FPR,ROC,P,N = computeROC_AUC(ZFGs_HDs_sim,HDs_keys,ZFGs_keys,m_ZFGenes,ZF_mouse_p, ZF_human_p,"o","fish")


    m_ROC_ls.append(ROC)
    m_CI_ls.append(2*(math. sqrt(ROC*(1-ROC)/(min([P,N])))))
    m_Genes_ls.append(len(ZFGenes))
    m_Associations_ls.append(P)
    m_H_ls[0].append(hitAt(TPR,P,FPR,N,1))
    m_H_ls[1].append(hitAt(TPR,P,FPR,N,10))
    m_H_ls[2].append(hitAt(TPR,P,FPR,N,100))
    m_P_ls[0].append(precisionAt(TPR,P,FPR,N,1))
    m_P_ls[1].append(precisionAt(TPR,P,FPR,N,10))
    m_P_ls[2].append(precisionAt(TPR,P,FPR,N,100))
    m_R_ls[0].append(recallAt(TPR,P,FPR,N,1))
    m_R_ls[1].append(recallAt(TPR,P,FPR,N,10))
    m_R_ls[2].append(recallAt(TPR,P,FPR,N,100))   
  
  
    TPR,FPR,ROC,P,N = computeROC_AUC(ZFGs_HDs_sim,HDs_keys,ZFGs_keys,ZFGenes,ZF_mouse_p, ZF_human_p,"h","fish")
    '''
    with open('TPR_FPR_Human_'+test_tag+"_"+sys.argv[1]+'.json', 'w') as fb:
        json.dump([TPR,FPR],fb)
    '''

    print("Fish phjenotypes results")
    print("Human Evaluation")
    print("ROC_AUC = ",ROC)
    print("CI = ",2*(math. sqrt(ROC*(1-ROC)/(min([P,N])))))
    print("Precision@1 @10 @50 @100 = ",precisionAt(TPR,P,FPR,N,1),precisionAt(TPR,P,FPR,N,10),precisionAt(TPR,P,FPR,N,50),precisionAt(TPR,P,FPR,N,100))
    print("Recall@1 @10 @50 @100 = ",recallAt(TPR,P,FPR,N,1),recallAt(TPR,P,FPR,N,10),recallAt(TPR,P,FPR,N,50),recallAt(TPR,P,FPR,N,100))
    print("Hits@1 @10 @50 @100 = ",hitAt(TPR,P,FPR,N,1),hitAt(TPR,P,FPR,N,10),hitAt(TPR,P,FPR,N,50),hitAt(TPR,P,FPR,N,100))
    print("Total positives = ", P)
    

    ROC_ls.append(ROC)
    CI_ls.append(2*(math. sqrt(ROC*(1-ROC)/(min([P,N])))))
    Genes_ls.append(len(ZFGenes))
    Associations_ls.append(P)
    H_ls[0].append(hitAt(TPR,P,FPR,N,1))
    H_ls[1].append(hitAt(TPR,P,FPR,N,10))
    H_ls[2].append(hitAt(TPR,P,FPR,N,100))
    P_ls[0].append(precisionAt(TPR,P,FPR,N,1))
    P_ls[1].append(precisionAt(TPR,P,FPR,N,10))
    P_ls[2].append(precisionAt(TPR,P,FPR,N,100))
    R_ls[0].append(recallAt(TPR,P,FPR,N,1))
    R_ls[1].append(recallAt(TPR,P,FPR,N,10))
    R_ls[2].append(recallAt(TPR,P,FPR,N,100))   


#fly
if(fly):
    FGenes = []
    F_human = {}
    F_human_p = {}
    F_mouse_p = {}
    F_OMIM = []
    m_FGenes = []


    FGs_keys = ""
    HDs_keys = ""
    FGs_HDs_sim = ""

    if(test_tag == "DL2Vec-PLR"):#
        with open("Data/"+test_tag+"fly.json",'r') as f: #HDs_keys.json
            sim_dic = json.load(f)

        items =  list(sim_dic.items())
        HDs_keys = [x[0] for x in items]
        FGs_HDs_sim = np.array([x[1] for x in items])


        with open("Data/"+test_tag+"flyGs_keys.json",'r') as f: #FGs_keys.json
            FGs_keys = json.load(f)

    if(test_tag == "OPA2Vec-Cosine"):#
        FGs_HDs_sim = np.load("Data/"+test_tag+"_FBGs_HDs_sim.npy")
        with open("Data/"+test_tag+"_HDs_keys.json",'r') as f:
            HDs_keys = json.load(f)
        with open("Data/"+test_tag+"_FBGs_keys.json",'r') as f:
            FGs_keys = json.load(f)


    
    if(test_tag == "OWL2Vec-Cosine"):#
        FGs_HDs_sim = np.load("Data/"+test_tag+"_FBGs_HDs_sim.npy")
        with open("Data/"+test_tag+"_HDs_keys.json",'r') as f:
            HDs_keys = json.load(f)
        with open("Data/"+test_tag+"_FBGs_keys.json",'r') as f:
            FGs_keys = json.load(f)

    if(test_tag == "Resnik_similarity"):#
        FGs_HDs_sim = np.load("Data/"+test_tag+"_FBGs_HDs_sim.npy")
        with open("Data/"+test_tag+"_HDs_keys.json",'r') as f:
            HDs_keys = json.load(f)
        with open("Data/"+test_tag+"_FBGs_keys.json",'r') as f:
            FGs_keys = json.load(f)

    if(test_tag == "DL2Vec-cosine"):
        FGs_HDs_sim = np.load("Data/"+test_tag+"_FBGs_HDs_sim.npy")
        with open("Data/"+test_tag+"_HDs_keys.json",'r') as f:
            HDs_keys = json.load(f)
        with open("Data/"+test_tag+"_FBGs_keys.json",'r') as f:
            FGs_keys = json.load(f)

    if(test_tag == "Resnik_IC"):#
        with open("Data/"+test_tag+"fly.json",'r') as f: #HDs_keys.json
            sim_dic = json.load(f)

        items =  list(sim_dic.items())
        HDs_keys = [x[0] for x in items]
        FGs_HDs_sim = np.array([x[1] for x in items])


        with open("Data/"+test_tag+"flyGs_keys.json",'r') as f: #FGs_keys.json
            FGs_keys = json.load(f)

    with open("Data/fly_human_mapping.json",'r') as f:
        F_human = json.load(f)

    for hg in evaluation_set:
        if hg in F_human.keys():
            for g in F_human[hg]:
                FGenes.append(g)
        else:
            print("not found",hg)

    with open("Data/human_positives.tsv") as f:
        content = f.readlines()
    for p in content:
        hg = p.split()[1]
        #print(hg)
        if hg in evaluation_set:
            if hg in F_human.keys():
                for g in F_human[hg]:
                    #FGenes.append(g)
                    for o in p.split()[0].split('|'):
                        om=o.strip()
                        if(om in F_human_p):
                            F_human_p[om].append(g)
                            F_OMIM.append(om)
                        else:
                            F_human_p[om]=[g]
                            F_OMIM.append(om)
    print("disease here Are", len(F_OMIM), len(list(set(F_OMIM))))

    
    #print(F_human_p)
    m_F_OMIM = []
    for om,mgList in m_positives.items():
        for mg in mgList:
            if(mg in F_mouse.keys()):
                for g in F_mouse[mg]:
                    if om in F_mouse_p:
                        F_mouse_p[om].append(g)                        
                    else:
                        F_mouse_p[om]=[g]
                        m_F_OMIM.append(om)
    

    print("genes",len(list(set(m_FGenes))),"diseases" ,len(list(set(m_F_OMIM))))

    

    TPR,FPR,ROC,P,N = computeROC_AUC(FGs_HDs_sim,HDs_keys,FGs_keys,FGenes,F_mouse_p, F_human_p,"h","fly")
    
    #with open('TPR_FPR_Human_'+test_tag+"_"+sys.argv[1]+'.json', 'w') as fb:
        #json.dump([TPR,FPR],fb)
    print("Fly phenotypes")
    print("Human Evaluation")
    print("ROC_AUC = ",ROC)
    print("CI = ",2*(math. sqrt(ROC*(1-ROC)/(min([P,N])))))
    print("Precision@1 @10 @50 @100 = ",precisionAt(TPR,P,FPR,N,1),precisionAt(TPR,P,FPR,N,10),precisionAt(TPR,P,FPR,N,50),precisionAt(TPR,P,FPR,N,100))
    print("Recall@1 @10 @50 @100 = ",recallAt(TPR,P,FPR,N,1),recallAt(TPR,P,FPR,N,10),recallAt(TPR,P,FPR,N,50),recallAt(TPR,P,FPR,N,100))
    print("Hits@1 @10 @50 @100 = ",hitAt(TPR,P,FPR,N,1),hitAt(TPR,P,FPR,N,10),hitAt(TPR,P,FPR,N,50),hitAt(TPR,P,FPR,N,100))
    print("Total positives = ", P)

    ROC_ls.append(ROC)
    CI_ls.append(2*(math. sqrt(ROC*(1-ROC)/(min([P,N])))))
    Genes_ls.append(len(FGenes))
    Associations_ls.append(P)
    H_ls[0].append(hitAt(TPR,P,FPR,N,1))
    H_ls[1].append(hitAt(TPR,P,FPR,N,10))
    H_ls[2].append(hitAt(TPR,P,FPR,N,100))
    P_ls[0].append(precisionAt(TPR,P,FPR,N,1))
    P_ls[1].append(precisionAt(TPR,P,FPR,N,10))
    P_ls[2].append(precisionAt(TPR,P,FPR,N,100))
    R_ls[0].append(recallAt(TPR,P,FPR,N,1))
    R_ls[1].append(recallAt(TPR,P,FPR,N,10))
    R_ls[2].append(recallAt(TPR,P,FPR,N,100))


    #----------------------------------------------------------
#yeast
if(y):
    FGenes = []
    F_human = {}
    F_human_p = {}
    F_mouse = Y_mouse
    F_mouse_p = {}
    m_FGenes = []
    F_OMIM = []

    

    FGs_keys = ""
    HDs_keys = ""
    FGs_HDs_sim = ""



    if(test_tag == "DL2Vec-PLR"):#
        with open("Data/"+test_tag+"yeast.json",'r') as f: #HDs_keys.json
            sim_dic = json.load(f)

        items =  list(sim_dic.items())
        HDs_keys = [x[0] for x in items]
        FGs_HDs_sim = np.array([x[1] for x in items])


        with open("Data/"+test_tag+"yeastGs_keys.json",'r') as f: #FGs_keys.json
            FGs_keys = json.load(f)

    if(test_tag == "OPA2Vec-Cosine"):#
        FGs_HDs_sim = np.load("Data/"+test_tag+"_YGs_HDs_sim.npy")
        with open("Data/"+test_tag+"_HDs_keys.json",'r') as f:
            HDs_keys = json.load(f)
        with open("Data/"+test_tag+"_YGs_keys.json",'r') as f:
            FGs_keys = json.load(f)
            

    if(test_tag == "OWL2Vec-Cosine"):#
        FGs_HDs_sim = np.load("Data/"+test_tag+"_YGs_HDs_sim.npy")
        with open("Data/"+test_tag+"_HDs_keys.json",'r') as f:
            HDs_keys = json.load(f)
        with open("Data/"+test_tag+"_YGs_keys.json",'r') as f:
            FGs_keys = json.load(f)

    if(test_tag == "Resnik_similarity"):#
        FGs_HDs_sim = np.load("Data/"+test_tag+"_YGs_HDs_sim.npy")
        with open("Data/"+test_tag+"_HDs_keys.json",'r') as f:
            HDs_keys = json.load(f)
        with open("Data/"+test_tag+"_YGs_keys.json",'r') as f:
            FGs_keys = json.load(f)


    if(test_tag == "DL2Vec-cosine"): 
        FGs_HDs_sim = np.load("Data/"+test_tag+"_YGs_HDs_sim.npy")
        with open("Data/"+test_tag+"_HDs_keys.json",'r') as f:
            HDs_keys = json.load(f)
        with open("Data/"+test_tag+"_YGs_keys.json",'r') as f:
            FGs_keys = json.load(f)


    if(test_tag == "Resnik_IC"):#
    if(test_tag == "DL2Vec-PLR"):#
        with open("Data/"+test_tag+"yeast.json",'r') as f: #HDs_keys.json
            sim_dic = json.load(f)

        items =  list(sim_dic.items())
        HDs_keys = [x[0] for x in items]
        FGs_HDs_sim = np.array([x[1] for x in items])


        with open("Data/"+test_tag+"yeastGs_keys.json",'r') as f: #FGs_keys.json
            FGs_keys = json.load(f)


    with open("Data/yeast_human_mapping.json",'r') as f:
        F_human = json.load(f)

    with open("Data/yeast_mouse_mapping.json",'r') as f:
        F_mouse = json.load(f)

    for hg in evaluation_set:
        if hg in F_human.keys():
            for g in F_human[hg]:
                FGenes.append(g)

    with open("Data/human_positives.tsv") as f:
        content = f.readlines()
    for p in content:
        hg = p.split()[1]
        #print(hg)
        if hg in evaluation_set:
            if hg in F_human.keys():
                for g in F_human[hg]:
                    #FGenes.append(g)
                    for om in p.split()[0].split('|'):
                        if(om in F_human_p):
                            F_human_p[om].append(g)
                            F_OMIM.append(om)
                        else:
                            F_human_p[om]=[g]
                            F_OMIM.append(om)



    TPR,FPR,ROC,P,N = computeROC_AUC(FGs_HDs_sim,HDs_keys,FGs_keys,FGenes,F_mouse_p, F_human_p,"h","yeast")
    
    #with open('TPR_FPR_Human_'+test_tag+"_"+sys.argv[1]+'.json', 'w') as fb:
        #json.dump([TPR,FPR],fb)
    print("Yeast phenotypes")
    print("Human Evaluation")
    print("ROC_AUC = ",ROC)
    print("CI = ",2*(math. sqrt(ROC*(1-ROC)/(min([P,N])))))
    print("Precision@1 @10 @50 @100 = ",precisionAt(TPR,P,FPR,N,1),precisionAt(TPR,P,FPR,N,10),precisionAt(TPR,P,FPR,N,50),precisionAt(TPR,P,FPR,N,100))
    print("Recall@1 @10 @50 @100 = ",recallAt(TPR,P,FPR,N,1),recallAt(TPR,P,FPR,N,10),recallAt(TPR,P,FPR,N,50),recallAt(TPR,P,FPR,N,100))
    print("Hits@1 @10 @50 @100 = ",hitAt(TPR,P,FPR,N,1),hitAt(TPR,P,FPR,N,10),hitAt(TPR,P,FPR,N,50),hitAt(TPR,P,FPR,N,100))
    print("Total positives = ", P)


    ROC_ls.append(ROC)
    CI_ls.append(2*(math. sqrt(ROC*(1-ROC)/(min([P,N])))))
    Genes_ls.append(len(FGenes))
    Associations_ls.append(P)
    H_ls[0].append(hitAt(TPR,P,FPR,N,1))
    H_ls[1].append(hitAt(TPR,P,FPR,N,10))
    H_ls[2].append(hitAt(TPR,P,FPR,N,100))
    P_ls[0].append(precisionAt(TPR,P,FPR,N,1))
    P_ls[1].append(precisionAt(TPR,P,FPR,N,10))
    P_ls[2].append(precisionAt(TPR,P,FPR,N,100))
    R_ls[0].append(recallAt(TPR,P,FPR,N,1))
    R_ls[1].append(recallAt(TPR,P,FPR,N,10))
    R_ls[2].append(recallAt(TPR,P,FPR,N,100))



#-Summary---------------------


print("------------------------(human evaluation data set )------------------------------")

print("ROC_AUC\t" + "\t".join(map(str, ROC_ls)))
print("CI\t"+"\t".join(map(str, CI_ls)) )

ls = []
for i in range(len(ROC_ls)):
    string_result = "$"+str(float("{:.3f}".format(ROC_ls[i])))+"\pm"+str(float("{:.3f}".format(CI_ls[i])))+"$"
    ls.append(string_result)
print("\t".join(ls))

ls = []
for i in range(len(ROC_ls)):
    string_result = str(float("{:.3e}".format(ROC_ls[i])))+"\pm"+str(float("{:.3e}".format(CI_ls[i])))
    ls.append(string_result)
print("\t".join(ls))

print("ganes\t" + "\t".join(map(str,Genes_ls)))

print("Associations\t" + "\t".join(map(str, Associations_ls)))

print("Hits@1\t" + "\t".join(map(str, H_ls[0])))
print("Hits@10\t" + "\t".join(map(str, H_ls[1])))
print("Hits@100\t" + "\t".join(map(str, H_ls[2])))

print("Precision@1\t" + "\t".join(map(str, P_ls[0])))
print("Precision@10\t" + "\t".join(map(str, P_ls[1])))
print("Precision@100\t" + "\t".join(map(str, P_ls[2])))

print("Recall@1\t" + "\t".join(map(str, R_ls[0])))
print("Recall@10\t" + "\t".join(map(str, R_ls[1])))
print("Recall@100\t" + "\t".join(map(str, R_ls[2])))

