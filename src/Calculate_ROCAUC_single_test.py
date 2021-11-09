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
test_tag = "resnik"
#test_tag = "resnik_IC"
#test_tag2 = "human_pheno_e"
#test_tag = "yeast_no_pheno_OPA"
#test_tag = "OPA_cosine"
#test_tag = "DL2Vec_cosine"
#test_tag = "OWL2Vec_cosine"
test_tag2 = "Pheno-e_"+sys.argv[1]
#test_tag2 = "uPheno"
#test_tag="DL2Vec"
#test_tag2="_1_balanced_model_upheno_"


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

print(test_tag,test_tag2)

def computeROC_AUC(OGs_HDs_sim,HDs_keys,OGs_keys,o_positives,h_positives,o_or_h):
    OGs_HDs_sim_t = OGs_HDs_sim.transpose()
    if(test_tag=="DL2Vec" or test_tag=="resnik_IC"):
        OGs_HDs_sim_t = OGs_HDs_sim
    ranks = np.argsort(OGs_HDs_sim_t, axis=1)
    TPR = [0]
    FPR = [0]
    prev= [0,0]
    P=0
    N=0
    positive_genes = set([])
    #print(h_positives)
    includedDiseases =  np.zeros(len(HDs_keys))
    print(OGs_HDs_sim_t.shape,len(HDs_keys),len(OGs_keys))

    positives_matrix = np.zeros([len(HDs_keys),len(OGs_keys)])
    for og in range(0,len(OGs_keys)):
        for hd in range(0,len(HDs_keys)):
            if o_or_h == 'o':
                if(HDs_keys[hd] in  o_positives):
                    if(OGs_keys[og] in o_positives[HDs_keys[hd]]):
                        positives_matrix[hd][og] = 1
                        includedDiseases[hd]=1
            if o_or_h == 'h':
                if(HDs_keys[hd] in  h_positives):
                    if(OGs_keys[og] in h_positives[HDs_keys[hd]]):
                        positives_matrix[hd][og] = 1
                        includedDiseases[hd]=1
                        positive_genes.add(OGs_keys[og])
    ranking_dic_disease={}
    ranking_dic_gene={}
    positives_ranks = []
    #np.save("tests_inputs/for_prof/similarity_matrex.npy",OGs_HDs_sim_t)
    #np.save("tests_inputs/for_prof/positives_matrex.npy",positives_matrix)
    #with open('tests_inputs/for_prof/h_positives.json', 'w') as fb:
    #    json.dump(h_positives,fb)

    for hd in range(0,len(HDs_keys)):
        if(includedDiseases[hd]==1):
            p = np.sum(positives_matrix[hd])
            P+= p
            N+= len(OGs_keys)-p
    print("p",P)
    #print("positive_genes",positive_genes)
    print("included_Disease", np.sum(includedDiseases))
    # r is the rank considered
    for r in range(0,len(OGs_keys)):
        #print("rank",r)
        TP = prev[0]
        FP = prev[1]
        for hd in range(0,len(HDs_keys)):
            if(includedDiseases[hd]==1):
                #loop throgh diseases
                g=len(OGs_keys)-r-1
                if o_or_h == 'o':
                    if (positives_matrix[hd,ranks[hd][g]] == 1):
                        TP+=1
                        positives_ranks.append(r+1)
                        if(HDs_keys[hd] not in ranking_dic_disease):
                            ranking_dic_disease[HDs_keys[hd]] = r+1
                            if (OGs_keys[g] not in ranking_dic_gene):
                                ranking_dic_gene[OGs_keys[g]]=[]
                            ranking_dic_gene[OGs_keys[g]].append(r+1)
                    else:
                        FP+=1
                elif o_or_h == 'h':
                    if (positives_matrix[hd,ranks[hd][g]] == 1):
                        TP+=1
                        positives_ranks.append(r+1)
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
        if(math.isnan(TP/N) or math.isnan(FP/N)):
            print(TP,FP,N,r,d,og,OGs_keys[og])
            return
        #print(r,TP/P,FP/N)

    with open('positive_ranks_'+test_tag+'_'+test_tag2+'.json', 'w') as fb:
        json.dump(positives_ranks,fb)    
    '''
    with open('ranking_dic_disease_'+test_tag+'_'+sys.argv[1]+test_tag2+'.json', 'w') as fb:
        json.dump(ranking_dic_disease,fb)
    with open('ranking_dic_genes_'+test_tag+'_'+sys.argv[1]+test_tag2+'.json', 'w') as fb:
        json.dump(ranking_dic_gene,fb)
    '''
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
if(sys.argv[1]=='h'):
    h = True
elif(sys.argv[1]=='m'):
    m = True
elif(sys.argv[1]=='z'):
    z = True
elif(sys.argv[1]=='y'):
    y = True
elif(sys.argv[1]=='f'):
    fly = True
else:
    print("please specify an organism")
'''
with open('test_set_for_human_all.json', 'r') as fb:
    human_eval = json.load(fb)

with open('test_set_for_mouse_all.json', 'r') as fb:
    mouse_eval = json.load(fb)
'''
#---------------------------------------------------------

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
    #print(hg)

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


# getting mice positive set
count=0
OMIMs=[]
with open("Data/mouse_positives.tsv") as f:
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


if(h):
    HGs_keys = ""
    HDs_keys = ""
    HGs_HDs_sim = ""

    HGs_HDs_sim = np.load(sys.argv[2])
    with open(sys.argv[3],'r') as f:
        HDs_keys = json.load(f)
    with open(sys.argv[4],'r') as f:
        HGs_keys = json.load(f)
    TPR,FPR,ROC,P,N = computeROC_AUC(HGs_HDs_sim,HDs_keys,HGs_keys,m_positives, h_positives, "h")
    '''
    with open('TPR_FPR_human_'+test_tag+"_"+sys.argv[1]+'.json', 'w') as fb:
        json.dump([TPR,FPR],fb)
    '''
    CI = 2*(math. sqrt(ROC*(1-ROC)/(min([P,N]))))
    print("ROC_AUC = ",ROC)
    print("CI = ",CI)
    print(str(float("{:.3f}".format(ROC)))+"\mp"+str(float("{:.3f}".format(CI))) )
    print("Precision@1 @10 @50 @100 = ",precisionAt(TPR,P,FPR,N,1),precisionAt(TPR,P,FPR,N,10),precisionAt(TPR,P,FPR,N,50),precisionAt(TPR,P,FPR,N,100))
    print("Recall@1 @10 @50 @100 = ",recallAt(TPR,P,FPR,N,1),recallAt(TPR,P,FPR,N,10),recallAt(TPR,P,FPR,N,50),recallAt(TPR,P,FPR,N,100))
    print("Hits@1 @10 @50 @100 = ",hitAt(TPR,P,FPR,N,1),hitAt(TPR,P,FPR,N,10),hitAt(TPR,P,FPR,N,50),hitAt(TPR,P,FPR,N,100))
    print("Total positives = ", P)
    



if(m):
    M_human_p = {}
    MGenes = []
    M_OMIM = []
    with open("Data/human_positives.tsv") as f:
        content = f.readlines()
    for p in content:
        hg = p.split()[1]
        #print(hg)
        if hg in m_human_mapping.keys():
            for g in m_human_mapping[hg]:
                #print('h_positive addition')
                MGenes.append(m_human_mapping[hg])
                for om in p.split()[0].split('|'):
                    if(om in M_human_p):
                        M_human_p[om].add(g)
                        M_OMIM.append(om)
                    else:
                        M_human_p[om]=set([g])
                        M_OMIM.append(om)
    positives = {}
    for om in M_human_p:
        positives[om] = list(M_human_p[om])
    #with open('m_positives.json', 'w') as fb:
    #    json.dump(positives,fb)

    MGs_keys = ""
    HDs_keys = ""
    MGs_HDs_sim = ""

    MGs_HDs_sim = np.load(sys.argv[2])
    with open(sys.argv[3],'r') as f:
        HDs_keys = json.load(f)
    with open(sys.argv[4],'r') as f:
        MGs_keys = json.load(f)

    

    TPR,FPR,ROC,P,N = computeROC_AUC(MGs_HDs_sim,HDs_keys,MGs_keys,m_positives, h_positives,"o")
    
       
    #with open('TPR_FPR_Mouse_'+test_tag+"_"+sys.argv[1]+'.json', 'w') as fb:
    #    json.dump([TPR,FPR],fb)
    print("Mouse Evaluation") 
    CI = 2*(math. sqrt(ROC*(1-ROC)/(min([P,N]))))
    print("ROC_AUC = ",ROC)
    print("CI = ",CI)
    print(str(float("{:.3f}".format(ROC)))+"\mp"+str(float("{:.3f}".format(CI))) )
    print("Precision@1 @10 @50 @100 = ",precisionAt(TPR,P,FPR,N,1),precisionAt(TPR,P,FPR,N,10),precisionAt(TPR,P,FPR,N,50),precisionAt(TPR,P,FPR,N,100))
    print("Recall@1 @10 @50 @100 = ",recallAt(TPR,P,FPR,N,1),recallAt(TPR,P,FPR,N,10),recallAt(TPR,P,FPR,N,50),recallAt(TPR,P,FPR,N,100))
    print("Hits@1 @10 @50 @100 = ",hitAt(TPR,P,FPR,N,1),hitAt(TPR,P,FPR,N,10),hitAt(TPR,P,FPR,N,50),hitAt(TPR,P,FPR,N,100))
    print("Total positives = ", P)

    

    TPR,FPR,ROC,P,N =  computeROC_AUC(MGs_HDs_sim,HDs_keys,MGs_keys,m_positives, M_human_p,"h")
    
    #with open('TPR_FPR_Human_'+test_tag+"_"+sys.argv[1]+'.json', 'w') as fb:
    #    json.dump([TPR,FPR],fb)
    print("Human Evaluation")
    CI = 2*(math. sqrt(ROC*(1-ROC)/(min([P,N]))))
    print("ROC_AUC = ",ROC)
    print("CI = ",CI)
    print(str(float("{:.3f}".format(ROC)))+"\mp"+str(float("{:.3f}".format(CI))) )    
    print("Precision@1 @10 @50 @100 = ",precisionAt(TPR,P,FPR,N,1),precisionAt(TPR,P,FPR,N,10),precisionAt(TPR,P,FPR,N,50),precisionAt(TPR,P,FPR,N,100))
    print("Recall@1 @10 @50 @100 = ",recallAt(TPR,P,FPR,N,1),recallAt(TPR,P,FPR,N,10),recallAt(TPR,P,FPR,N,50),recallAt(TPR,P,FPR,N,100))
    print("Hits@1 @10 @50 @100 = ",hitAt(TPR,P,FPR,N,1),hitAt(TPR,P,FPR,N,10),hitAt(TPR,P,FPR,N,50),hitAt(TPR,P,FPR,N,100))
    print("Total positives = ", P)

    #plt.scatter(FPR,TPR)
    #plt.ylabel('TPR')
    #plt.xlabel('FPR')
    #plt.savefig('ROC_figures/'+sys.argv[4]+'.png')


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
    ZF_human = {}
    ZF_human_p = {}
    ZF_mouse = {}
    ZF_mouse_p = {}
    ZF_OMIM = []

    with open('Data/human_1_n_fish_mapping.json','r') as f:
        ZF_human = json.load(f)
  
    with open("Data/mouse_1_n_fish_mapping.json",'r') as f:
        ZF_mouse = json.load(f)


    with open("Data/human_positives.tsv") as f:
        content = f.readlines()
    for p in content:
        hg = p.split()[1]
        if hg in ZF_human.keys():
            for g in ZF_human[hg]:
                ZFGenes.append(ZF_human[hg])
                for om in p.split()[0].split('|'):
                    if(om in ZF_human_p):
                        ZF_human_p[om].add(g)
                        ZF_OMIM.append(om)
                    else:
                        ZF_human_p[om]=set([g])
                        ZF_OMIM.append(om)
    positives = {}
    for om in ZF_human_p:
        positives[om] = list(ZF_human_p[om])
    #with open('z_positives.json', 'w') as fb:
    #    json.dump(positives,fb)

    for om,mgList in m_positives.items():
        for mg in mgList:
            if(mg in ZF_mouse.keys()):
                for g in ZF_mouse[mg]:
                    ZFGenes.append(ZF_mouse[mg])
                    if om in ZF_mouse_p:
                        ZF_mouse_p[om].append(g)
                        ZF_OMIM.append(om)
                    else:
                        ZF_mouse_p[om]=[g]
                        ZF_OMIM.append(om)
    ZFGs_keys = ""
    HDs_keys = ""
    ZFGs_HDs_sim = ""


    ZFGs_HDs_sim = np.load(sys.argv[2]) #ZFGs_HDs_sim.npy
    with open(sys.argv[3],'r') as f: #HDs_keys.json
        HDs_keys = json.load(f)
    with open(sys.argv[4],'r') as f: #ZFGs_keys.json
        ZFGs_keys = json.load(f)

    
    TPR,FPR,ROC,P,N = computeROC_AUC(ZFGs_HDs_sim,HDs_keys,ZFGs_keys,ZF_mouse_p, ZF_human_p,"o")


    #with open('TPR_FPR_Mouse_'+test_tag+"_"+sys.argv[1]+'.json', 'w') as fb:
    #    json.dump([TPR,FPR],fb)
    print("Mouse Evaluation")
    CI = 2*(math. sqrt(ROC*(1-ROC)/(min([P,N]))))
    print("ROC_AUC = ",ROC)
    print("CI = ",CI)
    print(str(float("{:.3f}".format(ROC)))+"\mp"+str(float("{:.3f}".format(CI))) )
    print("Precision@1 @10 @50 @100 = ",precisionAt(TPR,P,FPR,N,1),precisionAt(TPR,P,FPR,N,10),precisionAt(TPR,P,FPR,N,50),precisionAt(TPR,P,FPR,N,100))
    print("Recall@1 @10 @50 @100 = ",recallAt(TPR,P,FPR,N,1),recallAt(TPR,P,FPR,N,10),recallAt(TPR,P,FPR,N,50),recallAt(TPR,P,FPR,N,100))
    print("Hits@1 @10 @50 @100 = ",hitAt(TPR,P,FPR,N,1),hitAt(TPR,P,FPR,N,10),hitAt(TPR,P,FPR,N,50),hitAt(TPR,P,FPR,N,100))
    print("Total positives = ", P)
    

    TPR,FPR,ROC,P,N = computeROC_AUC(ZFGs_HDs_sim,HDs_keys,ZFGs_keys,ZF_mouse_p, ZF_human_p,"h")
    
    #with open('TPR_FPR_Human_'+test_tag+"_"+sys.argv[1]+'.json', 'w') as fb:
    #    json.dump([TPR,FPR],fb)
    print("Human Evaluation")
    CI = 2*(math. sqrt(ROC*(1-ROC)/(min([P,N]))))
    print("ROC_AUC = ",ROC)
    print("CI = ",CI)
    print(str(float("{:.3f}".format(ROC)))+"\mp"+str(float("{:.3f}".format(CI))) )
    print("Precision@1 @10 @50 @100 = ",precisionAt(TPR,P,FPR,N,1),precisionAt(TPR,P,FPR,N,10),precisionAt(TPR,P,FPR,N,50),precisionAt(TPR,P,FPR,N,100))
    print("Recall@1 @10 @50 @100 = ",recallAt(TPR,P,FPR,N,1),recallAt(TPR,P,FPR,N,10),recallAt(TPR,P,FPR,N,50),recallAt(TPR,P,FPR,N,100))
    print("Hits@1 @10 @50 @100 = ",hitAt(TPR,P,FPR,N,1),hitAt(TPR,P,FPR,N,10),hitAt(TPR,P,FPR,N,50),hitAt(TPR,P,FPR,N,100))
    print("Total positives = ", P)
    

#----------------------------------------------------------
#fly and yeast
if(fly or y):
    FGenes = []
    F_human = {}
    F_human_p = {}
    F_mouse = {}
    F_mouse_p = {}
    F_OMIM = []

    FGs_keys = ""
    HDs_keys = ""
    FGs_HDs_sim = ""
    
    if(fly):
        with open('Data/fly_human_mapping.json','r') as f:
            F_human = json.load(f)

        with open("Data/fly_mouse_mapping.json",'r') as f:
            F_mouse = json.load(f)
    if(y):
        with open("Data/yeast_human_mapping.json",'r') as f:
            F_human = json.load(f)

        with open("Data/yeast_mouse_mapping.json",'r') as f:
            F_mouse = json.load(f)





    

    with open("Data/human_positives.tsv") as f:
        content = f.readlines()
    for p in content:
        hg = p.split()[1]
        #print(hg)
        if hg in F_human.keys():
            for g in F_human[hg]:
                FGenes.append(g)
                for om in p.split()[0].split('|'):
                    if(om in F_human_p):
                        F_human_p[om].add(g)
                        F_OMIM.append(om)
                    else:
                        F_human_p[om]=set([g])
                        F_OMIM.append(om)
    positives = {}
    for om in F_human_p:
        positives[om] = list(F_human_p[om])




    FGs_HDs_sim = np.load(sys.argv[2]) #FGs_HDs_sim.npy
    with open(sys.argv[3],'r') as f: #HDs_keys.json
        HDs_keys = json.load(f)
    with open(sys.argv[4],'r') as f: #FGs_keys.json
        FGs_keys = json.load(f)


    #print(F_human_p)
    countf = 0
    countn = 0
    for om,mgList in m_positives.items():
        for mg in mgList:
            if(mg in F_mouse.keys()):
                #print("here!")
                for g in F_mouse[mg]:
                    if(g in FGs_keys):
                        countf+=1
                    else:
                        #print(g)
                        countn+=1
                    #print(g)
                    FGenes.append(g)
                    if om in F_mouse_p:
                        F_mouse_p[om].append(g)
                        
                        F_OMIM.append(om)
                        #m_index.append(p)
                    else:
                        F_mouse_p[om]=[g]
                        F_OMIM.append(om)
    
    print('mouse positive = ', len(F_mouse_p))
    #print(F_mouse_p)
    print('found ',countf ,' not found ', countn)
    
   
    
    TPR,FPR,ROC,P,N = computeROC_AUC(FGs_HDs_sim,HDs_keys,FGs_keys,F_mouse_p, F_human_p,"o")

    #with open('TPR_FPR_Mouse_'+test_tag+"_"+sys.argv[1]+'.json', 'w') as fb:
    #   json.dump([TPR,FPR],fb)
    print("Mouse Evaluation")
    CI = 2*(math. sqrt(ROC*(1-ROC)/(min([P,N]))))
    print("ROC_AUC = ",ROC)
    print("CI = ",CI)
    print(str(float("{:.3f}".format(ROC)))+"\mp"+str(float("{:.3f}".format(CI))) )
    print("Precision@1 @10 @50 @100 = ",precisionAt(TPR,P,FPR,N,1),precisionAt(TPR,P,FPR,N,10),precisionAt(TPR,P,FPR,N,50),precisionAt(TPR,P,FPR,N,100))
    print("Recall@1 @10 @50 @100 = ",recallAt(TPR,P,FPR,N,1),recallAt(TPR,P,FPR,N,10),recallAt(TPR,P,FPR,N,50),recallAt(TPR,P,FPR,N,100))
    print("Hits@1 @10 @50 @100 = ",hitAt(TPR,P,FPR,N,1),hitAt(TPR,P,FPR,N,10),hitAt(TPR,P,FPR,N,50),hitAt(TPR,P,FPR,N,100))
    print("Total positives = ", P)
    

    TPR,FPR,ROC,P,N = computeROC_AUC(FGs_HDs_sim,HDs_keys,FGs_keys,F_mouse_p, F_human_p,"h")
    
    
    #with open('TPR_FPR_Human_'+test_tag+"_"+sys.argv[1]+'.json', 'w') as fb:
    #    json.dump([TPR,FPR],fb)
    print("Human Evaluation")
    CI = 2*(math. sqrt(ROC*(1-ROC)/(min([P,N]))))
    print("ROC_AUC = ",ROC)
    print("CI = ",CI)
    print(str(float("{:.3f}".format(ROC)))+"\mp"+str(float("{:.3f}".format(CI))) )
    print("Precision@1 @10 @50 @100 = ",precisionAt(TPR,P,FPR,N,1),precisionAt(TPR,P,FPR,N,10),precisionAt(TPR,P,FPR,N,50),precisionAt(TPR,P,FPR,N,100))
    print("Recall@1 @10 @50 @100 = ",recallAt(TPR,P,FPR,N,1),recallAt(TPR,P,FPR,N,10),recallAt(TPR,P,FPR,N,50),recallAt(TPR,P,FPR,N,100))
    print("Hits@1 @10 @50 @100 = ",hitAt(TPR,P,FPR,N,1),hitAt(TPR,P,FPR,N,10),hitAt(TPR,P,FPR,N,50),hitAt(TPR,P,FPR,N,100))
    print("Total positives = ", P)

    
