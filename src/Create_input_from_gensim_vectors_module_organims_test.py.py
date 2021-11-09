# Sarah M. Alghamdi
#------------------------------------
# args
# sim_file disease_file genes_file outputs_prefix
#------------------------------------
import json
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
#------------------------------------
import json
import gensim
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
#------------------------------------
# separate the vectors

#init = str(sys.argv[1])

#from WV model:
#load mosule:



MGs, HDs, YGs, FBGs, ZFGs= {},{},{},{},{}
maps = [MGs, YGs, FBGs, ZFGs]
names = ["MGs.json", "YGs.json","FBGs.json","ZFGs.json"]
organism = ["mouse","yeast","fly","fish"]
annot=["MGI","SP","FBgn","ZDB-GENE"]

embedding_model = "inputs/to_compare_pheno_e"
word2vec_model=gensim.models.Word2Vec.load(embedding_model)
word_vocabs = word2vec_model.wv.vocab

for org in range(len(organism)):
    for key in word_vocabs:
        if(annot[org] in key and "http" not in key):
            maps[org][key] = word2vec_model[key].tolist()
HDs={}            
for key in word_vocabs:
    if(("OMIM" in key or "ORPHA" in key or "DECIPHER" in key) and ("http" not in key)):
        HDs[key] = word2vec_model[key].tolist()
       
i=0
for map_ in maps:
    with open(init+names[i], 'w') as fp:
        json.dump(map_, fp)
    i+=1

with open(init+"HDs.json", 'w') as fp:
    json.dump(HDs, fp)

HDs_keys = list(HDs.keys())

# compute the gene disease similarity matrices


HDs_vectors=[]
for key in HDs_keys:
    HDs_vectors.append(HDs[key])
with open(init+"HDs_keys.json", 'w') as fp:
    json.dump(HDs_keys, fp)
for og in range(len(maps)):   
    OGs_keys = list(maps[og].keys())
    with open(init+names[og].replace('.json','')+'_keys.json', 'w') as fp:
        json.dump(OGs_keys, fp)

    OGs_vectors=[]
    for key in OGs_keys:
        OGs_vectors.append(maps[og][key])

    OGs_HDs= cosine_similarity(np.array(OGs_vectors),np.array(HDs_vectors))

    print(names[og], 'OK')
    print(OGs_HDs[0][0])
    print(OGs_keys[0],HDs_keys[0])
    filename = names[og].replace('.json','')+'_HDs_sim'
    np.save(init+filename,OGs_HDs)
