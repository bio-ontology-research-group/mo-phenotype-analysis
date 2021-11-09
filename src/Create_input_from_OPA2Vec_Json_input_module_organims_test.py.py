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


#use this
init = 'test_name'

#with open('/home/alghsm0a/alghsm0a/Desktop/OPA2Vec/opa2vec-master/opa2vec-experimental/mouse_'+org+'_.json','r') as f:
with open('vectors.json','r') as f:
    all_vectors = json.load(f)

MGs = { key:value for key,value in all_vectors.items() if "MGI" in key and "http" not in key}
HDs = { key:value for key,value in all_vectors.items() if (("OMIM" in key or "ORPHA" in key or "DECIPHER" in key) and "http" not in key)}
YGs = { key:value for key,value in all_vectors.items() if "SP" in key and "http" not in key}
FBGs = { key:value for key,value in all_vectors.items() if "FBgn" in key and "http" not in key}
ZFGs = { key:value for key,value in all_vectors.items() if "ZDB-GENE" in key and "http" not in key}


print(len(MGs.keys()))
print(len(HDs.keys()))
print(len(YGs.keys()))
print(len(FBGs.keys()))
print(len(ZFGs.keys()))


maps = [MGs, YGs, FBGs,ZFGs]
#maps = [MGs, HDs, FBGs]
i = 0
names = ["MGs.json", "YGs.json", "FBGs.json", "ZFGs.json" ]#"ZFGs.json"]# "YGs.json", "FBGs.json"]#, 
#names = ["MGs.json", "HDs.json", "FBGs.json"]


for map in maps:
    with open(init+names[i], 'w') as fp:
        json.dump(map, fp)
    i+=1


HDs_vectors=[]
HDs_keys = list(HDs.keys())

for key in HDs_keys:
    HDs_vectors.append(HDs[key])
with open(init+'_opa_HDs_keys.json', 'w') as fp:
    json.dump(HDs_keys, fp)

for fil in names:
    with open(init+fil,'r') as f:
        OGs = json.load(f)

    OGs_keys = list(OGs.keys())
    with open(init+'_opa_'+fil.replace('.json','')+'_keys.json', 'w') as fp:
        json.dump(OGs_keys, fp)
    OGs_vectors=[]
    for key in OGs_keys:
        OGs_vectors.append(OGs[key])

    OGs_HDs= cosine_similarity(np.array(OGs_vectors),np.array(HDs_vectors))

    print(fil, 'OK')
    print(OGs_HDs[0][0])
    print(OGs_keys[0],HDs_keys[0])
    filename = init+'_opa_'+fil.replace('.json','')+'_sim.npy'
    np.save(filename,OGs_HDs)
