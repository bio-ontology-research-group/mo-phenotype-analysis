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


init = 'test_name'
org = 'yeast'
#OPA human
#get human gene human disease vectors
with open('vecors.json','r') as f:
    all_vectors = json.load(f)

HGs = { key:value for key,value in all_vectors.items() if key.isdigit()}
HDs = { key:value for key,value in all_vectors.items() if (("OMIM" in key or "ORPHA" in key or "DECIPHER" in key) and "http" not in key)}

with open('human_mouse_opa_hg.json', 'w') as fp:
        json.dump(HGs, fp)
with open('human_mouse_opa_hd.json', 'w') as fp:
        json.dump(HDs, fp)

# compute the gene disease similarity matrices
HDs_vectors=[]
HDs_keys = list(HDs.keys())
print(len(HDs_keys))
for key in HDs_keys:
    HDs_vectors.append(HDs[key])

with open('human_mouse_opa_hd_keys.json', 'w') as fp:
    json.dump(HDs_keys, fp)

HGs_vectors=[]
HGs_keys = list(HGs.keys())

for key in HGs_keys:
    HGs_vectors.append(HGs[key])
print(len(HGs_keys))
with open('human_mouse_opa_hg_keys.json', 'w') as fp:
    json.dump(HGs_keys, fp)

print(len(HGs_vectors), len(HDs_vectors))
OGs_HDs= cosine_similarity(np.array(HGs_vectors),np.array(HDs_vectors))
np.save('human_mouse_opa_sim.npy',OGs_HDs)
