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


embedding_model = sys.argv[2]
word2vec_model=gensim.models.Word2Vec.load(embedding_model)
word_vocabs = word2vec_model.wv.vocab

HGs={}
HDs={}

for key in word_vocabs:
    if(key.isdigit()):
        HGs[key] = word2vec_model[key].tolist()

    elif(("OMIM" in key) and ("http" not in key)):
        HDs[key] = word2vec_model[key].tolist()


with open(init+'HG.json', 'w') as fp:
        json.dump(HGs, fp)
with open(init+'HD.json', 'w') as fp:
        json.dump(HDs, fp)
# compute the gene disease similarity matrices
HDs_vectors=[]
HDs_keys = list(HDs.keys())
print(len(HDs_keys))
for key in HDs_keys:
    HDs_vectors.append(HDs[key])        

with open(init+'_d_keys.json', 'w') as fp:
    json.dump(HDs_keys, fp)

HGs_vectors=[]
HGs_keys = list(HGs.keys())

for key in HGs_keys:
    HGs_vectors.append(HGs[key])
print(len(HGs_keys))
with open(init+'_g_keys.json', 'w') as fp:
    json.dump(HGs_keys, fp)

print(len(HGs_vectors), len(HDs_vectors))
OGs_HDs= cosine_similarity(np.array(HGs_vectors),np.array(HDs_vectors))
np.save(init+'_dl_sim.npy',OGs_HDs)
