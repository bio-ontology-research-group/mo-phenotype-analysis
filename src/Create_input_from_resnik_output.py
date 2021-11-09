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

# for resnik's similarity

sim_scores = np.loadtxt(sys.argv[1], dtype = 'float32')
diseases = np.genfromtxt(sys.argv[2], dtype = 'str')
genes = np.genfromtxt(sys.argv[3], dtype = 'str')

sim_mat = sim_scores.reshape((len(genes),len(diseases)))

np.save(sys.argv[4]+"_resnik_sim",sim_mat)
with open(sys.argv[4]+'_resnik_g_keys.json', 'w') as fp:
    json.dump(genes.tolist(), fp)
with open(sys.argv[4]+'_resnik_d_keys.json', 'w') as fp:
    json.dump(diseases.tolist(), fp)

