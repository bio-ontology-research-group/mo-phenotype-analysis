
This repository contains all the scripts needed to test model organism phenotypes contribution in predicting gene disease associations using ontology-based phenotypic similarity


## Usage
To generate the disease gene phenotypic similarity we use two main approaches (Traditional semantic similarity, and Ontology embedding similarity) :

**1- Traditional semantic similarity**

Here we apply Resnik's semantic similarity using best maching approach of groupwise similarity

To generate this similarity on classified ontology run the code [Semantic similarity/fast_SimGDPairwise_Corpus.groovy](https://github.com/bio-ontology-research-group/mo-phenotype-analysis/blob/main/Semantic%20similarity/fast_SimGDPairwise_Corpus.groovy)

The input format of the phenotype annotations is a tab separated values, examples follows:

input-gene-phenotype_association_file

Example 1 for mouse gene-phenotypes:
```sh
MGI:104649    http://purl.obolibrary.org/obo/MP_0001544    http://purl.obolibrary.org/obo/MP_0004011
MGI:1096366    http://purl.obolibrary.org/obo/MP_0005092    http://purl.obolibrary.org/obo/MP_0004045    http://purl.obolibrary.org/obo/MP_0011704    http://purl.obolibrary.org/obo/MP_0005399    http://purl.obolibrary.org/obo/MP_0001732    http://purl.obolibrary.org/obo/MP_0000352
MGI:109352    http://purl.obolibrary.org/obo/MP_0000702    http://purl.obolibrary.org/obo/MP_0000691    http://purl.obolibrary.org/obo/MP_0004816    http://purl.obolibrary.org/obo/MP_0009796    http://purl.obolibrary.org/obo/MP_0002083    http://purl.obolibrary.org/obo/MP_0002023    http://purl.obolibrary.org/obo/MP_0003076    http://purl.obolibrary.org/obo/MP_0000693    http://purl.obolibrary.org/obo/MP_0000688    http://purl.obolibrary.org/obo/MP_0008412    http://purl.obolibrary.org/obo/MP_0002494    http://purl.obolibrary.org/obo/MP_0008498    http://purl.obolibrary.org/obo/MP_0004815    http://purl.obolibrary.org/obo/MP_0008943    http://purl.obolibrary.org/obo/MP_0000709
```

Example 2 for zebrafish gene-phenotypes:
```sh
ZDB-GENE-991019-6    http://aber-owl.net/phenotype.owl#PHENO:8    http://aber-owl.net/phenotype.owl#PHENO:6    http://aber-owl.net/phenotype.owl#PHENO:10    http://aber-owl.net/phenotype.owl#PHENO:1    http://aber-owl.net/phenotype.owl#PHENO:4    http://aber-owl.net/phenotype.owl#PHENO:3
ZDB-GENE-030131-9790    http://aber-owl.net/phenotype.owl#PHENO:21    http://aber-owl.net/phenotype.owl#PHENO:19    http://aber-owl.net/phenotype.owl#PHENO:12    http://aber-owl.net/phenotype.owl#PHENO:15    http://aber-owl.net/phenotype.owl#PHENO:17    http://aber-owl.net/phenotype.owl#PHENO:23
ZDB-GENE-040525-2    http://aber-owl.net/phenotype.owl#PHENO:30    http://aber-owl.net/phenotype.owl#PHENO:31    http://aber-owl.net/phenotype.owl#PHENO:27    http://aber-owl.net/phenotype.owl#PHENO:25
```
Example 3 for fly gene-phenotypes:
```sh
FBgn0020305    http://purl.obolibrary.org/obo/FBbtAB_00004893    http://purl.obolibrary.org/obo/FBcv_0000353    http://purl.obolibrary.org/obo/FBcv_0002015    http://purl.obolibrary.org/obo/FBcv_0002023
FBgn0000420    http://purl.obolibrary.org/obo/FBbtAB_00004729    http://purl.obolibrary.org/obo/FBcv_0000354
FBgn0040230    http://purl.obolibrary.org/obo/FBbtAB_00005116    http://purl.obolibrary.org/obo/FBbtAB_00004729    http://purl.obolibrary.org/obo/FBbtAB_00005837
FBgn0067779    http://purl.obolibrary.org/obo/FBbtAB_00004762    http://purl.obolibrary.org/obo/FBbtAB_00004765    http://purl.obolibrary.org/obo/FBbtAB_00005179    http://purl.obolibrary.org/obo/FBbtAB_00000046    http://purl.obolibrary.org/obo/FBbtAB_00000043    http://purl.obolibrary.org/obo/FBbtAB_00005169    http://purl.obolibrary.org/obo/FBbtAB_00004761
```

Example 4 for yeast gene-phenotype:
```sh
SPAC56F8.02    http://purl.obolibrary.org/obo/FYPO_0006680    http://purl.obolibrary.org/obo/FYPO_0000636    http://purl.obolibrary.org/obo/FYPO_0006930    http://purl.obolibrary.org/obo/FYPO_0000684    http://purl.obolibrary.org/obo/FYPO_0000088    http://purl.obolibrary.org/obo/FYPO_0006015    http://purl.obolibrary.org/obo/FYPO_0000085    http://purl.obolibrary.org/obo/FYPO_0003358    http://purl.obolibrary.org/obo/FYPO_0002061    http://purl.obolibrary.org/obo/FYPO_0000121
SPAC56F8.10    http://purl.obolibrary.org/obo/FYPO_0002698    http://purl.obolibrary.org/obo/FYPO_0003902    http://purl.obolibrary.org/obo/FYPO_0002697    http://purl.obolibrary.org/obo/FYPO_0000311    http://purl.obolibrary.org/obo/FYPO_0002061    http://purl.obolibrary.org/obo/FYPO_0000040
SPAC22A12.03c    http://purl.obolibrary.org/obo/FYPO_0000121
```


input_disease-phenotype_association_file

Example for human disease-phenotype:
```sh
OMIM:611771    http://purl.obolibrary.org/obo/HP_0000083    http://purl.obolibrary.org/obo/HP_0100820    http://purl.obolibrary.org/obo/HP_0012574    http://purl.obolibrary.org/obo/HP_0000093
OMIM:611773    http://purl.obolibrary.org/obo/HP_0030880    http://purl.obolibrary.org/obo/HP_0000107    http://purl.obolibrary.org/obo/HP_0005115    http://purl.obolibrary.org/obo/HP_0004944    http://purl.obolibrary.org/obo/HP_0001136    http://purl.obolibrary.org/obo/HP_0000006    http://purl.obolibrary.org/obo/HP_0000790    http://purl.obolibrary.org/obo/HP_0003394    http://purl.obolibrary.org/obo/HP_0000573    http://purl.obolibrary.org/obo/HP_0000083    http://purl.obolibrary.org/obo/HP_0001297    http://purl.obolibrary.org/obo/HP_0000112    http://purl.obolibrary.org/obo/HP_0002352
OMIM:611777    http://purl.obolibrary.org/obo/HP_0011705    http://purl.obolibrary.org/obo/HP_0001279    http://purl.obolibrary.org/obo/HP_0011712    http://purl.obolibrary.org/obo/HP_0001645    http://purl.obolibrary.org/obo/HP_0000006    http://purl.obolibrary.org/obo/HP_0012248    http://purl.obolibrary.org/obo/HP_0001663
```

After we calculate the semantic similarity we use [this script](https://github.com/bio-ontology-research-group/mo-phenotype-analysis/blob/main/src/Create_input_from_resnik_output.py) to generate the input of the gene-disease evaluation code


**2- Ontology embedding similarity**
Here we apply three different approaches to generate ontology embeddings, we used [OPA2Vec](https://github.com/bio-ontology-research-group/opa2vec), [OWL2Vec*](https://github.com/KRR-Oxford/OWL2Vec-Star) and [Dl2Vec](https://github.com/bio-ontology-research-group/DL2Vec)

We provide a modification of OWL2Vec* which allows the annotaions to be represented in the walks of the tool in Module_organism_phenotypes/OWL2Vec-modified

after genrating the vector representations of the entities annotated by the ontology using any of these tools we calcutate cosine similarities using the scripts 

[for vectors in json format](https://github.com/bio-ontology-research-group/mo-phenotype-analysis/blob/main/src/Create_input_from_OPA2Vec_Json_input_module_organims_test.py.py)

[for vectors in word2vec gensim format](https://github.com/bio-ontology-research-group/mo-phenotype-analysis/blob/main/src/Create_input_from_gensim_vectors_module_organims_test.py.py)


we also tested with supervised learning from the vectors representation code is available in [Simple_NN](https://github.com/bio-ontology-research-group/mo-phenotype-analysis/tree/main/Simple_NN)



## Phenotypes annotaions resources:

**Human disease annotations:**

we aquire this data from http://purl.obolibrary.org/obo/hp/hpoa/phenotype_annotation.tab downladed in Feb,2021

**Mouse phenotypes:**

we aquire this data from http://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt downladed in Feb,2021

**Zebrafish phenotypes:**

we aquire this data from https://zfin.org/downloads/phenoGeneCleanData_fish.txt downladed in Feb,2021

**Fruitfly phenotypes:**

we aquire this data from http://ftp.flybase.org/releases/FB2021_02/precomputed_files/alleles/allele_phenotypic_data_fb_2021_01.tsv.gz downladed in Feb,2021

**Yeast phenotypes:**

we aquire this data from https://www.pombase.org/downloads/phenotype-annotations downladed in Feb,2021




## Ontologies

Pheno-e can be downloaded from [here](http://aber-owl.net/ontology/Pheno-e/#/)

uPheno can be downloaded from [here](https://data.monarchinitiative.org/upheno2/current/upheno-release/all/upheno_all_with_relations.owl)


