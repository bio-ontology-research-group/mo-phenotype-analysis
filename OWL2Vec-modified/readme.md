This code is a modification to OWL2Vec* downloaded on 31 May 2021
original link is here

The addition we provide here is to add walks that starts from additional entities(like genes, diseases, protiens, ...) that are annotatated by some ontology classes
the additional annotaions will be in the following format (entity annotation_url) separated by space:

MGI:2673128 http://purl.obolibrary.org/obo/MP_0001262
MGI:2673128 http://purl.obolibrary.org/obo/MP_0001732
MGI:1913870 http://purl.obolibrary.org/obo/MP_0020069
MGI:1913870 http://purl.obolibrary.org/obo/MP_0003989
MGI:1913870 http://purl.obolibrary.org/obo/MP_0008948



you will then need to add to the defalt file this line:
annotation_extra=path_to_file
