//-----------------------------------------------------------
// Sarah M. Algamdi
//-----------------------------------------------------------
// This code makes the matrex calculation of ralating phenotypes super and subclasses and shows inter and intra connections
//----------------------------------------------------------- 


@Grapes([
 @Grab(group='org.slf4j', module='slf4j-simple', version='1.6.1'),
 @Grab(group = 'org.semanticweb.elk', module = 'elk-owlapi', version = '0.4.2'),
 @Grab(group = 'net.sourceforge.owlapi', module = 'owlapi-api', version = '4.2.5'),
 @Grab(group = 'net.sourceforge.owlapi', module = 'owlapi-apibinding', version = '4.2.5'),
 @Grab(group = 'net.sourceforge.owlapi', module = 'owlapi-impl', version = '4.2.5'),
 @Grab(group = 'net.sourceforge.owlapi', module = 'owlapi-parsers', version = '4.2.5'),
 @GrabConfig(systemClassLoader = true)
])


import org.semanticweb.owlapi.model.parameters.*
import org.semanticweb.elk.owlapi.ElkReasonerFactory;
import org.semanticweb.elk.owlapi.ElkReasonerConfiguration
import org.semanticweb.elk.reasoner.config.*
import org.semanticweb.owlapi.apibinding.OWLManager;
import org.semanticweb.owlapi.reasoner.*
import org.semanticweb.owlapi.reasoner.structural.StructuralReasoner
import org.semanticweb.owlapi.vocab.OWLRDFVocabulary;
import org.semanticweb.owlapi.model.*;
import org.semanticweb.owlapi.io.*;
import org.semanticweb.owlapi.owllink.*;
import org.semanticweb.owlapi.util.*;
import org.semanticweb.owlapi.search.*;
import org.semanticweb.owlapi.manchestersyntax.renderer.*;
import org.semanticweb.owlapi.reasoner.structural.*
import groovy.json.JsonOutput
import java.io.File 


OWLOntologyManager manager = OWLManager.createOWLOntologyManager()
OWLDataFactory fac = manager.getOWLDataFactory()
StructuralReasonerFactory f1 = new StructuralReasonerFactory()



ont =manager.loadOntologyFromOntologyDocument(new File(args[0]))


ConsoleProgressMonitor progressMonitor = new ConsoleProgressMonitor()
OWLReasonerConfiguration config = new SimpleConfiguration(progressMonitor)
ElkReasonerFactory f = new ElkReasonerFactory()
OWLReasoner reasoner = f.createReasoner(ont,config)

direct= false

def classToId(cl){
	cl.toString().replaceAll('<','').replaceAll('>','')

}

/*
class_type = ["HP","MP", "PHENO", "FB","FYPO"]
organism = ["Human","Mouse","Fish","Fly", "Yeast"] 
//subclass_matrix = [[0]*5]*5
subclass_matrix = [[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]
superclass_matrix = [[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]

class_type.eachWithIndex{type, index ->

	type_count = 0
	ont.getClassesInSignature(true).each {cl ->

		clID = classToId(cl)

		//print(clID)

		if (clID.indexOf(type)>-1 ){
			type_count++



			ls=reasoner.getSubClasses(cl, direct).getFlattened()
			
		   	flag = [false]*5
		   	class_type.eachWithIndex{type2, index2 ->
			   ls.each{ child ->
			   	child_id = classToId(child)
			   		if (child_id.indexOf(type2)>-1 ){
			   			flag[index2]=true
				   		}
					}
				

				}
			flag.eachWithIndex{o,i ->
				if(o){
					subclass_matrix[index][i]+=1

					}
				}

				
			
			ls=reasoner.getSuperClasses(cl, direct).getFlattened()
			
		   	flag = [false]*5
		   	class_type.eachWithIndex{type2, index2 ->
			   ls.each{ parent ->
			   	parent_id = classToId(parent)
			   		if (parent_id.indexOf(type2)>-1 ){
			   			flag[index2]=true
				   		}
					}
				

				}
			flag.eachWithIndex{o,i ->
				if(o){
					superclass_matrix[index][i]+=1

					}
				}
		}



	}
	println(type)
	println(type_count)

}

*/

class_type = ["HP","MP", "PHENO", "FB","FYPO"]
organism = ["Human","Mouse","Fish","Fly", "Yeast"] 
//subclass_matrix = [[0]*5]*5
subclass_matrix = [[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]
superclass_matrix = [[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]




	type_count = 0
	ont.getClassesInSignature(true).each {cl ->

		clID = classToId(cl)


		//if (clID.indexOf("UPHENO_")>-1){



			ls=reasoner.getSubClasses(cl, direct).getFlattened()
			
		   	class_type.eachWithIndex{type, index ->
		   		if(clID.indexOf(type)==-1){
				   	class_type.eachWithIndex{type2, index2 ->
				   		if(clID.indexOf(type2)>-1){
					   	   flag1 = false
					   	   flag2 = false
						   ls.each{ child ->
						   	child_id = classToId(child)
						   	if (child_id.indexOf(type)==-1 ){
						   			flag1 = true
							   		}
								
						    if (child_id.indexOf(type2)>-1 ){
						   			flag2 = true
							   		}
								}

							if(flag1 && flag2){
								subclass_matrix[index][index2]++

							}
						}
					}
				}
			}
			

				




	//}
}
	










println("subclass_matrix")
println(class_type.join("\t"))
subclass_matrix.eachWithIndex{line,i->
    print(class_type[i]+"\t")
	println(line.join("\t"))
}



println("superclass_matrix")
println(class_type.join("\t"))
superclass_matrix.eachWithIndex{line,i->
    print(class_type[i]+"\t")
	println(line.join("\t"))
}
