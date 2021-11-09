
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



//ont =manager.loadOntologyFromOntologyDocument(IRI.create("http://purl.obolibrary.org/obo/hp.owl"))
ont =manager.loadOntologyFromOntologyDocument(new File("/home/alghsm0a/Ontologies/Pheno-e.owl"))
init = "Pheno_e"

ConsoleProgressMonitor progressMonitor = new ConsoleProgressMonitor()
OWLReasonerConfiguration config = new SimpleConfiguration(progressMonitor)
ElkReasonerFactory f = new ElkReasonerFactory()
OWLReasoner reasoner = f.createReasoner(ont,config)

all_parents =[:]
direct_parents = [:]
direct_children = [:]
all_chlidren = [:]

def classToId(cl){
	cl.toString().replaceAll('<http://purl.obolibrary.org/obo/','').replaceAll('<http://aber-owl.net/phenotype.owl#','').replaceAll('>','').replaceAll('_',':')
}


ont.getClassesInSignature(true).each {cl ->
	clID = classToId(cl)
	direct_children[clID] = []
	all_chlidren[clID] = []
	all_parents[clID] = []
	direct_parents[clID] = []
/*
	reasoner.getSubClasses(cl, true).each { n-> n.each{ child->
	    direct_children[clID].add(classToId(child))
	} }
	reasoner.getSubClasses(cl, false).each { n-> n.each{ child->
	    all_chlidren[clID].add(classToId(child))
	} }
*/
	reasoner.getSuperClasses(cl, false).each { n-> n.each{ parent->
	    all_parents[clID].add(classToId(parent))
	} }
/*	
	reasoner.getSuperClasses(cl, true).each { n-> n.each{ parent->
	    direct_parents[clID].add(classToId(parent))
	} }
*/	
}


/*

def output = JsonOutput.toJson(direct_children)
new File(init+'_direct_children.json').withWriter('utf-8') { 
         writer -> writer.writeLine output
}

output = JsonOutput.toJson(all_chlidren)
new File(init+'_all_chlidren.json').withWriter('utf-8') { 
         writer -> writer.writeLine output
}
*/
output = JsonOutput.toJson(all_parents)
new File(init+'_all_parents.json').withWriter('utf-8') { 
         writer -> writer.writeLine output
}
/*
output = JsonOutput.toJson(direct_parents)
new File(init+'_direct_parents.json').withWriter('utf-8') { 
         writer -> writer.writeLine output
}
*/