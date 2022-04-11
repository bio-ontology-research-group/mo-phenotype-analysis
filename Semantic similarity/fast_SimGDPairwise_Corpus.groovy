// -----------------------------------------------
// Sarah M. Alghamdi
// -----------------
// this implementation of semantic similarity allows for the use of equivilant class axioms,
// therefore realize the direct relations between classes and deal with them as one class, this is done on the annotaion level...
// arguuments are :
// output_files_prefix  input-gene-phenotype_association_file  input_disease-phenotype_association_file ontology_path
// output :
// 
// -----------------------------------------------

@Grab(group='com.github.sharispe', module='slib-sml', version='0.9.1')
@Grab(group='org.codehaus.gpars', module='gpars', version='1.1.0')
@Grab(group='net.sourceforge.owlapi', module='owlapi-api', version='4.2.5')
@Grab(group='net.sourceforge.owlapi', module='owlapi-apibinding', version='4.2.5')
@Grab(group='net.sourceforge.owlapi', module='owlapi-impl', version='4.2.5')



@Grab(group='org.slf4j', module='slf4j-simple', version='1.6.1')
@Grab(group='net.sourceforge.owlapi', module='owlapi-parsers', version='4.2.5')


import org.semanticweb.owlapi.model.parameters.*

import org.semanticweb.elk.reasoner.config.*
import org.semanticweb.owlapi.reasoner.*
import org.semanticweb.owlapi.reasoner.structural.StructuralReasoner
import org.semanticweb.owlapi.vocab.OWLRDFVocabulary;
import org.semanticweb.owlapi.model.*;
import org.semanticweb.owlapi.io.*;
import org.semanticweb.owlapi.owllink.*;
import org.semanticweb.owlapi.search.*;
import org.semanticweb.owlapi.manchestersyntax.renderer.*;
import org.semanticweb.owlapi.reasoner.structural.*
import java.io.File;
import org.semanticweb.owlapi.util.OWLOntologyMerger;





import java.net.*
import org.openrdf.model.vocabulary.*
import org.openrdf.sail.memory.model.MemURI
import slib.sglib.io.loader.*
import slib.sml.sm.core.metrics.ic.utils.*
import slib.sml.sm.core.utils.*
import slib.sglib.io.loader.bio.obo.*
import org.openrdf.model.URI
import slib.graph.algo.extraction.rvf.instances.*
import slib.sglib.algo.graph.utils.*
import slib.utils.impl.Timer
import slib.graph.algo.extraction.utils.*
import slib.graph.model.graph.*
import slib.graph.model.repo.*
import slib.graph.model.impl.graph.memory.*
import slib.sml.sm.core.engine.*
import slib.graph.io.conf.*
import slib.graph.model.impl.graph.elements.*
import slib.graph.algo.extraction.rvf.instances.impl.*
import slib.graph.model.impl.repo.*
import slib.graph.io.util.*
import slib.graph.io.loader.*
import groovyx.gpars.GParsPool
import org.semanticweb.owlapi.apibinding.OWLManager;
import org.semanticweb.owlapi.model.*;

System.setProperty("jdk.xml.entityExpansionLimit", "0");
System.setProperty("jdk.xml.totalEntitySizeLimit", "0");


// Parameters 
grouping = "BMA" // choose between "BMM" or "BMA"
ic_calculation_based_on_genes_only = true


OWLOntologyManager manager = OWLManager.createOWLOntologyManager()
OWLDataFactory fac = manager.getOWLDataFactory()


ConsoleProgressMonitor progressMonitor = new ConsoleProgressMonitor()
OWLReasonerConfiguration config = new SimpleConfiguration(progressMonitor)

//load ontology
def ontology = manager.loadOntologyFromOntologyDocument(new File( args[3]))

def clean_classes = {cl->cl
    return cl.toString().replace(">","").replace("<","")
}

class2equiv = [:]


ontology.getClassesInSignature(true).each {
 cl ->
 cl_url = clean_classes(cl)
  EntitySearcher.getEquivalentClasses(cl ,ontology).each {
   cl2 ->

   
   cl2_url = clean_classes(cl)


    if(!class2equiv.containsKey(cl_url) && !class2equiv.containsKey(cl2_url)){
      class2equiv[cl_url] = cl_url
      class2equiv[cl2_url] = cl_url
    }else if(!class2equiv.containsKey(cl_url)){
      class2equiv[cl_url] = class2equiv[cl2_url]

    }else if(!class2equiv.containsKey(cl2_url)){
      class2equiv[cl2_url] = class2equiv[cl_url]

    }

  }

  if(!class2equiv.containsKey(cl_url)){
    class2equiv[cl_url] = cl_url
  }


}



def factory = URIFactoryMemory.getSingleton()
OWLThing = factory.getURI("http://www.w3.org/2002/07/owl#Thing")
class Gene {

  int id
  Set annotations
  MemURI owlThing = URIFactoryMemory.getSingleton().getURI("http://www.w3.org/2002/07/owl#Thing")

  public Gene(id_i, annotations_i) {
    id = id_i
    annotations = annotations_i
  }

  void addAnnotation(annotation) {
    annotations.add(annotation);
  }

  Set getAnnotations() {
    if(annotations.size()!=0){annotations}
    else{
      println(id)
      annotations.add(owlThing)
      annotations
    }
  }
}

//a = new Gene(1,[1,2,3])
//println(a.getAnnotations())

URI graph_uri = factory.getURI("http://purl.obolibrary.org/obo/")
G graph = new GraphMemory(graph_uri)

def getGeneOntology = {
  // Load OBO file to graph "go.obo"
  GDataConf goConf = new GDataConf(GFormat.RDF_XML, args[3])
  GraphLoaderGeneric.populate(goConf, graph)

  // Add virtual root for 3 subontologies__________________________________
  URI virtualRoot = factory.getURI("http://purl.obolibrary.org/obo/virtualRoot")
  graph.addV(virtualRoot)
  GAction rooting = new GAction(GActionType.REROOTING)
  rooting.addParameter("root_uri", virtualRoot.stringValue())
  GraphActionExecutor.applyAction(factory, rooting, graph)
  phenomNet_classes = graph.getV()


  return graph
}

// BMA+Resnik, BMA+Schlicker2006, BMA+Lin1998, BMA+Jiang+Conrath1997,
// DAG-GIC, DAG-NTO, DAG-UI

graph = getGeneOntology()


new File(args[1]).splitEachLine('\t') { items ->
            String id = items[0];
            URI idURI = factory.getURI("http://purl.obolibrary.org/obo/" + id);
           //id = URLEncoder.encode(id);
           for (int j=1; j < items.size(); j++){
             String pheno = items[j];
             URI phenoURI = factory.getURI(class2equiv[pheno])//;
             Edge e = new Edge(idURI, RDF.TYPE, phenoURI);
             //println e.toString()
             graph.addE(e);
           }
         }
if(!ic_calculation_based_on_genes_only){

    new File(args[2]).splitEachLine('\t') { items ->
            String id = items[0];
            //println(id)
            URI idURI = factory.getURI("http://purl.obolibrary.org/obo/" + id);
           //id = URLEncoder.encode(id);
           for (int j=1; j < items.size(); j++){
             String pheno = items[j];
             //println(pheno)
             URI phenoURI = factory.getURI(class2equiv[pheno])//;
             Edge e = new Edge(idURI, RDF.TYPE, phenoURI);
             //println e.toString()
             graph.addE(e);
           }
         }



}


def getURIfromGO = { go ->
  return factory.getURI("http://purl.obolibrary.org/obo/" + go)
}

def gene_names = []
def getGenes = {
  def genes = []
  def i = 0
  new File(args[1]).splitEachLine('\t') { items ->
    def s = 0
    gene_names.add(items[0])
    genes.push(new Gene(i, new LinkedHashSet()))
    for (int j = 1; j < items.size(); j++) {
      if(ontology.containsClassInSignature(IRI.create(class2equiv[items[j]]))) {
          genes[i].addAnnotation(factory.getURI(class2equiv[items[j]]))
      }
       
    }
    i++
  }
  return genes
}

def dis_names = []
def getDiseases = {
  def dis = []
  def i = 0
  new File(args[2]).splitEachLine('\t') { items ->
    def s = 0
    dis_names.add(items[0])
    dis.push(new Gene(i, new LinkedHashSet()))
    for (int j = 1; j < items.size(); j++) {
       //println(items[j])
      // dis[i].addAnnotation(factory.getURI("http://purl.obolibrary.org/obo/" + items[j]))
       if(ontology.containsClassInSignature(IRI.create(class2equiv[items[j]]))) {
          dis[i].addAnnotation(factory.getURI(class2equiv[items[j]]))
       }
    }
    i++
  }
  return dis
}

genes = getGenes()
diseases = getDiseases()


//GraphActionExecutor.applyAction(factory, rooting, graph);
ICconf icConf = new IC_Conf_Corpus("ResnikIC", SMConstants.FLAG_IC_ANNOT_RESNIK_1995_NORMALIZED);

//this.smConfPairwise = new SMconf("Jaccard", SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_JACCARD_IC  );
//this.smConfPairwise = new SMconf("Resnik", SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_LIN_1998  );
this.smConfPairwise = new SMconf("Resnik", SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_RESNIK_1995 );

if(grouping == "BMA"){
 this.smConfGroupwise = new SMconf("BMA", SMConstants.FLAG_SIM_GROUPWISE_BMA); 
}
if(grouping == "BMM"){
 this.smConfGroupwise = new SMconf("BMM", SMConstants.FLAG_SIM_GROUPWISE_BMM); 
}
smConfPairwise.setICconf(icConf);


def result = new Double[genes.size() * diseases.size()]
for (i = 0; i < result.size(); i++) {
  result[i] = 0
}


SM_Engine engine = new SM_Engine(graph)

def c = 0

GParsPool.withPool {
  result.eachParallel { val ->
    def i = val.toInteger()
    def x = i.intdiv(diseases.size())
    def y = i % diseases.size()
    def genlist = genes[x].getAnnotations()
    def dislist =  diseases[y].getAnnotations()
    result[i] = engine.compare(
            smConfGroupwise,
            smConfPairwise,
            genlist,
            dislist)
    if (c % 100000 == 0)
      println c
    c++
  }
}


def fout = new PrintWriter(new BufferedWriter(
  new FileWriter(args[0]+"sim_scores_resnik.txt")))
for (i = 0; i < result.size(); i++) {
  def x = i.intdiv(diseases.size()) //
  def y = i % diseases.size()
  fout.println(result[i])

}
fout.flush()
fout.close()

def disout = new PrintWriter(new BufferedWriter(
  new FileWriter(args[0]+'diseases_resnik.txt')))
for (i = 0; i < diseases.size(); i++){
  disout.println(dis_names[i])
}
disout.flush()
disout.close()

def genout = new PrintWriter(new BufferedWriter(
  new FileWriter(args[0]+'genes_resnik.txt')))
  for (i = 0; i < genes.size(); i++){
    genout.println(gene_names[i].toString())
  }
genout.flush()
genout.close()
