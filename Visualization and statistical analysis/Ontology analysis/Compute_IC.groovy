
// -----------------------------------------------
// Sarah M. Alghamdi
// -----------------
// arguuments are :
// ontology_path  input-gene-phenotype_association_file

// -----------------------------------------------

@Grab(group='com.github.sharispe', module='slib-sml', version='0.9.1')
@Grab(group='org.codehaus.gpars', module='gpars', version='1.1.0')
@Grab(group='net.sourceforge.owlapi', module='owlapi-api', version='4.2.5')
@Grab(group='net.sourceforge.owlapi', module='owlapi-apibinding', version='4.2.5')
@Grab(group='net.sourceforge.owlapi', module='owlapi-impl', version='4.2.5')

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
import slib.graph.algo.accessor.GraphAccessor;

System.setProperty("jdk.xml.entityExpansionLimit", "0");
System.setProperty("jdk.xml.totalEntitySizeLimit", "0");


// Parameters 
grouping = "BMA" // choose between "BMM" or "BMA"
ic_calculation_based_on_genes_only = true





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
  GDataConf goConf = new GDataConf(GFormat.RDF_XML, args[0])
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
//println(args)

// BMA+Resnik, BMA+Schlicker2006, BMA+Lin1998, BMA+Jiang+Conrath1997,
// DAG-GIC, DAG-NTO, DAG-UI
OWLOntologyManager manager = OWLManager.createOWLOntologyManager()
ontology = manager.loadOntologyFromOntologyDocument(IRI.create(new File ( args[0]) ))
graph = getGeneOntology()


new File(args[1]).splitEachLine('\t') { items ->
            String id = items[0];
            URI idURI = factory.getURI("http://purl.obolibrary.org/obo/" + id);
           //id = URLEncoder.encode(id);
           for (int j=1; j < items.size(); j++){
             String pheno = items[j];
             URI phenoURI = factory.getURI(pheno);
             Edge e = new Edge(idURI, RDF.TYPE, phenoURI);
             //println e.toString()
             graph.addE(e);
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
      if(ontology.containsClassInSignature(IRI.create(items[j]))) {
          genes[i].addAnnotation(factory.getURI(items[j]))
      }
       
    }
    i++
  }
  return genes
}


genes = getGenes()


//GraphActionExecutor.applyAction(factory, rooting, graph);
ICconf icConf = new IC_Conf_Corpus("ResnikIC", SMConstants.FLAG_IC_ANNOT_RESNIK_1995);//_NORMALIZED );

SM_Engine engine = new SM_Engine(graph)


Set<URI> goTerms = GraphAccessor.getClasses(graph);
        System.out.println("GO terms : " + goTerms.size());

 goTerms.each{goTerm->
  if(goTerm.toString() != "http://www.w3.org/2002/07/owl#Nothing"){
    println(goTerm.toString() + "\t" + engine.getIC(icConf,goTerm).toString());
  }

    }