package data_parsing;

import com.github.phenomics.ontolib.formats.hpo.HpoDiseaseAnnotation;
import com.github.phenomics.ontolib.io.base.TermAnnotationParserException;
import com.github.phenomics.ontolib.io.obo.hpo.HpoDiseaseAnnotationParser;
import ontologizer.association.Association;
import ontologizer.association.AssociationContainer;
import ontologizer.association.Gene2Associations;
import ontologizer.go.*;
import ontologizer.types.ByteString;
import sonumina.math.graph.DirectedGraph;
import sonumina.math.graph.SlimDirectedGraphView;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.Map;

import static test.NewRefinedBOQATest.getOboParser;

//Phenotype annotation tab parser
public class PATParser {
    HashMap<Term, HashMap<ByteString, Integer>> pheno_disease_freq;
    public void initializeHashmap(TermContainer tc)
    {
        pheno_disease_freq = new HashMap<>();
        for (Term t : tc)
        {
            pheno_disease_freq.put(t, new HashMap<ByteString, Integer>());
        }
    }

    public Map<Term, Integer> getHPOToFreqMapping(TermContainer tc)
    {
        //check if the term is already in there!
//        tc.get().
        return null;

    }

    public PATParser(){

    }

    public void doParse()  throws OBOParserException, IOException, URISyntaxException
    {

        SlimDirectedGraphView<Term> slim;
        OBOParser hpoParser = getOboParser();
        AssociationContainer assocs;
        TermContainer tc = new TermContainer(hpoParser.getTermMap(), hpoParser.getFormatVersion(), hpoParser.getDate());
        Ontology ontology = new Ontology(tc);
        Map<Term, Integer> hpo2freq = getHPOToFreqMapping(tc);
        initializeHashmap(tc);
        //need to add all the diseases in... the issue is that: we don't have a group by!
        //hence, we should check the diseaseContainer if it is already there


        assocs = new AssociationContainer();

        System.out.println("Working Directory = " +
                System.getProperty("user.dir"));
        File inputFile = new File("C:\\Users\\johnp\\Desktop\\git_stuff\\boqa\\new_boqa\\resources\\phenotype_annotation.tab");
        System.out.println(inputFile.toString());
        try {
            HpoDiseaseAnnotationParser parser = new HpoDiseaseAnnotationParser(inputFile);
            while (parser.hasNext()) {
                HpoDiseaseAnnotation anno = parser.next();

                if (anno.getDb().equals(OMIM_name))
                {
                    System.out.println(anno);
                }

                ByteString item = new ByteString(anno.getDbName());

//                Term t = slim.getVertex()
//                Term t  = new Term(new TermID(anno.getTermId().getIdWithPrefix()));

                TermID t = new TermID(anno.getTermId().getIdWithPrefix());

                //This is the solution!
                Term tx = tc.get(t);


                Association a = new Association(item, tx.getIDAsString());

                //it will auto group it if its already in there!
//                if (assocs.containsGene(item))
//                {
//                    assocs.
                //here we are simply required to remember what TID and temr was.
                //we want to be able to update the indices, based on the vertex2ancestor info from before, and
                //we CAN do that!
                //since for example, the interface between
                //let us make it a mapping between terms, and items and frequencies
                pheno_disease_freq.get(t).put(item, hpo2freq.get(tc.get(anno.getFrequencyModifier())));



                assocs.addAssociation(a);

            }
        } catch (IOException e) {
            System.err.println("Problem reading from file.");
        } catch (TermAnnotationParserException e) {
            System.err.println("Problem parsing file.");
        }
    }

    static final String OMIM_name = "OMIM"; //Constant
    public static void main(String[]args) throws OBOParserException, IOException, URISyntaxException {
        PATParser p = new PATParser();
        p.doParse();

    }

}
