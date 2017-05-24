package sonumina.boqa.tests;

/**
 * Created by johnchen on 18/05/17.
 */
import java.net.URI;
import ontologizer.association.Association;
import ontologizer.association.AssociationContainer;
import ontologizer.benchmark.Datafiles;
import ontologizer.go.*;
import ontologizer.types.ByteString;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import sonumina.boqa.benchmark.Benchmark;
import sonumina.boqa.calculation.BOQA;
import sonumina.boqa.calculation.ReducedBoqa;
import sonumina.boqa.calculation.Observations;
import sonumina.math.graph.AbstractGraph;
import sonumina.math.graph.SlimDirectedGraphView;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.*;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.HashSet;
import java.util.Random;
import java.util.logging.Level;

import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ontologizer.association.Association;
import ontologizer.association.AssociationContainer;
import ontologizer.benchmark.Datafiles;
import ontologizer.go.OBOParser;
import ontologizer.go.OBOParserException;
import ontologizer.go.Ontology;
import ontologizer.go.Term;
import ontologizer.go.TermContainer;
import ontologizer.types.ByteString;
import sonumina.boqa.benchmark.Benchmark;
import sonumina.boqa.calculation.BOQA;
import sonumina.boqa.calculation.BOQA.Result;
import sonumina.boqa.calculation.Observations;
import sonumina.math.graph.AbstractGraph.DotAttributesProvider;
import sonumina.math.graph.SlimDirectedGraphView;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class ReducedBOQATest
{
    private static Datafiles hpo;

    private Logger logger = LoggerFactory.getLogger(BOQATest.class);

    @BeforeClass
    public static void loadHPO() throws InterruptedException, IOException, URISyntaxException
    {
//        hpo = new Datafiles(
//                new File(ClassLoader.getSystemResource("human-phenotype-ontology.obo.gz").toURI()).getCanonicalPath(),
//                new File(ClassLoader.getSystemResource("phenotype_annotation.omim.gz").toURI()).getCanonicalPath());

    }
    //Generates num diseases and associates them with terms in the graph
    //Diseases get 2 to 18 annotations, mirroring real life.
    public static AssociationContainer generateAnnotations (int num, SlimDirectedGraphView<Term> slim)
    {
        //Set the random to a seed
        Random rnd = new Random(2);
        AssociationContainer assocs = new AssociationContainer();
        for (int i = 0; i<num; i++)
        {

            ByteString item = new ByteString("item" + i);

            // Association a = new Association(item,slim.getVertex(10).getIDAsString());
            // assocs.addAssociation(a);

            for (int j = 0; j < rnd.nextInt(16) + 2; j++) {
                Term t;
                do {
                    t = slim.getVertex(rnd.nextInt(slim.getNumberOfVertices())); //randomly select a vertex
                    //keeps doing this til it gets a non-obsolete vertex
                } while (t.isObsolete());

                Association a = new Association(item, t.getIDAsString());
                //System.err.println(a.toString());
                //print(a);
                assocs.addAssociation(a); //this seems to not hve any effect on BOQA... (nvm, it is used inb boqa.setup)
            }
        }
        return assocs;
    }
    @Test
    public void testLargeNumberOfItems() throws IOException, OBOParserException, URISyntaxException
    {


        final ReducedBoqa boqa = new ReducedBoqa();
        //get the file, then get its canonical path
        AssociationContainer assocs;


        URL resource = ClassLoader.getSystemResource("hp.obo.gz");
        if (resource==null)
        {

            throw new NullPointerException("Couldn't find it!");
        }
        URI resourceURI = resource.toURI();
        File hpo_file = new File(resourceURI);
        String final_path = hpo_file.getCanonicalPath();
        //.toURI();

        OBOParser hpoParser = new OBOParser(
                final_path);
        hpoParser.doParse();

        //blackbox: it gets all the terms (in the HPO)
        TermContainer tc = new TermContainer(hpoParser.getTermMap(), hpoParser.getFormatVersion(), hpoParser.getDate());
        Ontology ontology = new Ontology(tc);
        SlimDirectedGraphView<Term> slim = ontology.getSlimGraphView();
        assocs = generateAnnotations(10, slim);

        //pseudo:
        //boqa.setup
        //boqa.assignMarginals (get best score)
        //inference step: do some sampling to see which might be best (like in Monte Carlo tree search)

        //Run BOQA once to get the initial guesses.
        ArrayList<String> initial_guesses = null;

        boqa.setup(ontology, assocs);

        Observations o = new Observations();
        o.observations = new boolean[boqa.getOntology().getNumberOfTerms()];
        long start = System.nanoTime();
        this.logger.info("Calculating");
        boqa.assignMarginals(o, false, 1);
        long end = System.nanoTime();

        this.logger.info(((end - start) / 1000 / 1000) + "ms");
    }
}

