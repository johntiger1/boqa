package sonumina.boqa.tests;

/**
 * Created by johnchen on 18/05/17.
 */
import java.net.URI;

import ontologizer.association.*;
import ontologizer.benchmark.Datafiles;
import ontologizer.enumeration.GOTermEnumerator;
import ontologizer.enumeration.ItemEnumerator;
import ontologizer.go.*;
import ontologizer.types.ByteString;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import sonumina.boqa.BOQABenchmark;
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
import sun.font.TrueTypeFont;

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
    public static AssociationContainer generateAnnotations (int num, SlimDirectedGraphView<Term> slim
                                                            )
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

    static public ArrayList<String> getTopDiseasesAsByteStrings(final ReducedBoqa.Result res)
    {
        // All of this is sorting diseases by marginals
        Integer[] order = new Integer[res.size()];
        for (int i = 0; i < order.length; i++) {
            order[i] = i;
        }

        //we should be able to get index2term

        Arrays.sort(order, new Comparator<Integer>()
        {
            @Override
            public int compare(Integer o1, Integer o2)
            {
                if (res.getScore(o1) < res.getScore(o2)) {
                    return 1;
                }
                if (res.getScore(o1) > res.getScore(o2)) {
                    return -1;
                }
                return 0;
            }
        }); //order[i] will be sorted according to the comparator, so [1..n] will become [3,2,67,1,..]

        // Get top 20 results
        ArrayList<String> results = new ArrayList<String>();
        //for (int i = 0; i < 20 && i<order.length; i++) {
        for (int i = 0; i<order.length/2; i++) {
            int id = order[i]; //presumably, order[i] is now in order from lowest to msot score
            results.add( "item" + id ); //bytestrings can be immediately constructed from this
                    //"Disease "+ id + "\t"  + "Probs"  + res.getScore(id) ); //all amrginals are the same...
        }

        return results;
    }

    @Test
    //this test
    //in reality, it seems like we really only need one giant test of correctness
    //the test needs to simulate the physician entering phenotypes in and getting the diseases out
    //we can simulate it with our test, always telling the truth:
    //"inverted boqa", enter the diseases and see the symptoms? but that is already the canonical way of diseases

    //"lying akinator" you are allowed to tell one lie (or total # of answers * lie_rate)
    //
    public void testConvergence() throws IOException, OBOParserException, URISyntaxException
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
        int num = 1000;
        assocs = generateAnnotations(num, slim);

        ByteString item = new ByteString("item" + num);
        Random rnd = new Random(3); //this is our true disease
        for (int j = 0; j < rnd.nextInt(16) + 2; j++) {
            Term t;
            do {
                t = slim.getVertex(rnd.nextInt(slim.getNumberOfVertices())); //randomly select a vertex
                //keeps doing this til it gets a non-obsolete vertex
            } while (t.isObsolete());

            //how it works, is we store key:value pairs, regardless of whether key already exists
            //(hence not like a dictionary)
            Association trueDiseasePhenotype = new Association(item, t.getIDAsString());
            System.out.println(trueDiseasePhenotype.getEvidence());

            //we have a huge array of size (#items by #terms large)
            //alternatively just maintain an arraylist of <item, arraylist <terms>>
            //we already have that, we just want to reconstruct new ByteStrings from
            //the list of diseases we get back
            //we just need to get the top numbers, from the result.
            //take the top n/2 results from the Result

            //assocs.get(item).getAssociations();
            //System.err.println(a.toString());
            //print(a);
            assocs.addAssociation(trueDiseasePhenotype); //this seems to not hve any effect on BOQA... (nvm, it is used inb boqa.setup)
        }
        //save one of the things
        //Term t1 = boqa.termEnumerator.getAnnotatedGenes(new TermID());
        ByteString disease_2 = new ByteString(("item" + num));//otherwise we can just STORE ALL OF THESE


        Gene2Associations s = assocs.get(item); //should do the inverse op and get the annotations
        Gene2Associations s2 = assocs.get(disease_2); //should do the inverse op and get the annotations
        //we can make membership queries to the g2a
        //use the termcontainers etc. to make queries

        //Note that most of our job is actually building the REVERSE thing:
        s.getAssociations(); //gets all the associations
        //we can just maintain an array (of all the terms) and build it up as we go through the
        //associations for EACH disease.

        //

        // we need to obtain termIDs somehiow, from the BOQA, or the graph
        //perhaps from the ranked list that we get back
        //result
        //s.containsID(91);




        //now we know which is the true disease, as well as what it annotates


        //pseudo:
        //boqa.setup
        //boqa.assignMarginals (get best score)
        //inference step: do some sampling to see which might be best (like in Monte Carlo tree search)

        //Run BOQA once to get the initial guesses.
        ArrayList<String> initial_guesses = null;

        boqa.setup(ontology, assocs);
        //provides a reverse mapping (from HPO term to disease)
        GOTermEnumerator x = boqa.termEnumerator;
        //x.

        //however, we

        Observations o = new Observations();
        o.observations = new boolean[boqa.getOntology().getNumberOfTerms()];
//        for (int i = 0; i < o.observations.length; i++)
//        {
//
//            Random rand = new Random(1);
//            o.observations[i] = (rand.nextInt(2) > 0) ? true :false;
//            //(rand.nextInt(2) > 0) ? o.observations[i] = true : o.observations[i] = false;
//        }
        //o.observations[10] = true;  //has no effect

        long start = System.nanoTime();
        this.logger.info("Calculating");
        ReducedBoqa.Result res = boqa.assignMarginals(o, false, 1);
        //this returns an int array[], where each elt is the prob of item with that index
        //we can introconvert if we have the index2term for example
        //it is interesting since they almost exclusively interface with the int id representations
        //in the BOQA, yet here use termIds and such
        //the key bridge is itemEnumerator
//        ItemEnumerator itemEnumerator = ItemEnumerator.createFromTermEnumerator(this.termEnumerator);
//
//        itemEnumerator.getTermsAnnotatedToTheItem(item);

        System.out.println(java.util.Arrays.toString(res.marginals)); //doesn't res store ONE thing tho?
        System.out.println(java.util.Arrays.toString(res.scores));

        double max = -Double.MAX_VALUE;
        int max_ind = 0;
        for (int i = 0 ; i < res.marginals.length; i++)
        {
            if (res.marginals[i] > max){

                max = res.marginals[i];
                max_ind = i;
            }

        }

        System.out.println("max_ind is " + max_ind+ " max is " + max);
        System.out.println(getTopDiseasesAsByteStrings(res));
        initial_guesses = getTopDiseasesAsByteStrings(res); //we have essentially the top ids now
        //from the ids, we can get the mappings they have
        Integer [] termcounts = new Integer[boqa.getOntology().getNumberOfTerms()];
        for (int i = 0; i < termcounts.length; i++){

            termcounts[i] = 0;
        }
        for (String str: initial_guesses)
        {
            Gene2Associations temp = assocs.get(new ByteString(str)); // we could also have memoized from before
            //temp.getAssociations(). //we would like an iterator
//            for (TermID a : temp.getAssociations())
//            {
//                termcounts[a.id]++; //consider using the termcontainer,
//                //just like how I am using the assocs
//            }

            for (TermID tid : temp.getAssociations()) {
                System.out.println("vale is " + termcounts[boqa.slimGraph.getVertexIndex(boqa.getOntology().getTerm(tid))]);
                termcounts[boqa.slimGraph.getVertexIndex(boqa.getOntology().getTerm(tid))]++ ;
            }

            //this must be recursive to get ALL THE PHENOTYPES

            //termcounts[0]++;
            Arrays.sort(termcounts, Collections.<Integer>reverseOrder());
            //our own hashing function:

        }

        //most likely the terms are also just ints and so can be
        HashMap<Association, Integer> termCounts;



        //get input from physician, and update the observations object
        o.recordObs(10, true);
        //o.setValue();
    }


    //Reduces the BOQA result accordingly,
    //one thing we could do is return only a slice of the array.

    public ReducedBoqa.Result reduceDiseases(final ReducedBoqa.Result res)
    {

        // All of this is sorting diseases by marginals
        Integer[] order = new Integer[res.size()];
        for (int i = 0; i < order.length; i++) {
            order[i] = i;
        }

        //we should be able to get index2term
        //sorts the [1..n] according to how their getScores(i) are
        //unfortunately this does not persist!
        Arrays.sort(order, new Comparator<Integer>()
        {
            @Override
            public int compare(Integer o1, Integer o2)
            {
                if (res.getScore(o1) < res.getScore(o2)) {
                    return 1;
                }
                if (res.getScore(o1) > res.getScore(o2)) {
                    return -1;
                }
                return 0;
            }
        }); //order[i] will be sorted according to the comparator, so [1..n] will become [3,2,67,1,..]

        //res.

        return res;

    }
    //updates the internal state of BOQA with the new observations, etc.
    public void updateBOQA(int index)
    {
        //TODO: add an internal "true observed" state to BOQA



    }

    //associations are between terms and items (diseases)
    //roughly, they map HPO term #s to an Bytestrign (which includes an integer)
    public Term getBestPhenotype(ReducedBoqa.Result res)
    {


        return null;
    }

    @Test
    public void testLargeNumberOfItems() throws IOException, OBOParserException, URISyntaxException
    {

        //Testing framework:
        //first, generate the annotations and diseases
        //then pick a random disease, and save its state
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
        assocs = generateAnnotations(25, slim);

        //pseudo:
        //boqa.setup
        //boqa.assignMarginals (get best score)
        //inference step: do some sampling to see which might be best (like in Monte Carlo tree search)

        //Run BOQA once to get the initial guesses.
        ArrayList<String> initial_guesses = null;

        boqa.setup(ontology, assocs);

        Observations o = new Observations();
        o.observations = new boolean[boqa.getOntology().getNumberOfTerms()];
//        for (int i = 0; i < o.observations.length; i++)
//        {
//
//            Random rand = new Random(1);
//            o.observations[i] = (rand.nextInt(2) > 0) ? true :false;
//            //(rand.nextInt(2) > 0) ? o.observations[i] = true : o.observations[i] = false;
//        }
        //o.observations[10] = true;  //has no effect

        long start = System.nanoTime();
        this.logger.info("Calculating");
        ReducedBoqa.Result res = boqa.assignMarginals(o, false, 1);
        System.out.println(java.util.Arrays.toString(res.marginals)); //doesn't res store ONE thing tho?
        System.out.println(java.util.Arrays.toString(res.scores));

        double max = -Double.MAX_VALUE;
        int max_ind = 0;
        for (int i = 0 ; i < res.marginals.length; i++)
        {
            if (res.marginals[i] > max){

                max = res.marginals[i];
                max_ind = i;
            }

        }

        System.out.println("max_ind is " + max_ind+ " max is " + max);
        System.out.println(getTopDiseasesAsByteStrings(res));
        //for (double t: res.)
        //write a method that keeps track of the top 10 scores
        //use concept of lower and upper bound
        //instead: use selection algorithm
        //however, it should be online (based on the api exported)
        //insert the first 10 unconditionally
        //then, for each element, check if it should be put in or not, --this is a linear time algorithm
        //but we must keep track of the max and min (i.e. go through the array and update hte max, min indices each time...)

        //find the 10th largest number, using quickselect
        //then, we shall have the 10 larger numbers on one side and we can just return that
        //easier way is to just sort the array and then take the top n elements, except we need a reference
        //to previous. This issue is also a problem in using quickselect too.

        //we can just use parallel arrays though. for example, lookup[i] = pos_in__sorted_array


        long end = System.nanoTime();

        this.logger.info(((end - start) / 1000 / 1000) + "ms");
    }


    static public ArrayList<String> getTopDiseases(final BOQA.Result res)
    {
        // All of this is sorting diseases by marginals
        Integer[] order = new Integer[res.size()];
        for (int i = 0; i < order.length; i++) {
            order[i] = i;
        }
        System.out.println("this is what order has" + java.util.Arrays.toString(order));
        Arrays.sort(order, new Comparator<Integer>()
        {
            @Override
            public int compare(Integer o1, Integer o2)
            {
                if (res.getMarginal(o1) < res.getMarginal(o2)) {
                    return 1;
                }
                if (res.getMarginal(o1) > res.getMarginal(o2)) {
                    return -1;
                }
                return 0;
            }
        }); //order[i] will be sorted according to the comparator, so [1..n] will become [3,2,67,1,..]


        System.out.println("this is what order has" + java.util.Arrays.toString(order));
        // Get top 20 results
        ArrayList<String> results = new ArrayList<String>();
        for (int i = 0; i < 20 && i<order.length; i++) {
            int id = order[i];
            results.add( "Disease "+ id + "\t"  + "Probs"  + res.getMarginal(id) ); //all amrginals are the same...
        }

        return results;
    }

    @Test
    public void vanillaTestLargeNumberOfItems() throws IOException, OBOParserException, URISyntaxException
    {


        final BOQA boqa = new BOQA();
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
        assocs = generateAnnotations(25, slim);

        //pseudo:
        //boqa.setup
        //boqa.assignMarginals (get best score)
        //inference step: do some sampling to see which might be best (like in Monte Carlo tree search)

        //Run BOQA once to get the initial guesses.
        ArrayList<String> initial_guesses = null;
        boqa.setConsiderFrequenciesOnly(false);
        boqa.setup(ontology, assocs);

        Observations o = new Observations();
        o.observations = new boolean[boqa.getOntology().getNumberOfTerms()];

        long start = System.nanoTime();
        this.logger.info("Calculating");
        BOQA.Result res = boqa.assignMarginals(o, false, 1);
        System.out.println(getTopDiseases(res));
        //for (double t: res.)
        //write a method that keeps track of the top 10 scores
        //use concept of lower and upper bound
        //instead: use selection algorithm
        //however, it should be online (based on the api exported)
        //insert the first 10 unconditionally
        //then, for each element, check if it should be put in or not, --this is a linear time algorithm
        //but we must keep track of the max and min (i.e. go through the array and update hte max, min indices each time...)

        //find the 10th largest number, using quickselect
        //then, we shall have the 10 larger numbers on one side and we can just return that
        //easier way is to just sort the array and then take the top n elements, except we need a reference
        //to previous. This issue is also a problem in using quickselect too.

        //we can just use parallel arrays though. for example, lookup[i] = pos_in__sorted_array


        long end = System.nanoTime();

        this.logger.info(((end - start) / 1000 / 1000) + "ms");
    }
}

