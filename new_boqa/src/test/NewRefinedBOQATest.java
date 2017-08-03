package test;

import com.sun.javafx.image.IntPixelGetter;
import ontologizer.association.Association;
import ontologizer.association.AssociationContainer;
import ontologizer.association.Gene2Associations;
import ontologizer.enumeration.GOTermEnumerator;
import ontologizer.go.*;
import ontologizer.types.ByteString;
import calculation_john_pack.Observations;
import calculation_john_pack.ReducedBoqa;
import org.junit.jupiter.api.Test;
import sonumina.math.graph.SlimDirectedGraphView;

import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.*;

/**
 * Created by johnp on 2017-06-29.
 */
public class NewRefinedBOQATest {
    public boolean getObservation(int index, ReducedBoqa rb) {
        ; //as long as only Terns are vertices, we will be fine.

        //However, these results must correspond with

        //rb.getOntology().getSlimGraphView().get;
        //Randomly generate a number. If it is less than or equal to the probability given, then we report that the
        //phenotype was observed (return True). Else return false.
        //mapping is between bytestring and associations
        //given a termid, we want to get the term out

        Term t;
        Term target_pheno_to_check = rb.slimGraph.getVertex(index);
        double probs = 0;
        for (TermID tt: trueDiseaseMapping.getAssociations())
        {

            t = rb.getOntology().getTerm(tt);
            //check if the term is a descendant of the input index term
            //Because you are your own ancestor, we do not need to deal with the self case specially.
            if (rb.slimGraph.isAncestor(target_pheno_to_check,t ))
            {
                //DO NOT use this nakedly! (it will sum to > 1!)
                probs+= pheno_disease_freq.get(t).get(trueDiseaseMapping.name()) ;
            }
            ///This just has the frequency class (1-5)
            //for now, let us use it directly as P(ph|D)




        }
        Random r = new Random(21);
        //this heavily rewards those with long parental chains. however we assume that
        //we only have the most specific. however, that is not justified!
        return r.nextDouble() < probs;
    }

    public boolean trueDiseaseInTopNDiseases(String target, List<String> top) {

        for (String s : top) {
            if (target.equals(s))
                return true;
        }
        return false;
    }

    //Generates num diseases and associates them with terms in the graph
    //Diseases get 2 to 18 annotations, mirroring real life.
    public AssociationContainer generateAnnotations(int num, SlimDirectedGraphView<Term> slim
    ) {
        pheno_disease_freq = new HashMap<>();

        //initialize all the inner hashmaps:

        //Set the random to a seed
        Random rnd = new Random(2);
        AssociationContainer assocs = new AssociationContainer();
        for (int i = 0; i < num; i++) {

            ByteString item = new ByteString("item" + i);

            for (int j = 0; j < rnd.nextInt(16) + 2; j++) {
                Term t;
                do {
                    t = slim.getVertex(rnd.nextInt(slim.getNumberOfVertices())); //randomly select a vertex
                    //keeps doing this til it gets a non-obsolete vertex
                } while (t.isObsolete());
                Association a = new Association(item, t.getIDAsString());
                //here we are simply required to remember what TID and temr was.
                //we want to be able to update the indices, based on the vertex2ancestor info from before, and
                //we CAN do that!
                //since for example, the interface between
                //let us make it a mapping between terms, and items and frequencies
                if (pheno_disease_freq.containsKey(t)) {
                    pheno_disease_freq.get(t).put(item, rnd.nextInt(freq_categories.length)); //TODO, make this vary based on the length of the freq array

                } else {
                    pheno_disease_freq.put(t, new HashMap<ByteString, Integer>());
                    pheno_disease_freq.get(t).put(item, rnd.nextInt(freq_categories.length)); //these correspond to the frequency classes

                }

                assocs.addAssociation(a);
            }
        }
        return assocs;
    }

    static public ArrayList<String> getTopDiseasesAsStrings(final ReducedBoqa.Result res) {
        // All of this is sorting diseases by marginals
        Integer[] order = new Integer[res.size()];
        for (int i = 0; i < order.length; i++) {
            order[i] = i;
        }

        //we should be able to get index2term

        Arrays.sort(order, new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
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

        //order.length/2
        for (int i = 0; i < 10; i++) {
            int id = order[i]; //presumably, order[i] is now in order from lowest to msot score
            results.add("item" + id); //bytestrings can be immediately constructed from this
            //"Disease "+ id + "\t"  + "Probs"  + res.getScore(id) ); //all amrginals are the same...
        }

        return results;
    }

    /***
     * Postcondition: we expect obs[getBestPhenotype] to be set to true later
     * @param rb
     * @param phenotype_frequencies
     * @return
     */
    public int getBestPhenotype(ReducedBoqa rb, double[] phenotype_frequencies) {
        //why not just maintain the observed array from before?
        //int [] top_phenotypes; //TODO keep the n largest phenotypes
        //we can use a min heap here! we want to maintain the largest k elements in a stream
        //(simply check the min elemnt, and if its larger, kick out one of the elements (replace)

        //we have the index, but not about the TermID
        int best_phenotype_index = 0;
        double best_phenotype_value = 0;
        double temp = 0;
        ReducedBoqa.Result res;
        ArrayList<Integer> turnedOn;
        for (int i = 0; i < rb.o.observations.length; i++)//(int termInd : rb.o.observations)
        {
            if (!rb.o.observations[i]) {

                //actually, we must set ALL the above nodes to true too!
                //actually this is almost NEVER true (since negatives don't tell us anything more)
                rb.o.real_observations[i] = true;
                //we need this next line so that assignMarginals makes sense
                rb.o.observations[i] = true;
                turnedOn = setAncestors(rb, i);
                System.out.println("starting assignmarginals");
                long start = System.nanoTime();

                res = rb.assignMarginals(rb.o, false, 1);
                System.out.println("done assignmarginals. Took" + (System.nanoTime()-start));

                //after assignMarginals, the scores array is updated
                //boqa roll back: undoes the last one

                //if our current best worse than this new one, update it
                System.out.println("starting scoring and misc");
                start = System.nanoTime();

                //Should be memoized (in order to use the previous value (from the previous iteration), we need to
                //write it down somewhere!
                if (best_phenotype_value <
                        (temp = scoringFunction(res, rb) * phenotype_frequencies[i])) //pass in rb in case we need ot infer other things
                {       //weight on the phenotype_frequency! (how likely it is to be true)
                    best_phenotype_index = i;
                    best_phenotype_value = temp;

                    //index to termID is not well supported.
                    //since these things only make sense wrt a graph
                    //also termid and term are not well related, we really just want
                    //term,since it encapsulates everything
                    //hence we should use 3termcontainer to send it into

                }
                System.out.println("done assignmarginals. Took" + (System.nanoTime()-start));

                //undoes the setting action
                rb.o.real_observations[i] = false;
                rb.o.observations[i] = false;
                rollback(rb, turnedOn);
            }

        }
        //rb.getOntology().getSlimGraphView().getVertex()//this can recover the Term if need be
        return best_phenotype_index;
    }

    public double scoringFunctionOnArray(double [] freqs) {
        double score = 0;

        for (double marg : freqs) {
            score += marg * marg;
        }

        return score;
    }

    //pheno rows, disease cols
    public int multiGetBestPhenotype(double [][] phenoDiseaseDist, ReducedBoqa rb)
    {
        int best_phenotype_index = 0;
        double best_phenotype_value = Double.NEGATIVE_INFINITY;
        double temp;
        double val;
        for (int i = 0; i<phenoDiseaseDist.length; i++)
        {
            //assert rb.o.observations.length =phenoDiseaseDist.length
            //We cannot return pick a phenotype twice
            if (!rb.o.observations[i]) {
                if (best_phenotype_value < (temp = phenotype_frequencies[i] * scoringFunctionOnArray(phenoDiseaseDist[i]))) {
                    best_phenotype_index = i;
                    best_phenotype_value = temp;
                }
            }

        }
        return best_phenotype_index;
    }

    //Represents the disease-phenotype frequency annotation data.
    //I1: Disease
    //I2: Phenotype
    //I3: Frequency Category
    double[] phi_phenotype_frequencies;
    double[] phenotype_frequencies;
    double[] disease_frequencies; //this is actually just BOQA's marginals
    HashMap<Term, HashMap<ByteString, Integer>> pheno_disease_freq;
    double [] freq_categories = {0.2,0.4,0.6,0.8,1};
    ByteString trueDisease;
    Set<Term> trueDiseasePhentoypes; //perhaps an association container might have been best
    //AssociationContainer;
    Gene2Associations gx = new Gene2Associations(new ByteString("aa"));
    Gene2Associations trueDiseaseMapping;

    /***
     * TODO skip phenotypes that have already been observed
     */
    //it would need to compute the intrinsics separately

    public void computePhiPhenotypeFrequencies(ReducedBoqa rb) {
        long start = System.nanoTime();
//        System.out.println("starting compute of Phi");

        double temp;
        double counter = 0;
        int pheno;
//        System.out.println("size is" + pheno_disease_freq.size());

        SlimDirectedGraphView<Term> graph = rb.getOntology().getSlimGraphView();

        for (Term pheno_term : pheno_disease_freq.keySet()) {
            counter++;

            temp = 0;
            //alternatively:
            //rb.slimGraph.getVertexIndex()

            ByteString bs;
            //ideally, we have a bytestring->index
            //if we are just doing phenotype to disease, then we can directly use these elements
            for (Map.Entry annotation : pheno_disease_freq.get(pheno_term).entrySet()) {

                //does P(D)*P(I|D)
                temp += disease_frequencies[rb.item2Index.get( annotation.getKey())] * freq_categories[(Integer) annotation.getValue()];
                //now, we need to use the index of the Bytestring now!


            }
            pheno = graph.getVertexIndex(pheno_term);
            //pheno = rb.getOntology().getSlimGraphView().getVertexIndex(pheno_term);
            phi_phenotype_frequencies[pheno] = temp;

        }
//        System.out.println("finish compute of Phi");
//        System.out.println("done, took " + (System.nanoTime()-start));

    }
    //computes it from the REST of the array.
    public void computePhenotypeFrequencies(ReducedBoqa rb)
    {
//        long start = System.nanoTime();
//        System.out.println("starting compute of phen-freq");
        double temp;
        int pheno;
        SlimDirectedGraphView<Term> graph = rb.getOntology().getSlimGraphView();
        //either fill it left to right, or using the order found in the hashmap
        for (Term pheno_term : pheno_disease_freq.keySet()) {
            temp = 0;
            pheno = graph.getVertexIndex(pheno_term);
            //The self component
            temp+= phi_phenotype_frequencies[pheno];

            //The descendant component
            for (Integer i: rb.term2Descendants[pheno])
            {
                //actually don't add in these ones (we only want the COMPONENTS)
                //to get the components would be to just call this function again!
                temp+= phi_phenotype_frequencies[i];

            }


            phenotype_frequencies[pheno] =temp;
        }
//        System.out.println("finish compute of phen-freq");
//        System.out.println("done, took " + (System.nanoTime()-start));





    }

    //sets all ancestors of a node in the boqa instance observations to on
    //only call this function if you are recording a positive, as when checking for
    //next best phenotype, or recording a positive phenotype

    //note that is like intersection of linked lists! when we have a common node,
    //we can just stop!

    //either return an array of modified (safe)
    //alternatievely, is rollback we can compute the intersection of nodes
    //however, this has the possiblity of turning OFF nodes that were turnt on by others
    //for example, the root node will be turned on by anything, but rolling back in this
    //way will turn it off (undesired), this violates true path rule etc.
    public static ArrayList<Integer> setAncestors(ReducedBoqa rb, int index) {
        //list of terms that were turned on by setting the index to 1 (or true)
        ArrayList<Integer> turnedOnTerms = new ArrayList<>();
        for (int anc : rb.term2Ancestors[index]) {
            //essentially, return the set that WAS changed by the setting of the nodes.
            //only set unobs false ancestors
            if (!rb.o.observations[anc] && !rb.o.real_observations[anc]) {


                rb.o.observations[anc] = true; //understated line tbh. -- this means we treat propagated obs as TRUE obs
                rb.o.real_observations[anc] = true;
                //only do this if we propagate positives up too
                //do NOT propagate negatives up (since this is not how the true path rule works)


                //there is no way to do the following now


                //if its not observed, is there any chance that it was already set to true?
                //Yes, if it was inferred. Note that, the wya our code works is that it will still be added as a real obs even
                //if it was jsut the kid that was checked
                //hence only add it if it was also false, which will be the case msotly

                    //the only way it could be true is if it is observed, or one of its children was observed
                    //hence the first if condition (checking that it wasn't observed) is good enough
                    //Hence we guarantee the state was off to start off with
                    //essentially, if if it is observed, there is only one case, where it is true
                turnedOnTerms.add(anc);//otherwise this would add index a lot!

            }


        }

        return turnedOnTerms; //this will be input to the rollback

    }

    public void rollback(ReducedBoqa rb, ArrayList<Integer> turntOn) {
        for (Integer i : turntOn) {
            rb.o.observations[i] = false;
            rb.o.real_observations[i] = false;
            //this is valid since condition for entry to turntOn is false in BOTH

        }

    }

    public double scoringFunction(ReducedBoqa.Result result, ReducedBoqa rb) {
        double score = 0;

        for (double marg : result.marginals) {
            score += marg * marg;
        }

        return score;
    }


    public void generateTrueDisease(SlimDirectedGraphView<Term> slim, AssociationContainer assocs)
    {
        Random rnd = new Random(41); //this is our true disease
        for (int j = 0; j < rnd.nextInt(16) + 2; j++) {
            Term t;
            do {
                t = slim.getVertex(rnd.nextInt(slim.getNumberOfVertices())); //randomly select a vertex
                //keeps doing this til it gets a non-obsolete vertex
            } while (t.isObsolete());

            //how it works, is we store key:value pairs, regardless of whether key already exists
            //(hence not like a dictionary)
            Association trueDiseasePhenotype = new Association(trueDisease, t.getIDAsString());

            //we have a huge array of size (#items by #terms large)
            //alternatively just maintain an arraylist of <item, arraylist <terms>>
            //we already have that, we just want to reconstruct new ByteStrings from
            //the list of diseases we get back
            //we just need to get the top numbers, from the result.
            //take the top n/2 results from the Result

            //assocs.get(item).getAssociations();
            //System.err.println(a.toString());
            //print(a);

            if (pheno_disease_freq.containsKey(t)) {
                pheno_disease_freq.get(t).put(trueDisease, rnd.nextInt(freq_categories.length));

            } else {
                pheno_disease_freq.put(t, new HashMap<ByteString, Integer>());
                pheno_disease_freq.get(t).put(trueDisease, rnd.nextInt(freq_categories.length)); //these correspond to the frequency classes

            }
            trueDiseaseMapping.add(trueDiseasePhenotype);


            assocs.addAssociation(trueDiseasePhenotype); //this seems to not hve any effect on BOQA... (nvm, it is used inb boqa.setup)
        }

        System.out.println(trueDiseaseMapping.getAssociations());


    }

    //approximateArrayList by printing the indices of where it is true only
    public List<Integer> getIndicesOfTrue(boolean [] arr)
    {
        List<Integer> l = new ArrayList<>();
        for (int i = 0; i < arr.length; i++){

            if (arr[i]){
                l.add(i);
            }
        }

        return  l;
    }

    public void print_find_ancestors_of_trueDisease(ReducedBoqa rb, TermContainer tc)
    {
        for (TermID ti :trueDiseaseMapping.getAssociations())
        {
            System.out.println("the ancestors of " + ti + " are "
                    + rb.slimGraph.getAncestors(tc.get(ti)));

        }
    }

    @Test
    public void testConvergence() throws IOException, OBOParserException, URISyntaxException {
        int num = 10000;
        final ReducedBoqa boqa = new ReducedBoqa();
        //boqa.getOntology().
        //boqa.getOntology().getTerm() //FROM THE TERMID, we can recover the terms, and also recover the indexes?
        //yes, we are sure that the ints produced are the same ints as used in the Boqa.java (since, we get the
        //vertex2ancestors just immediately from boqa)

        //get the file, then get its canonical path
        AssociationContainer assocs;


        URL resource = ClassLoader.getSystemResource("hp.obo.gz");
        ClassLoader cl = ClassLoader.getSystemClassLoader();

        URL[] urls = ((URLClassLoader)cl).getURLs();

        for(URL url: urls){
            System.out.println(url.getFile());
        }
        if (resource == null) {
            throw new NullPointerException("Couldn't find it!");
        }
        URI resourceURI = resource.toURI();
        File hpo_file = new File(resourceURI);
        String final_path = hpo_file.getCanonicalPath();
        //.toURI();

        OBOParser hpoParser = new OBOParser(
                final_path);
        hpoParser.doParse();
        //this line illustrates the 301 lesson: getting an object, then being able to modify it
        //we should return a copy of it or something
        long start = System.nanoTime();
        removeObsoleteTerms(hpoParser, start);

        //blackbox: it gets all the terms (in the HPO)
        //getTermMap returns a list of all terms!
        TermContainer tc = new TermContainer(hpoParser.getTermMap(), hpoParser.getFormatVersion(), hpoParser.getDate());


        Ontology ontology = new Ontology(tc);
        SlimDirectedGraphView<Term> slim = ontology.getSlimGraphView();


        assocs = generateAnnotations(num, slim);

        trueDisease = new ByteString("item" + num);
        trueDiseaseMapping = new Gene2Associations(trueDisease);

        generateTrueDisease(slim, assocs);


        Observations o = new Observations();
        int numberOfTerms = ontology.getNumberOfTerms();
        o.observations = new boolean[numberOfTerms];
        o.real_observations = new boolean[numberOfTerms];
        //Run BOQA once to get the initial guesses.
        ArrayList<String> initial_guesses = null;

        boqa.setup(ontology, assocs);
        boqa.setO(o);

        int steps = 0;
        double increment = 0.01;
        boolean discovered = false;
        phenotype_frequencies = new double[numberOfTerms]; //alternatively, just copy over the
        //array length from the item2ancestors for example
        phi_phenotype_frequencies = new double[numberOfTerms];


        //initalization/first step stuff
        ReducedBoqa.Result res=new ReducedBoqa.Result();
        steps++;
        boqa.setInitial_beta(increment * steps);
        res = boqa.assignMarginals(o, false, 1);
        disease_frequencies = res.marginals;
        long total = System.nanoTime();
        print_find_ancestors_of_trueDisease(boqa, tc);
        while (!discovered) {

            if (steps >1)
            {
                for (int i  = 0; i < boqa.multiDiseaseDistributions.length-1; i++)
                {
                    for (int j = i+1; j < boqa.multiDiseaseDistributions.length; j++)
                    {
                        if (Arrays.equals(boqa.multiDiseaseDistributions[i],
                                boqa.multiDiseaseDistributions[j]))
                        {
                            System.out.println("These arrays are equal " + i + "  " + j);
                        }

                    }
                }

            }
            total = System.nanoTime();
            System.out.println("this is step" + steps);
            System.out.println("These are the ones checked");
            System.out.println(getIndicesOfTrue(boqa.o.observations));
            System.out.println("These are the one present");
            System.out.println(getIndicesOfTrue(boqa.o.real_observations));
            //boqa.setInitial_beta(boqa.getInitial_beta()-boqa.getInitial_beta()/30);
            //Alternatively, we could jsut have the difference too (inital beta-experimental beta)


            //assign marginals with the new o.


            System.out.println("starting multishot");
            start = System.nanoTime();
            boqa.assignMultiShotMarginals(o,false,8);

            //now that we have the matrix of probabiltiies, we can look up the pheno-freqs and weight it accordingly
            //this process will probably take some time as well (+20s)

            System.out.println("done multishot. Took" + (System.nanoTime()-start));
//            System.exit(0);
            //TODO doesn't need to be called on every loop

            //this returns an int array[], where each elt is the prob of item with that index
            //we can introconvert if we have the index2term for example
            //it is interesting since they almost exclusively interface with the int id representations
            //in the BOQA, yet here use termIds and such
            //the key bridge is itemEnumerator
//        ItemEnumerator itemEnumerator = ItemEnumerator.createFromTermEnumerator(this.termEnumerator);
//
//        itemEnumerator.getTermsAnnotatedToTheItem(item);




            computePhiPhenotypeFrequencies(boqa);
            computePhenotypeFrequencies(boqa);

            //update with the results of the new boqa run
            System.out.println("starting pheno check");
            start = System.nanoTime();
            int phenotype_to_check = multiGetBestPhenotype(boqa.multiDiseaseDistributions,boqa); //in here we do all the phenotype checks

            System.out.println("done pheno check. Took" + (System.nanoTime()-start));
            //This allows us to go from TermID->index, but what about the other way>
            //int index =boqa.slimGraph.getVertexIndex(boqa.getOntology().getTerm(phenotype_to_check));

            int index = phenotype_to_check; //much simpler
            System.out.println("this is the one im checking" + index);
            System.out.println("this is the HPO im checking" + boqa.slimGraph.getVertex(index));

            if (steps>=24){
                Term tempt = boqa.slimGraph.getVertex(index);
                boqa.slimGraph.getAncestors(tempt);
                System.out.println(boqa.slimGraph.getAncestors(tempt));
                System.out.println(boqa.slimGraph.getDescendants(tempt));
            }

            boolean present_or_not = getObservation(index,boqa);

            //get input from physician, and update the observations object
            //ALL ancestors must be updated as well!
            //o doesn't have this information, so we need to full info from boqa:
            if (present_or_not){
                for (int anc : boqa.term2Ancestors[index]) {
                    o.observations[anc] = true; //@TODO we assume observing hte child is the same as observing the parent
                    o.real_observations[anc] = true; //only do this if we propagate positives up too
                    //do NOT propagate negatives up (since this is not how the true path rule works)


                    //if its negative, should we

                }}
            //this should all be abstracted to another function!
            o.real_observations[index] = present_or_not;



            o.observations[index] = true; //recall the new meanings: observations means whether it was
            //checked, while the arraylist determines whether it was true or not
            //this has been deprecated

            //repeating this process should segregate everything
            //however, this can and SHOULD happen normally (probability of seeing something false
            //repeatedly is low. this is not a healing love, this is a wicked fantasy
            //test shoudl work:

            //updating steps for the next one:
            steps++;
            boqa.setInitial_beta(increment * steps);
            res = boqa.assignMarginals(o, false, 1);
            disease_frequencies = res.marginals;


//            System.out.println(java.util.Arrays.toString(res.marginals));
//            System.out.println(java.util.Arrays.toString(res.scores));
            printTopDisease(res);

            //sorts the array, by getScore and takes the top N
            System.out.println(getTopDiseasesAsStrings(res));
            initial_guesses = getTopDiseasesAsStrings(res); //we have essentially the top ids now
            //from the ids, we can get the mappings they have

            //now, we recompute the marginals.
            //o.setValue()if ()
            if (trueDiseaseInTopNDiseases(trueDisease.toString(), initial_guesses)) {
                discovered = true;
                System.out.println("we are finishedd!");
            }

            System.out.println("done loop iter. Took" + (System.nanoTime()-total));
        }
    }

    private void printTopDisease(ReducedBoqa.Result res) {
        //Computes max element and the index it occurs at
        double max = -Double.MAX_VALUE;
        int max_ind = 0;
        for (int i = 0; i < res.marginals.length; i++) {
            if (res.marginals[i] > max) {

                max = res.marginals[i];
                max_ind = i;
            }

        }

        System.out.println("max_ind is " + max_ind + " max is " + max);
    }

    private void removeObsoleteTerms(OBOParser hpoParser, long start) {
        Set<Term> terms = hpoParser.getTermMap();
        Iterator<Term> iter= terms.iterator();

        while (iter.hasNext())
        {
            Term t = iter.next();

            if (t.isObsolete())
            {
                //System.out.println(t);
                iter.remove();
            }
        }
        System.out.println("Removal of obsolete took" + (System.nanoTime()-start));
    }
}
