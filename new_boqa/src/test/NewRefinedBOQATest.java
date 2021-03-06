package test;

import calculation_john_pack.ReducedConfiguration;
import com.github.phenomics.ontolib.formats.hpo.HpoDiseaseAnnotation;
import com.github.phenomics.ontolib.io.base.TermAnnotationParserException;
import com.github.phenomics.ontolib.io.obo.hpo.HpoDiseaseAnnotationParser;
import data_parsing.PATParser;
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
import weka.classifiers.Classifier;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.*;
import java.util.Map.Entry;

import static data_parsing.PATParser.ORPHA_name;
import static data_parsing.PATParser.getHPOToFreqMappingHardCoded;

/**
 * Created by johnp on 2017-06-29.
 */
public class NewRefinedBOQATest {

    public int  getFreeObs( ReducedBoqa rb) {
        Term t;
        int ind = 0;
        Random r = new Random ();
        for (TermID tt : trueDiseaseMapping.getAssociations()) {

            t = rb.getOntology().getTerm(tt);
            //check if the term is a descendant of the input index term
            //Because you are your own ancestor, we do not need to deal with the self case specially.
            ind = rb.getOntology().getSlimGraphView().getVertexIndex(t);
            if (r.nextDouble()>0.5)
            {
                System.out.println("I am giving you this one for free " + t);
                System.out.println("The index corresponding is " + ind);
                break;
            }

        }
        return ind;
    }

    public boolean getObservation(int index, ReducedBoqa rb, boolean noise) {
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
//                probs+= freq_categories[pheno_disease_freq.get(t).get(trueDiseaseMapping.name())] ;
//                TODO: we should probably add up NORMALIZED probabilities, and probably have a normalized freqeuncy*probability table
                probs+= pheno_disease_freq.get(t).get(trueDiseaseMapping.name()) ;
            }
            ///This just has the frequency class (1-5)
            //for now, let us use it directly as P(ph|D)



        }

        //standardize randomness
        Random r = new Random(120);
        if (noise)
        {
            //modify these probabilities
            //if present: make it not
            //
            if (probs==0)
            {
                System.out.println("THIS IS --FP!!!!!");
                return r.nextDouble() < rb.getALPHA_GRID()[0];
                //with some probability, observe it!


            }

//            OK, I see: by fluke, this means that (truthfully) at least one disease phenotype is a descendant of this term
            else
            {
                //with some probability don't observe;
                boolean temp = r.nextDouble() < rb.getExperimental_beta();

                if (temp)
                {
                    System.out.println("THIS IS NOISE--FN!!");
                    return false;

                }

                else {

                    return true;
                } //IOW, return complement of r.nextDouble()< ...
                //note this treats it as being associated = being present,
                //disregarding the frequencies!

            }



        }
        //if noise:



        //this heavily rewards those with long parental chains. however we assume that
        //we only have the most specific. however, that is not justified!
        return r.nextDouble() < 0.5;
    }

    public boolean trueDiseaseInTopNDiseases(String target, List<String> top) {

        for (String s : top) {
            if (target.equals(s))
                return true;
        }
        return false;
    }

    public void initializeHashmap(TermContainer tc)
    {
        pheno_disease_freq = new HashMap<>();
        for (Term t : tc)
        {
            pheno_disease_freq.put(t, new HashMap<ByteString, Integer>());
        }
    }

    //Generates num diseases and associates them with terms in the graph
    //Diseases get 2 to 18 annotations, mirroring real life.
    public AssociationContainer generateAnnotations(int num, SlimDirectedGraphView<Term> slim,
                                                    TermContainer tc
    )
    {




        initializeHashmap(tc);
        //initialize all the inner hashmaps:


        //Set the random to a seed
        Random rnd = new Random();
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
                    pheno_disease_freq.get(t).put(item, rnd.nextInt(freq_categories.length));

                } else {
                    pheno_disease_freq.put(t, new HashMap<ByteString, Integer>());
                    pheno_disease_freq.get(t).put(item, rnd.nextInt(freq_categories.length)); //these correspond to the frequency classes

                }

                assocs.addAssociation(a);
            }
        }

        return assocs;
    }

//    public static void printParallelSorted(double [] array,
//                                           TermContainer tc,
//                                           GOTermEnumerator gte,
//                                           SlimDirectedGraphView
//                                                   <Term> slim)
//    {
//
//        Integer[] order = new Integer[array.length];
//        for (int i = 0; i < order.length; i++) {
//            order[i] = i;
//        }
//        Arrays.sort(order, new Comparator<Integer>() {
//            @Override
//            public int compare(Integer o1, Integer o2) {
//                if (gte.getAnnotatedGenes(slim.getVertex(o1).getID()).totalAnnotated < res.getScore(o2)) {
//                    return 1;
//                }
//                if (res.getScore(o1) > res.getScore(o2)) {
//                    return -1;
//                }
//                return 0;
//            }
//        });
//
//    }

     public ArrayList<String> getTopDiseasesAsStrings(final ReducedBoqa.Result res,
                                                      ReducedBoqa rb) {

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
//"is now at"
        for (int i = 0; i < order.length; i++) {
            if (rb.allItemList.get(order[i]) .equals(trueDisease))
            {
                System.out.println("the actual disease is at position"
                + i); //(since the last disease is the right one)
                break;
            }

            }


        // Get top 20 results
        ArrayList<String> results = new ArrayList<String>();
        //for (int i = 0; i < 20 && i<order.length; i++) {

        //order.length/2\
         int THRESHOLD = 10;
        for (int i = 0; i < THRESHOLD ; i++) {
            int id = order[i]; //presumably, order[i] is now in order from lowest to msot score

            results.add(rb.allItemList.get(id).toString());
            //must use item2index here
//            results.add("item" + id); //bytestrings can be immediately constructed from this
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

//                QJK: why don;t we add noise to this?
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

    
    public double GetEntropy(ReducedBoqa.Result result)
    {
        
        double entropy = 0;
        for (double marg : result.marginals) {
            entropy += marg * marg;
        }
        return entropy;
    }
    
//    computes the entropies given an array of probabilities 
    public double GetEntropy(double [] probabilities)
    {
        return 0;
    }

//    QJK: Essentially, we WERE computing underneath the split before!
    public int ObsoleteQJKGetBestPhenotype(ReducedBoqa rb, double[] phenotype_frequencies) {
        //why not just maintain the observed array from before?
        //int [] top_phenotypes; //TODO keep the n largest phenotypes
        //we can use a min heap here! we want to maintain the largest k elements in a stream
        //(simply check the min elemnt, and if its larger, kick out one of the elements (replace)

        //we have the index, but not about the TermID

        int best_phenotype_index = -1;
        double best_phenotype_value = Double.POSITIVE_INFINITY;
        double temp = 0;
        ReducedBoqa.Result res;
        ArrayList<Integer> turnedOn;
        double val = 0;
        for (int i = 0; i < rb.o.observations.length; i++)//(int termInd : rb.o.observations)
        {
            if (!rb.o.observations[i]) {

                //actually, we must set ALL the above nodes to true too!
                //actually this is almost NEVER true (since negatives don't tell us anything more)

//                QJK: why don;t we add noise to this?

//                QJK: we should compute two universes: one where real_obs is True and one where real_obs is false
                rb.o.real_observations[i] = true;
                //we need this next line so that assignMarginals makes sense
                rb.o.observations[i] = true;
                turnedOn = setAncestors(rb, i);
//                System.out.println("starting assignmarginals");
                long start = System.nanoTime();

                res = rb.assignMarginals(rb.o, false, 1);
//                System.out.println("done assignmarginals. Took" + (System.nanoTime()-start));

                //after assignMarginals, the scores array is updated
                //boqa roll back: undoes the last one

                //if our current best worse than this new one, update it
//                System.out.println("starting scoring and misc");
                start = System.nanoTime();
                val = QJKscoringFunctionOnArray(res.marginals);

                //Should be memoized (in order to use the previous value (from the previous iteration), we need to
                //write it down somewhere!
                if (best_phenotype_value >=
                        (

                                temp = -1*

//                                    this term is large when the phenotype is likely and the entropy is HIGH
//                                    if we are picking specificially on this characteristic, then
//                                    we pick when this term is small since we get the negative, and then take it if it dcreases the value!
                                        (
                                                phenotype_frequencies[i] * (val)

//                                            This term will be large for those who are unlikely
//                                            + (1-phenotype_frequencies[i]) * val

                                        )



                                )) //pass in rb in case we need ot infer other things


                    {       //weight on the phenotype_frequency! (how likely it is to be true)

//                    QJK: process the result:
//                    We should compute the P(ph). We should compute the Entropy underneath this distribution
//                    function: calcualte DiseaseEntropy
//                    then, we should have a weighting . I guess, we already have the existing entropy, so it doesn't play an effect.
//                    actually, NO, we have a probability of the phenotype, and we NOT OBSERVED is NOT the same as not present and present! 

                    
//                    
                    best_phenotype_index = i;
                    best_phenotype_value = temp;

                    //index to termID is not well supported.
                    //since these things only make sense wrt a graph
                    //also termid and term are not well related, we really just want
                    //term,since it encapsulates everything
                    //hence we should use 3termcontainer to send it into

                }
//                System.out.println("done assignmarginals. Took" + (System.nanoTime()-start));

                //undoes the setting action
                rb.o.real_observations[i] = false;
                rb.o.observations[i] = false;
                rollback(rb, turnedOn);
            }

        }
        //rb.getOntology().getSlimGraphView().getVertex()//this can recover the Term if need be
        return best_phenotype_index;
    }


    public double[] normalizeMarginals(double [] freqs){

        double sum = 0;
        for (double marg : freqs) {
            sum += marg;
        }
        double [] normalized = Arrays.copyOf(freqs, freqs.length);
        for (int i =0; i <freqs.length; i++)
        {
            normalized[i] /= sum;
        }


        return normalized;

    }
    public double QJKscoringFunctionOnArray(double [] freqs) {

        double entropy = 0;

        for (double marg : normalizeMarginals(freqs)) {
//            double intermediate = Math.log(marg)*marg;
            entropy += Math.log(marg)*marg;
        }

        return entropy*-1;
    }

    public double scoringFunctionOnArray(double [] freqs) {

        return 1;
//        double entropy = 0;
//        for (double marg : normalizeMarginals(freqs)) {
//            double intermediate = Math.log(marg)*marg;
//            entropy += Math.log(marg)*marg;
//        }
//
//        return entropy*-1;
//        double score = 0;
//
//        for (double marg : freqs) {
//            score += marg * marg;
//        }
//
//        return -score;
    }

    public int QJKBestPhenotype(double [][] phenoDiseaseDist, ReducedBoqa rb)
    {
        int best_phenotype_index = -1;
        double best_phenotype_value = Double.POSITIVE_INFINITY;
        double temp = 0;
        double val = 0;
        for (int i = 0; i<phenoDiseaseDist.length; i++)
        {
            //assert rb.o.observations.length =phenoDiseaseDist.length
            //We cannot return pick a phenotype twice
            if (!rb.o.observations[i]) {
                if (phenotype_frequencies[i] != 0)

    val = QJKscoringFunctionOnArray(phenoDiseaseDist[i]);
//                    QJK: as it stands right now, we simply ask for the most likely phenotype!
                    if (best_phenotype_value >= (

                            temp = -1*

//                                    this term is large when the phenotype is likely and the entropy is HIGH
//                                    if we are picking specificially on this characteristic, then
//                                    we pick when this term is small since we get the negative, and then take it if it dcreases the value!
                                    (
                                            phenotype_frequencies[i] * (1.0/val)

//                                            This term will be large for those who are unlikely
//                                            + (1-phenotype_frequencies[i]) * val

                                    )
                                    )) {
                        best_phenotype_index = i;
                        best_phenotype_value = temp;
                    }
            }

        }
        return best_phenotype_index;
    }

    //pheno rows, disease cols
    public int multiGetBestPhenotype(double [][] phenoDiseaseDist, ReducedBoqa rb)
    {
        int best_phenotype_index = 0;
        double best_phenotype_value = Double.NEGATIVE_INFINITY;
        double temp = 0;
        double val;
        for (int i = 0; i<phenoDiseaseDist.length; i++)
        {
            //assert rb.o.observations.length =phenoDiseaseDist.length
            //We cannot return pick a phenotype twice
            if (!rb.o.observations[i]) {
                if (phenotype_frequencies[i] != 0)

//                    QJK: as it stands right now, we simply ask for the most likely phenotype!
                if (best_phenotype_value <= (

                        temp = phenotype_frequencies[i] * scoringFunctionOnArray(phenoDiseaseDist[i])

                )) {
                    best_phenotype_index = i;
                    best_phenotype_value = temp;
                }
            }

        }
        return best_phenotype_index;
    }
    public void computeVeniness(GOTermEnumerator gte, TermContainer tc,
    SlimDirectedGraphView<Term> slim)
    {


        for (Term t: tc)
        {
            if (gte.getAnnotatedGenes(t.getID()).totalAnnotatedCount() > 5000)
            System.out.println("This term " + t + "which has index" +
                            slim.getVertexIndex(t) + "has this many annotations" +

                    gte.getAnnotatedGenes(t.getID()).totalAnnotatedCount());
        }

    }
    //Represents the disease-phenotype frequency annotation data.
    //I1: Disease
    //I2: Phenotype
    //I3: Frequency Category
    double[] phi_phenotype_frequencies;
    double[] phenotype_frequencies;
    double[] disease_frequencies; //this is actually just BOQA's marginals
    HashMap<Term, HashMap<ByteString, Integer>> pheno_disease_freq;
    double [] freq_categories = {1, 0.9, 0.55, 0.175,0.02,0,0.30};
    ByteString trueDisease;
    Set<Term> trueDiseasePhentoypes; //perhaps an association container might have been best
    //AssociationContainer;
    Gene2Associations gx = new Gene2Associations(new ByteString("aa"));
    Gene2Associations trueDiseaseMapping;
    List<ByteString> all_diseases = new ArrayList<>();
    ByteString[] index2item;


    public int getIndexFromTermId()
    {
        return 101; //index of phenotypic abnormality
    }
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

            //go over the disease-frequency pairings
            //ideally, we have a bytestring->index
            //if we are just doing phenotype to disease, then we can directly use these elements
            for (Entry annotation : pheno_disease_freq.get(pheno_term).entrySet()) {

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

//        QJK normalization step:
        phi_phenotype_frequencies = normalizeMarginals(phi_phenotype_frequencies);
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

//QJK: normalization
        phenotype_frequencies = normalizeMarginals(phenotype_frequencies);


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

//    probability of a phenotype, is already known! we can do P(ph) * log p(ph)
//    we already summed up over all the diseases
//    find the average entropy information
//    Info Gain = 0.94- I(S)
//    note that our highest entropy is not necessarily the best we want! we want the difference
//    actually, maybe the lowest info IS the best as well

//    one is entropy based, the other can be info gain based


    public double scoringFunction(ReducedBoqa.Result result, ReducedBoqa rb) {
        double score = 0;

        for (double marg : result.marginals) {
            score += marg * marg;
        }

        return score;
    }

    public void QJKgenerateTrueDisease( AssociationContainer assocs)
    {

        int limit = all_diseases.size();

//        THIS MUST BE RANDOMIZED!!!
        Random rnd = new Random();
        int disease = rnd.nextInt(limit);
        trueDisease = all_diseases.get(disease);
        trueDiseaseMapping = assocs.get(trueDisease);



    }

    public void generateTrueDisease(SlimDirectedGraphView<Term> slim, AssociationContainer assocs)
    {
        Random rnd = new Random(); //this is our true disease
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
                pheno_disease_freq.get(t).put(trueDisease, rnd.nextInt(freq_categories.length)); //these correspond to the frequency classes

            }
            trueDiseaseMapping.add(trueDiseasePhenotype);


            assocs.addAssociation(trueDiseasePhenotype); //this seems to not hve any effect on BOQA... (nvm, it is used inb boqa.setup)
        }

        System.out.println(trueDiseaseMapping.getAssociations());


    }

    public void checkScoresAreEqual(ReducedConfiguration rc1,
                                  ReducedConfiguration rc2,
                                    ReducedBoqa rb)
    {



        System.out.println("The beta values are " + rb.getInitial_beta()
        + "for initial/naive beta" + rb.getExperimental_beta()
                );
        if (rc1.getScore(rb.getALPHA_GRID()[0],rb.getInitial_beta(),
                rb.getExperimental_beta()) ==
                rc2.getScore(rb.getALPHA_GRID()[0],rb.getInitial_beta(),
                rb.getExperimental_beta()))
        {
            System.out.println("they are the same score;");
        }

        else
        {
            System.out.println("they are not the same score");
        }

    }

    public void checkStatsAreSame(ReducedConfiguration rc1,
                                  ReducedConfiguration rc2)
    {
        for (int i = 0; i<rc1.stats.length;i++)
        {
            if (rc1.stats[i]!=rc2.stats[i])
            {
                System.out.println("the configs are not equal");
                return;
            }
        }

        System.out.println("The configs are equal;");
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

    public ReducedConfiguration printInfoAboutDisease(ReducedBoqa rb, int index)
    {
        System.out.println("This is the disease " +
                rb.allItemList.get(index));
        System.out.println("This is their annotations" + Arrays.toString(rb.items2Terms[index]));
        ReducedConfiguration rc = new ReducedConfiguration();

        int numTerms = rb.slimGraph.getNumberOfVertices();
        boolean [] hidden = new boolean[numTerms];
        for (int x: rb.items2Terms[index])
        {
            hidden[x] = true;
        }
        for (int i = 0; i < numTerms; i++) {
            ReducedConfiguration.NodeCase c = rb.getNodeCase(i, hidden, rb.o);
            rc.increment(c); //increment the case that c is in
        }

        System.out.println("this is its stats" + rc);

        return rc;
    }

//    public void getBatchNumber()
//    {
//        List<String> readSmallTextFile(String aFileName) throws IOException {
//        Path path = Paths.get(aFileName);
//        return Files.readAllLines(path, ENCODING);
//    } catch (IOException e) {
//
//            e.printStackTrace();
//
//        }
//    }

    @Test
    public void testConvergenceWrapper() throws IOException, OBOParserException, URISyntaxException
    {
        double sum = 0;
        int NUM_TESTS = 10;
//        try (BufferedWriter bw = new BufferedWriter(new FileWriter(FILENAME))) {
//
//            String content = "This is the content to write into file\n";
//
//            bw.write(content);
//
//            // no need to close it.
//            //bw.close();
//
//            System.out.println("Done");
//
//        } catch (IOException e) {
//
//            e.printStackTrace();
//
//        }

        for (int i = 0; i < NUM_TESTS; i++)
        {
            sum += testConvergence();
            System.out.println("Done test "+ i );
            System.out.printf( "ran %d tests, avg is %f\n", i, sum/i);
        }

        System.out.printf("DONE RUNNING\n" +
                "ran %d tests, avg is %f\n", NUM_TESTS, sum/NUM_TESTS);
        //calls testConvergence repeatedly and tracks how long it took
        //(idelaly, writes this info to a file)
        //and computes avg, std. dev and such
    }


    public ByteString getTrueDisease(AssociationContainer assocs)
    {
        Random r = new Random();
//        int i =0;
        int selected = r.nextInt(assocs.getAllAnnotatedGenes().size());
//        while ( i < selected)
//        {
//
//        }
        Set<ByteString> allDisease = assocs.getAllAnnotatedGenes();
        Iterator<ByteString> iter = allDisease.iterator();
        ByteString tD = null;
        for (int i = 0; i < selected; i++)
        {
            tD = iter.next();
            i++;

        }
        System.out.println("there are " + assocs.getAllAnnotatedGenes().size() + " diseases");
        return tD;
//        return assocs.get(tD);
//        return assocs.get(assocs.getAllAnnotatedGenes().)r.nextInt(assocs.getAllAnnotatedGenes().size());
//        assocs.getAllAnnotatedGenes().size()


    }

    //QJK: pretty sure this is how we get diseases now...
    public AssociationContainer getAnnotations(TermContainer tc)  throws OBOParserException, IOException, URISyntaxException
    {
        Map<Term, Integer> hpo2freq = getHPOToFreqMappingHardCoded(tc);
        initializeHashmap(tc);
        //need to add all the diseases in... the issue is that: we don't have a group by!
        //hence, we should check the diseaseContainer if it is already there


        AssociationContainer assocs = new AssociationContainer();

        System.out.println("Working Directory = " +
                System.getProperty("user.dir"));
        File inputFile = new File("C:\\Users\\johnp\\Desktop\\git_stuff\\boqa\\new_boqa\\resources\\phenotype_annotation.tab");
//        System.out.println(inputFile.toString());
        try {
            HpoDiseaseAnnotationParser parser = new HpoDiseaseAnnotationParser(inputFile);
            while (parser.hasNext()) {
                HpoDiseaseAnnotation anno = parser.next();

                if (anno.getDb().equals(ORPHA_name))
                {
//                    System.out.println(anno);
                    ByteString item = new ByteString(anno.getDbName());

//                Term t = slim.getVertex()
//                Term t  = new Term(new TermID(anno.getTermId().getIdWithPrefix()));

                    TermID t = new TermID(anno.getTermId().getIdWithPrefix());

                    //This is the solution!
                    Term tx = tc.get(t);



                    Association a = new Association(item, tx.getIDAsString());
                    all_diseases.add(item);
                    String freq_mod = anno.getFrequencyModifier();
                    if (freq_mod.equals(""))
                    {
                        pheno_disease_freq.get(tx).put(item, 6); //special 0.5 case
                    }

                    else{
//                        System.out.println(freq_mod);
//                        System.out.println(anno);
                        Term freq_term = tc.get(freq_mod);

                        pheno_disease_freq.get(tx).put(item, hpo2freq.get(freq_term));
                    }


//                pheno_disease_freq.get(t).put(item, hpo2freq.get(freq_term));



                    assocs.addAssociation(a);
                    //pick randomly using this
//                    assocs.getAllAnnotatedGenes().size()
                    //then: update the freq dict (so it goes from 0 to 6, 7 indices), and then:
                    //then: why can't we run BOQA?
                    //possibly some more tweaking with the false pos/neg on the diseases
                    //probably need to instantiate this in the other class, so that we can still access it

                }



            }
        } catch (IOException e) {
            System.err.println("Problem reading from file.");
        } catch (TermAnnotationParserException e) {
            System.err.println("Problem parsing file.");
        }
        System.out.println("done");
        return assocs;
    }


    public int QJKtestConvergence() throws IOException, OBOParserException, URISyntaxException {


        boolean noise = true;
        boolean give_free = false;
        int num = 10000;
        final ReducedBoqa boqa = new ReducedBoqa();
        //boqa.getOntology().
        //boqa.getOntology().getTerm() //FROM THE TERMID, we can recover the terms, and also recover the indexes?
        //yes, we are sure that the ints produced are the same ints as used in the Boqa.java (since, we get the
        //vertex2ancestors just immediately from boqa)

        //get the file, then get its canonical path
        OBOParser hpoParser = getOboParser();
        removeObsoleteTerms(hpoParser);
        AssociationContainer assocs;
        long start;


        //blackbox: it gets all the terms (in the HPO)
        //getTermMap returns a list of all terms!
        TermContainer tc = new TermContainer(hpoParser.getTermMap(), hpoParser.getFormatVersion(), hpoParser.getDate());
        Ontology ontology = new Ontology(tc);
        SlimDirectedGraphView<Term> slim = ontology.getSlimGraphView();


//        assocs = generateAnnotations(num, slim, tc);
        assocs = getAnnotations(tc);
//        trueDisease = new ByteString("item" + num);
//        trueDiseaseMapping = new Gene2Associations(trueDisease);

//        generateTrueDisease(slim, assocs);

        QJKgenerateTrueDisease(assocs);
        System.out.println("THIS IS THE TRUE DISEASE" + trueDisease);

        Observations o = new Observations();
        int numberOfTerms = ontology.getNumberOfTerms();
        o.observations = new boolean[numberOfTerms];
        o.real_observations = new boolean[numberOfTerms];

        //Run BOQA once to get the initial guesses.
        ArrayList<String> initial_guesses = null;

        boqa.setup(ontology, assocs);
        boqa.setO(o);
//        index2item = buildReverseArrayMapping(boqa.item2Index);

        int steps = 0;
        double increment = 0.00; //using no unobs neg!:
        boolean discovered = false;
        phenotype_frequencies = new double[numberOfTerms]; //alternatively, just copy over the
        //array length from the item2ancestors for example
        phi_phenotype_frequencies = new double[numberOfTerms];

        //initalization/first step stuff
        ReducedBoqa.Result res=new ReducedBoqa.Result();
        boqa.setInitial_beta(boqa.getInitial_beta()-increment * steps);
        res = boqa.assignMarginals(o, false, 1);
        disease_frequencies = res.marginals;
        long total = System.nanoTime();
        print_find_ancestors_of_trueDisease(boqa, tc);

        if (give_free)
        {
            int free = getFreeObs(boqa);
            //should set the free to true as well.
            setAncestors(boqa, free);
        }

//        o.observations[free] = true;
//        o.real_observations[free] = true;

        computeVeniness(boqa.termEnumerator,tc,slim);


//        IndexToTermPrinter.printMapping(slim.vertex2Index);
        while (!discovered) {
            ReducedBoqa.iteration++;



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

//            QJK: discounted for some reason:
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

//            call the multiGetBest one...
            int phenotype_to_check = QJKBestPhenotype(boqa.multiDiseaseDistributions,boqa); //in here we do all the phenotype checks

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

            boolean present_or_not = getObservation(index,boqa, noise);

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
//            TODO: semantics here are off: if it's IN the patient, then we should MISOBSERVE with a probability
//            TODO: let the getObservation handle all this: and have the case based misclassifications
            o.real_observations[index] = present_or_not;



            o.observations[index] = true; //recall the new meanings: observations means whether it was
            //checked, while the arraylist determines whether it was true or not
            //this has been deprecated

            //repeating this process should segregate everything
            //however, this can and SHOULD happen normally (probability of seeing something false
            //repeatedly is low. this is not a healing love, this is a wicked fantasy
            //test shoudl work:

            //updating steps for the next one:

            res = boqa.assignMarginals(o, false, 1);
            disease_frequencies = res.marginals;


//            System.out.println(java.util.Arrays.toString(res.marginals));
//            System.out.println(java.util.Arrays.toString(res.scores));

            int max_ind = printTopDisease(res, boqa);
            ReducedConfiguration rc1 = printInfoAboutDisease(boqa,max_ind);
            ReducedConfiguration rc2 = printInfoAboutDisease(boqa,boqa.item2Index.get(trueDisease));
            checkStatsAreSame(rc1,rc2);
            checkScoresAreEqual(rc1,rc2, boqa);
            //sorts the array, by getScore and takes the top N

            initial_guesses = getTopDiseasesAsStrings(res,boqa); //we have essentially the top ids now
            System.out.println(initial_guesses);
            //from the ids, we can get the mappings they have
//            System.out.println();
            //now, we recompute the marginals.
            //o.setValue()if ()
            if (trueDiseaseInTopNDiseases(trueDisease.toString(), initial_guesses)) {
                discovered = true;
                System.out.println("we are finished! took " + steps
                        + " guesses");
            }
            steps++;
            boqa.setInitial_beta(boqa.getInitial_beta()-increment);
            System.out.println("done loop iter. Took" + (System.nanoTime()-total));

        }
        return steps+1; //to stay consistent with the previous estimates
    }

    public int testConvergence() throws IOException, OBOParserException, URISyntaxException {
        boolean noise = true;
        boolean give_free = false;
        int num = 10000;
        final ReducedBoqa boqa = new ReducedBoqa();
        //boqa.getOntology().
        //boqa.getOntology().getTerm() //FROM THE TERMID, we can recover the terms, and also recover the indexes?
        //yes, we are sure that the ints produced are the same ints as used in the Boqa.java (since, we get the
        //vertex2ancestors just immediately from boqa)

        //get the file, then get its canonical path
        OBOParser hpoParser = getOboParser();
        removeObsoleteTerms(hpoParser);
        AssociationContainer assocs;
        long start;


        //blackbox: it gets all the terms (in the HPO)
        //getTermMap returns a list of all terms!
        TermContainer tc = new TermContainer(hpoParser.getTermMap(), hpoParser.getFormatVersion(), hpoParser.getDate());
        Ontology ontology = new Ontology(tc);
        SlimDirectedGraphView<Term> slim = ontology.getSlimGraphView();


//        assocs = generateAnnotations(num, slim, tc);
        assocs = getAnnotations(tc);
//        trueDisease = new ByteString("item" + num);
//        trueDiseaseMapping = new Gene2Associations(trueDisease);

//        generateTrueDisease(slim, assocs);
        QJKgenerateTrueDisease(assocs);
//        trueDisease = getTrueDisease(assocs);
//        trueDiseaseMapping = assocs.get(trueDisease);


        Observations o = new Observations();
        int numberOfTerms = ontology.getNumberOfTerms();
        o.observations = new boolean[numberOfTerms];
        o.real_observations = new boolean[numberOfTerms];

        //Run BOQA once to get the initial guesses.
        ArrayList<String> initial_guesses = null;

        boqa.setup(ontology, assocs);
        boqa.setO(o);
//        index2item = buildReverseArrayMapping(boqa.item2Index);

        int steps = 0;
        double increment = 0.00; //using no unobs neg!:
        boolean discovered = false;
        phenotype_frequencies = new double[numberOfTerms]; //alternatively, just copy over the
        //array length from the item2ancestors for example
        phi_phenotype_frequencies = new double[numberOfTerms];


        //initalization/first step stuff
        ReducedBoqa.Result res=new ReducedBoqa.Result();
        boqa.setInitial_beta(boqa.getInitial_beta()-increment * steps);
        res = boqa.assignMarginals(o, false, 1);
        disease_frequencies = res.marginals;
        long total = System.nanoTime();
        print_find_ancestors_of_trueDisease(boqa, tc);

        if (give_free)
        {
            int free = getFreeObs(boqa);
            //should set the free to true as well.
            setAncestors(boqa, free);
        }

//        o.observations[free] = true;
//        o.real_observations[free] = true;

        computeVeniness(boqa.termEnumerator,tc,slim);

//        IndexToTermPrinter.printMapping(slim.vertex2Index);
        while (!discovered) {
            ReducedBoqa.iteration++;



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

//            QJK: discounted for some reason:
//            boqa.assignMultiShotMarginals(o,false,8);

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

//            call the multiGetBest one...

//            QJK; note that multiDiseaseDistributions is NOT computed! Hence, we are doing passing in misleading stuff here!!

            int phenotype_to_check = ObsoleteQJKGetBestPhenotype(boqa, phenotype_frequencies); //in here we do all the phenotype checks

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

            boolean present_or_not = getObservation(index,boqa, noise);

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

            res = boqa.assignMarginals(o, false, 1);
            disease_frequencies = res.marginals;


//            System.out.println(java.util.Arrays.toString(res.marginals));
//            System.out.println(java.util.Arrays.toString(res.scores));

            int max_ind = printTopDisease(res, boqa);
            ReducedConfiguration rc1 = printInfoAboutDisease(boqa,max_ind);
            ReducedConfiguration rc2 = printInfoAboutDisease(boqa,boqa.item2Index.get(trueDisease));
            checkStatsAreSame(rc1,rc2);
            checkScoresAreEqual(rc1,rc2, boqa);
            //sorts the array, by getScore and takes the top N

            initial_guesses = getTopDiseasesAsStrings(res,boqa); //we have essentially the top ids now
            System.out.println(initial_guesses);
            //from the ids, we can get the mappings they have
//            System.out.println();
            //now, we recompute the marginals.
            //o.setValue()if ()
            if (trueDiseaseInTopNDiseases(trueDisease.toString(), initial_guesses)) {
                discovered = true;
                System.out.println("we are finished! took " + steps
                + " guesses");
            }
            steps++;
            boqa.setInitial_beta(boqa.getInitial_beta()-increment);
            System.out.println("done loop iter. Took" + (System.nanoTime()-total));

        }
        return steps+1; //to stay consistent with the previous estimates
    }

    public static OBOParser getOboParser() throws URISyntaxException, IOException, OBOParserException {
        AssociationContainer assocs;


        URL resource = ClassLoader.getSystemResource("hp.obo");
        ClassLoader cl = ClassLoader.getSystemClassLoader();

//        URL[] urls = ((URLClassLoader)cl).getURLs();
//
//        for(URL url: urls){
//            System.out.println(url.getFile());
//        }
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

        return hpoParser;
    }

    private int printTopDisease(ReducedBoqa.Result res,
                                 ReducedBoqa rb) {
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
        System.out.println("this corresponds to " + rb.allItemList.get(max_ind));
        return max_ind;
    }

    //we could also provide an array implementation too
    public Map<Integer, ByteString> buildReverseMapping(
            Map<ByteString,Integer> item2index)
    {
        Map<Integer, ByteString> hm = new HashMap<>();
        for (Entry <ByteString,Integer> e : item2index.entrySet())
        {
            hm.put(e.getValue(),e.getKey());

        }

        return hm;

    }

    public ByteString[] buildReverseArrayMapping(
            Map<ByteString,Integer> item2index)
    {
        ByteString[] bs = new ByteString[item2index.size()];
        for (Entry <ByteString,Integer> e : item2index.entrySet())
        {
            bs[e.getValue()] = e.getKey();

        }
        //assert no elements are null
        for (int i =0; i < bs.length; i++)
        {
            assert bs[i]!=null;
        }
        return bs;

    }

    private void removeObsoleteTerms(OBOParser hpoParser) {
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
    }
}
