package sonumina.boqa.calculation;

import ontologizer.association.AssociationContainer;
import ontologizer.enumeration.GOTermEnumerator;
import ontologizer.enumeration.ItemEnumerator;
import ontologizer.go.Ontology;
import ontologizer.go.Term;
import ontologizer.go.TermID;
import ontologizer.set.PopulationSet;
import ontologizer.types.ByteString;
import sonumina.math.graph.SlimDirectedGraphView;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

/**
 * Created by johnchen on 18/05/17.
 */
public class ReducedBoqa {




    //Most likely will be an array of all the frequencies for phenotype
    //Note that: this array should always sum to 1. (should it?)-- ask brudno about this
    public Object PhenotypeFrequencyDistributions;

    //todo: note that this already exists--it is the marginals
    //an array of all the frequencies for disease. This can be re-computed at each stage.
    public Object DiseaseFrequencyDistributions;

    //However, it is a good question: is it the marginals being recomputed each time, or just them getting updated?
    //this is more of a implementation/engineering issue
    //TODO: use a lambda to pass in arbitrary function at call
    public double scoringFunction(double [] phen_arr, double [] dis_arr )
    {
        double sum = 0;
        //for example, just compute sum of squares here.
        for (double elt: dis_arr)
        {
            sum+=elt;
        }

        return sum;
//        for (int i = 0; i < phen_arr.length; i++)
//        {
//
//        }

    }

    public void scoringFunctionWrapper(double [] phen_arr, double [] dis_arr )
    {

        //Calls scoringFunction each time
        //Picking one on or off would cause
        for (int i = 0; i < phen_arr.length; i++)
        {

        }
    }

    //really i want a dictionary
    //DIctioanry of observations: Term-Index (vertex) to Boolean
    //Another choice could have been Term-ID to Boolean
    //public HashMap<Integer,Boolean> registered_obserations; //we still need this to store whether observed
    public Observations o;

    private double initial_beta = 0.3 ; //this will approach the experimental beta
    private double experimental_beta = 0.01 ;

    public void setInitial_beta(double new_beta)
    {
        if (new_beta >= experimental_beta)
        {
            initial_beta = new_beta;
        }

    }

    public double getInitial_beta()
    {
        return initial_beta;
    }


    Ontology graph;
    AssociationContainer assoc;
    public ArrayList<ByteString> allItemList;

    /** Contains all the ancestors of the terms */
    public int[][] term2Ancestors;
    //these are all JAGGED ARRAYS

    /** Contains the parents of the terms */
    public int[][] term2Parents;

    /** Contains the children of the term */
    public int[][] term2Children;

    /** Contains the descendants of the (i.e., children, grand-children, etc.) */
    public int[][] term2Descendants;

    /** Map items to their index */
    public HashMap<ByteString, Integer> item2Index;

    public GOTermEnumerator termEnumerator;

    public SlimDirectedGraphView<Term> slimGraph;


    private long timeDuration;

    public int[][] items2Terms;

    /**
     * For each item, contains the term ids which need to be switched on, if the previous item was on.
     * this is misleading: i believe this is saying which terms need to be on if the disease item is on
     * No, actually it is saying exactly what it says: an order-dependent way of quickly turing things on
     */
    public int[][] diffOnTerms;

    /**
     * Same as diffOnTerms but for switching off terms.
     */
    public int[][] diffOffTerms;

    /**
     * Similar to diffOnTerms but each adjacent frequency-implied state
     */
    public int[][][] diffOnTermsFreqs;

    /**
     * Similar to diffOffTerms but each adjacent frequency-implied state
     */
    public int[][][] diffOffTermsFreqs;



    private double ALPHA_GRID[] = new double[] {  0.001 };

    private double BETA_GRID[] = new double[] { 0.01 };

    public void setup(Ontology ontology, AssociationContainer assocs)
    {
        //this.registered_obserations = new HashMap<>();
        //this.o = new Observations();

        this.assoc = assocs;
        this.graph = ontology;

        this.term2Parents = this.graph.getSlimGraphView().vertexParents; //recall this was an array of arrays
        this.term2Children = this.graph.getSlimGraphView().vertexChildren;
        this.term2Ancestors = this.graph.getSlimGraphView().vertexAncestors;
        this.term2Descendants = this.graph.getSlimGraphView().vertexDescendants; //presumably, this is a vertex induced subgraph

        this.slimGraph = graph.getSlimGraphView();

        HashSet<ByteString> allItemsToBeConsidered = new HashSet<ByteString>(this.assoc.getAllAnnotatedGenes());


        PopulationSet allItems = new PopulationSet("all");
        allItems.addGenes(allItemsToBeConsidered);
        this.termEnumerator =  //since a.getEvidence is always false
                allItems.enumerateGOTerms(this.graph, this.assoc, null);
        //a
        ItemEnumerator itemEnumerator = ItemEnumerator.createFromTermEnumerator(this.termEnumerator);



        this.allItemList = new ArrayList<ByteString>();
        this.item2Index = new HashMap<ByteString, Integer>(); //these are the diseases
        int i = 0;
        for (ByteString item : itemEnumerator) { //we need this to get the actual terms!!
            this.allItemList.add(item);
            this.item2Index.put(item, i);
            i++;
        }



        this.items2Terms = new int[this.allItemList.size()][];
        i = 0;
        for (ByteString item : itemEnumerator) {
            int j = 0;
            //these all the terms connected via the induced subgraph (i.e. following the annotation propagation rule)
            ArrayList<TermID> tids = itemEnumerator.getTermsAnnotatedToTheItem(item);
            this.items2Terms[i] = new int[tids.size()]; //jagged arrays

            for (TermID tid : tids) {
                this.items2Terms[i][j++] = this.slimGraph.getVertexIndex(this.graph.getTerm(tid));
            }

            Arrays.sort(this.items2Terms[i]);
            i++;
        }

        createDiffVectors(); //need this crucial line from provideGlobals()
    }

    boolean areFalsePositivesPropagated()
    {

        return false;
    }

    boolean areFalseNegativesPropagated()
    {

        return false;
    }

    private ReducedConfiguration.NodeCase getNodeCase(int node, boolean[] hidden,Observations o)
    {
        if (areFalsePositivesPropagated()) {
            /* Here, we consider that false positives are inherited.
             * J: hence if my child is on, then no questions asked, i am on
              * (that is, when I am also on, I consider the case to be because of INHERIT_TRUE
              * and that correspondingly sets my calculations)*/
            for (int i = 0; i < this.term2Children[node].length; i++) {
                int chld = this.term2Children[node][i]; //TODO use a lambda here, or some python syn
                if (o.real_observations.get(chld)) { //this explains false positive propagation
                    //translation: if both are observed has no real analogue anymore
                    //actually, if we actually do do true propagation, then
                    //we can just check that both are in
                    //however, this implies that we should have "truly observed"
                    //since, it could just be on solely because one of its kids was accidentally on


                    //if both are unobserved, it DOES have an analogue though! Both are not in
                    //the registered observations

                    if (o.real_observations.get(node)) { //in the observed layer, parent points to child, so
                        //if child is on, then parenbt is on [contrary to how BN semantrics normally
                        //work, where really only parent affects child]
                        return ReducedConfiguration.NodeCase.INHERIT_TRUE; //the causality seems backwards
                        //it (the parent) is only on because its child was als oon
                        //however, this could have been handled immediately in a different graph
                        //preprocessing step

                        //we assume that it became true via inheritance
                    } else {
                        /* NaN */
                        return ReducedConfiguration.NodeCase.FAULT;
                    }
                }
            }
        }

        if (areFalseNegativesPropagated()) {
            /* Here, we consider that false negatives are inherited */
            for (int i = 0; i < this.term2Parents[node].length; i++) {
                int parent = this.term2Parents[node][i];
                //handle both false and null case
                //TODO fix this so that null is checked FIRST (to avoid !null errors)
                if (!o.real_observations.get(parent) || o.real_observations.get(parent) == null) {
                    if (!o.real_observations.get(node) || o.real_observations.get(node) == null) {
                        return ReducedConfiguration.NodeCase.INHERIT_FALSE; //wasn't actually false but inherited it..
                    } else {
                        /* NaN */
                         return ReducedConfiguration.NodeCase.FAULT;
                    }
                }
            }
        }

        if (hidden[node]) {
            /* Term is truly on */ //only 4 possible cases;
            // this only takes into account the hidden and the query variable
            //unlike the paper which takes into account the parents
            //hidden and observed are both on?? even though they are different nodes
            //i mean doesnt each node in the BN get its own id??
            //1-t-1 correspondence, I'll allow it... (since these are "slices" of the BN, unoike
            //previously where we were iterating over an array of size ALL of the ITEMS)

            //when they are NOT inherited, what happens? a node CANNOT be made into a false
            //by another node--hence it IS conclusively false positve/negative
            //observed just means whether it was actually looked at or not
            if (o.real_observations.containsKey(node) && o.real_observations.get(node)) {
                return ReducedConfiguration.NodeCase.TRUE_OBSERVED_POSITIVE;
                //ensure this is based on the correct beta being used (when being scored)
                //however, it depends. Say we go with annotation propagation. Then, there is a
                //chance that it could be on due to its children. Then we can reward still checking
                //by being able to take advanatage of the smaller beta (i.e. a True true positive)
                //OTOH we can also consider that a node only ever turns on if it was observed. Hence,
                //we can automatically use the small beta
                //the issue with the "lookup beta" even for the positives approach is that in general,
                //we do not condone checking a parent of a node you already checked.
                //however, under our current semantics, the probabilities no longer add to 1!
                //since (1-e_beta) + (n_beta) > 1
                //this probably has an effect on the probability distributions though
            } else {

                //really, this should just be false
                if (o.real_observations.containsKey(node))
                {
                    return ReducedConfiguration.NodeCase.FALSE_OBSERVED_NEGATIVE;

                }
                return ReducedConfiguration.NodeCase.FALSE_UNOBSERVED_NEGATIVE; //this implicitly had two cases before
                //hidden[node] AND the thing is NOT in the realo observations
                //note that one case is messed up, since in order for it to have a value of true, it MUST be in the real_observations,
                //i.e. AN observation must have been made

            }
        } else {
            //Note this will cause errors!: !o.real_observations.get(node)
            /* Term is truly off */
            if (o.real_observations.get(node) == null|| !o.real_observations.get(node)) {
                return ReducedConfiguration.NodeCase.TRUE_NEGATIVE;
            } else {
                //however, this doesn't really make sense anymore: how can it be fp if indeed
                //it was not in hidden yet also
                return ReducedConfiguration.NodeCase.FALSE_POSITIVE;
            }
        }
    }
    //either need an instance variable or another way of maintaining the actual obs
    private void determineCases(Observations o, boolean[] hidden, ReducedConfiguration stats)
    {
        //total number of nodes
        int numTerms = this.slimGraph.getNumberOfVertices();

        for (int i = 0; i < numTerms; i++) {
            ReducedConfiguration.NodeCase c = getNodeCase(i, hidden, o);
            stats.increment(c); //increment the case that c is in
        }
    }


    static public class Result
    {
        /** Contains the marginal probability for each item */
        public double[] marginals;

        /** Contains the marginal probability for each item */
        public double[] marginalsIdeal;
        //hence we can lookup how likely a disease is (and sort it later)


        public double[] scores;



        /** Some statistics for each item (number of false-positives, etc. ) */
        ReducedConfiguration[] stats;

        /**
         * Get the score of the given item.
         *
         * @param i
         * @return
         */
        public double getScore(int i)
        {
            return this.scores[i];
        }

        /**
         * Get the marginal probability of the given item.
         *
         * @param i
         * @return
         */
        public double getMarginal(int i)
        {
            return this.marginals[i];
        }

        public double getMarginalIdeal(int i)
        {
            return this.marginalsIdeal[i];
        }

        public ReducedConfiguration getStats(int i)
        {
            return this.stats[i];
        }

        public int size()
        {
            return this.marginals.length;
        }
    }
    public void setO(Observations o)
    {
        this.o = o;
    }
    public Ontology getOntology()
    {
        return this.graph;
    }

    private static int[] setDiff(int[] a, int[] b)
    {
        int[] c = new int[a.length];
        int cc = 0; /* current c */

        /* Obviously, this could be optimized to linear time if a and b would be assumed to be sorted */
        for (int element : a) {
            boolean inB = false;

            for (int element2 : b) {
                if (element == element2) {
                    inB = true;
                    break;
                }
            }

            if (!inB) {
                c[cc++] = element;
            }
        }
        int[] nc = new int[cc];
        for (int i = 0; i < cc; i++) {
            nc[i] = c[i];
        }
        return nc;
    }

    private void createDiffVectors()
    {
        int i;

        long sum = 0;
        /* Fill diff matrix */
        this.diffOnTerms = new int[this.allItemList.size()][];
        this.diffOffTerms = new int[this.allItemList.size()][];
        this.diffOnTerms[0] = this.items2Terms[0]; /* For the first step, all terms must be activated */
        //this makes sense, since the one "behind" it is the empty set; hence all must be taken!


        //all terms annotated to the first item are diffOnTerms for that item as well

        //
        this.diffOffTerms[0] = new int[0];   //this is an empty array
        for (i = 1; i < this.allItemList.size(); i++) {
            int prevOnTerms[] = this.items2Terms[i - 1];      //items2terms[0] on first iteration
            int newOnTerms[] = this.items2Terms[i];

            this.diffOnTerms[i] = setDiff(newOnTerms, prevOnTerms);
            this.diffOffTerms[i] = setDiff(prevOnTerms, newOnTerms);

            sum += this.diffOnTerms[i].length + this.diffOffTerms[i].length;
        }

    }



    private WeightedConfigurationList determineCasesForItem(int item, Observations o,
                                                            boolean takeFrequenciesIntoAccount, boolean[] previousHidden, ReducedConfiguration previousStats)
    {
        int numTerms = this.slimGraph.getNumberOfVertices();

        if (previousHidden == null && previousStats != null) {
            throw new IllegalArgumentException();
        }
        if (previousHidden != null && previousStats == null) {
            throw new IllegalArgumentException();
        } //presumably they just cannot be equal to one another (or maybe only these specific
        //cases are forbidden)

        long now = System.nanoTime();

        /* Tracks the hidden state ReducedConfiguration that matches the observed state best */
        // double bestScore = Double.NEGATIVE_INFINITY;
        // boolean [] bestTaken = new boolean[numTermsWithExplicitFrequencies];

        WeightedConfigurationList statsList = new WeightedConfigurationList();

        boolean[] hidden;
        ReducedConfiguration stats;

        if (previousHidden == null) {
            hidden = new boolean[numTerms];
        } else {
            hidden = previousHidden;
        }

        if (previousStats == null) {
            stats = new ReducedConfiguration();
        } else {
            stats = previousStats;
        }
        //we don't need diffOn, since we have the actual items as well as ancestors etc. we can just
        //compute everything on the fly


        if (!takeFrequenciesIntoAccount) {
            /* New */
            int[] diffOn = this.diffOnTerms[item];
            int[] diffOff = this.diffOffTerms[item];

            /* Decrement config stats of the nodes we are going to change */
            //J: we will change them back (update) later!
            //presumably this is JUST like items2terms, except it takes into account some existing frequency stuff
            ReducedConfiguration decrement_config = new ReducedConfiguration();
            ReducedConfiguration increment_config = new ReducedConfiguration();
            for (int element : diffOn) {
                decrement_config.increment(getNodeCase(element, hidden, o));
                stats.decrement(getNodeCase(element, hidden, o));
                //these 3 uniquely can identify a state
            }
            for (int element : diffOff) {
                decrement_config.increment(getNodeCase(element, hidden, o));
                stats.decrement(getNodeCase(element, hidden, o)); //lookup the hidden and observed too
            }

            /* Change nodes states */ //why is hidden[0] always on, esp. when we set observed to be on only
            for (int i = 0; i < diffOn.length; i++) {
                hidden[diffOn[i]] = true; //each element in the diffOn array, will change
                //(regarde, diffon is fixed for the item)
            }
            for (int i = 0; i < diffOff.length; i++) {
                hidden[diffOff[i]] = false; //what we need to switch off to get to disease i
            }

            /* Increment config states of nodes that we have just changed */
            for (int element : diffOn) {
                increment_config.increment(getNodeCase(element, hidden, o));
                stats.increment(getNodeCase(element, hidden, o));
            }
            for (int element : diffOff) {
                increment_config.increment(getNodeCase(element, hidden, o));
                stats.increment(getNodeCase(element, hidden, o));
            } //this winds up exactly undoing what we just earlier did. however, apparently, the state of hidden
            //must have changed.

            statsList.add(stats.clone(), 0);
        }

        this.timeDuration += System.nanoTime() - now;

        return statsList;
    }

    public Result assignMarginals(final Observations observations, final boolean takeFrequenciesIntoAccount,
                                       final int numThreads)
    {

        int i;

        final Result res = new Result();
        res.scores = new double[this.allItemList.size()]; //its an array size of the entire BN
        res.marginals = new double[this.allItemList.size()];
        res.marginalsIdeal = new double[this.allItemList.size()];
        res.stats = new ReducedConfiguration[this.allItemList.size()]; // an ARRAY of ReducedConfigurations..
        //perhaps they are saying how for each node in THE BN, there is a cofig?


        //this is just initializing each stats[i] to contain something
        for (i = 0; i < res.stats.length; i++) {
            res.stats[i] = new ReducedConfiguration(); //each ITEM gets its own ReducedConfiguration
        }
        for (i = 0; i < res.scores.length; i++) {
            res.scores[i] = Double.NEGATIVE_INFINITY;
        }

        //J for all of this
        //n^3 matrix, one dimension each for parameters (alpha and beta) and one dimension for all the scores
        //i.e. at this alpha, beta, what are all the scores for all the nodes?
        final double[][][] scores = new double[this.allItemList.size()][this.ALPHA_GRID.length][this.BETA_GRID.length];
        final double[] idealScores = new double[this.allItemList.size()]; //this is likely a "target" whom we wish to approximate
        //as much as possible
        //interestingly we also have ideal marginals, which I am not sure how the two differ by



        final ExecutorService es;
        if (numThreads > 1) {
            es = Executors.newFixedThreadPool(numThreads);
        } else {
            es = null;
        }

        final boolean[] previousHidden = new boolean[this.slimGraph.getNumberOfVertices()];
        final ReducedConfiguration previousStat = new ReducedConfiguration();

        //ReducedConfiguration is not null, but: it also was not initalized to anything //previousstat contains all the info from runnign determinecases--do we run it first just for the multithreading?; part of the diffOn, etc.?
        determineCases(observations, previousHidden, previousStat);
        //this pops back with all the cases; and having incremented the particular stats
        //(i.e previousStat)

        ArrayList<Future<?>> futureList = new ArrayList<Future<?>>();

        for (i = 0; i < this.allItemList.size(); i++) {
            final int item = i;

            Runnable run = new Runnable()
            {

                @Override
                public void run() //since this is an inner, class, we need final int item
                        //we cannot change it inside this scope!
                {

                    //for 1 item, we get a WCL
                    //then we check it against all values of alpha, beta
                    //then we ge tthe score of the WCL
                    //there is a list of ReducedConfigurations, and we get
                    //
                    //this call does the differntials, and also returns the list to work on
                    //(with all the nodes set); but it does NOT determine the cases!
                    //i believe the stats does not get updated
                    WeightedConfigurationList stats =
                            determineCasesForItem(item, observations, takeFrequenciesIntoAccount,
                                    numThreads > 1 ? null : previousHidden, numThreads > 1 ? null : previousStat);

                    //stats only has 1 element
                    //moreover, we can get a new element each time using the iterator
                    //System.out.print("aaa");
                    for (WeightedConfiguration wc : stats) //this is a list of weighted configs, which itself is simply a list of
                    {
                        System.out.println("Itme " + item);
                        System.out.println(wc.stat);



                    }
                    //since multithreading with the differentials would be too difficult,
                    //we only sequentially set flags in the stats, (when we only have one process)

                    //PERHAPS: we should assign stats to some of the res.
                    //In general, res.stats is not updated again, sadly

                    System.out.println(stats.toString());
                    //if numthreads > 1 then null, else previousHidden

                    //J: the scoring function is critical; also this seems to be finding the best
                    //rate for alpha and beta.at a particular alpha, beta value
                    //this puts the score in for an item;

                    //int temp_terms [] = new int[] ; // since this is inside a void run inside an interface
                    //might be tricky...

                    for (int a = 0; a < ReducedBoqa.this.ALPHA_GRID.length; a++) {
                        for (int b = 0; b < ReducedBoqa.this.BETA_GRID.length; b++) {
                            //This scor is a weighted sum!
                            scores[item][a][b] = stats.score(ReducedBoqa.this.ALPHA_GRID[a],
                                    experimental_beta,initial_beta);

                            //normally we maximize this
                            res.scores[item] = Util.logAdd(res.scores[item], scores[item][a][b]);
                        }
                    }

                    System.out.println(observations.observationStats != null);
                    if (observations.observationStats != null) {
                        double fpr = observations.observationStats.falsePositiveRate();
                        if (fpr == 0) {
                            fpr = 0.0000001; //get the false positive rate, which for some reason
                            //cannot be exactly 0 or 1 (this is so that we do not get a 0 term in the
                            //product which will just wipe our probabilities

                        } else if (fpr == 1.0) {
                            fpr = 0.999999;
                        } else if (Double.isNaN(fpr)) {
                            fpr = 0.5;
                        }

                        double fnr = observations.observationStats.falseNegativeRate();
                        if (fnr == 0) {
                            fnr = 0.0000001;
                        } else if (fnr == 1) {
                            fnr = 0.999999;
                        } else if (Double.isNaN(fnr)) {
                            fnr = 0.5;
                        }

                        idealScores[item] = stats.score(fpr, experimental_beta,initial_beta);
                        //run it again with the **actual** scores
                    }
                }
            };

            if (es != null) {
                futureList.add(es.submit(run));
            } else {
                run.run();
            }
        }

        if (es != null) {
            es.shutdown();

            for (Future<?> f : futureList) {
                try {
                    f.get();
                } catch (Exception e) {
                    throw new RuntimeException(e);
                }
            }

            try {
                while (!es.awaitTermination(10, TimeUnit.SECONDS)) {
                    ;
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        double normalization = Math.log(0);
        double idealNormalization = Math.log(0);
        //already at this point, why are all the scores the same?
        for (i = 0; i < this.allItemList.size(); i++) {
            normalization = Util.logAdd(normalization, res.scores[i]);
            idealNormalization = Util.logAdd(idealNormalization, idealScores[i]);
        }

        for (i = 0; i < this.allItemList.size(); i++) {
            res.marginals[i] = Math.min(Math.exp(res.scores[i] - normalization), 1); //no, also the scores are same
            res.marginalsIdeal[i] = Math.min(Math.exp(idealScores[i] - idealNormalization), 1);

            // System.out.println(i + ": " + idealScores[i] + " (" + res.getMarginalIdeal(i) + ") " + res.scores[i] +
            // " (" + res.getMarginal(i) + ")");
            // System.out.println(res.marginals[i] + " " + res.marginalsIdeal[i]);
        }

        if (res.marginalsIdeal[observations.item] < res.marginals[observations.item]) {
            for (i = 0; i < this.allItemList.size(); i++) {
                res.marginalsIdeal[i] = res.marginals[i];
            }
        }

        // System.out.println(idealNormalization + " " + normalization);
        // if (exitNow)
        // System.exit(10);
        return res;
    }
}
