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

    Ontology graph;
    AssociationContainer assoc;
    public ArrayList<ByteString> allItemList;

    /** Contains all the ancestors of the terms */
    public int[][] term2Ancestors;

    /** Contains the parents of the terms */
    public int[][] term2Parents;

    /** Contains the children of the term */
    public int[][] term2Children;

    /** Contains the descendants of the (i.e., children, grand-children, etc.) */
    public int[][] term2Descendants;

    /** Map items to their index */
    public HashMap<ByteString, Integer> item2Index;

    GOTermEnumerator termEnumerator;

    SlimDirectedGraphView<Term> slimGraph;


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
        this.assoc = assocs;
        this.graph = ontology;
        //without graph inducing--how long will the computation take?
        //why is inducement a valid thing?
        //is it because of the only one valid configuration for a disease type thing ? => 0's in the Hidden
        //still mess it up since they could be 1 in the observed (False Positive). Hence
        //not entirely sure of rationale (if rational) to break it up
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
        this.item2Index = new HashMap<ByteString, Integer>();
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
            this.items2Terms[i] = new int[tids.size()];

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

        return true;
    }

    boolean areFalseNegativesPropagated()
    {

        return true;
    }

    private Configuration.NodeCase getNodeCase(int node, boolean[] hidden, boolean[] observed)
    {
        if (areFalsePositivesPropagated()) {
            /* Here, we consider that false positives are inherited.
             * J: hence if my child is on, then no questions asked, i am on
              * (that is, when I am also on, I consider the case to be because of INHERIT_TRUE
              * and that correspondingly sets my calculations)*/
            for (int i = 0; i < this.term2Children[node].length; i++) {
                int chld = this.term2Children[node][i]; //TODO use a lambda here, or some python syn
                if (observed[chld]) {
                    if (observed[node]) { //in the observed layer, parent points to child, so
                        //if child is on, then parenbt is on [contrary to how BN semantrics normally
                        //work, where really only parent affects child]
                        return Configuration.NodeCase.INHERIT_TRUE; //the causality seems backwards
                        //it (the parent) is only on because its child was als oon
                        //however, this could have been handled immediately in a different graph
                        //preprocessing step

                        //we assume that it became true via inheritance
                    } else {
                        /* NaN */
                        return Configuration.NodeCase.FAULT;
                    }
                }
            }
        }

        if (areFalseNegativesPropagated()) {
            /* Here, we consider that false negatives are inherited */
            for (int i = 0; i < this.term2Parents[node].length; i++) {
                int parent = this.term2Parents[node][i];
                if (!observed[parent]) {
                    if (!observed[node]) {
                        return Configuration.NodeCase.INHERIT_FALSE; //wasn't actually false but inherited it..
                    } else {
                        /* NaN */
                         return Configuration.NodeCase.FAULT;
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
            if (observed[node]) {
                return Configuration.NodeCase.TRUE_POSITIVE;
            } else {
                return Configuration.NodeCase.FALSE_NEGATIVE;
            }
        } else {
            /* Term is truly off */
            if (!observed[node]) {
                return Configuration.NodeCase.TRUE_NEGATIVE;
            } else {
                return Configuration.NodeCase.FALSE_POSITIVE;
            }
        }
    }

    private void determineCases(boolean[] observedTerms, boolean[] hidden, Configuration stats)
    {
        //total number of nodes
        int numTerms = this.slimGraph.getNumberOfVertices();

        for (int i = 0; i < numTerms; i++) {
            Configuration.NodeCase c = getNodeCase(i, hidden, observedTerms);
            stats.increment(c); //increment the case that c is in
        }
    }


    static public class Result
    {
        /** Contains the marginal probability for each item */
        private double[] marginals;

        /** Contains the marginal probability for each item */
        private double[] marginalsIdeal;
        //hence we can lookup how likely a disease is (and sort it later)


        private double[] scores;

        /** Some statistics for each item (number of false-positives, etc. ) */
        Configuration[] stats;

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

        public Configuration getStats(int i)
        {
            return this.stats[i];
        }

        public int size()
        {
            return this.marginals.length;
        }
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



    private WeightedConfigurationList determineCasesForItem(int item, boolean[] observed,
                                                            boolean takeFrequenciesIntoAccount, boolean[] previousHidden, Configuration previousStats)
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

        /* Tracks the hidden state configuration that matches the observed state best */
        // double bestScore = Double.NEGATIVE_INFINITY;
        // boolean [] bestTaken = new boolean[numTermsWithExplicitFrequencies];

        WeightedConfigurationList statsList = new WeightedConfigurationList();

        boolean[] hidden;
        Configuration stats;

        if (previousHidden == null) {
            hidden = new boolean[numTerms];
        } else {
            hidden = previousHidden;
        }

        if (previousStats == null) {
            stats = new Configuration();
        } else {
            stats = previousStats;
        }
        //we don't need diffOn, since we have the actual items as well as ancestors etc. we can just
        //compute everything on the fly

        //
        takeFrequenciesIntoAccount = false;
        if (!takeFrequenciesIntoAccount) {
            /* New */
            int[] diffOn = this.diffOnTerms[item];
            int[] diffOff = this.diffOffTerms[item];

            /* Decrement config stats of the nodes we are going to change */
            //J: we will change them back (update) later!
            //presumably this is JUST like items2terms, except it takes into account some existing frequency stuff

            for (int element : diffOn) {
                stats.decrement(getNodeCase(element, hidden, observed));
                //these 3 uniquely can identify a state
            }
            for (int element : diffOff) {
                stats.decrement(getNodeCase(element, hidden, observed));
            }

            /* Change nodes states */
            for (int i = 0; i < diffOn.length; i++) {
                hidden[diffOn[i]] = true;
            }
            for (int i = 0; i < diffOff.length; i++) {
                hidden[diffOff[i]] = false;
            }

            /* Increment config states of nodes that we have just changed */
            for (int element : diffOn) {
                stats.increment(getNodeCase(element, hidden, observed));
            }
            for (int element : diffOff) {
                stats.increment(getNodeCase(element, hidden, observed));
            }

            statsList.add(stats.clone(), 0);
        }

        this.timeDuration += System.nanoTime() - now;

        return statsList;
    }

    public Result assignMarginals(final Observations observations, final boolean takeFrequenciesIntoAccount,
                                       final int numThreads)
    {
//        System.out.println("woot");
//        System.out.println("woot");
//        System.out.println("woot");
//        System.out.println("woot");
//        System.out.println("woot");

        int i;

        final Result res = new Result();
        res.scores = new double[this.allItemList.size()]; //its an array size of the entire BN
        res.marginals = new double[this.allItemList.size()];
        res.marginalsIdeal = new double[this.allItemList.size()];
        res.stats = new Configuration[this.allItemList.size()]; // an ARRAY of configurations..
        //perhaps they are saying how for each node in THE BN, there is a cofig?


        //this is just initializing each stats[i] to contain something
        for (i = 0; i < res.stats.length; i++) {
            res.stats[i] = new Configuration(); //each ITEM gets its own configuration
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

        //Score each item:

        for (i = 0; i < this.allItemList.size(); i++) {
            final int item = i;

            WeightedConfigurationList stats =
                    determineCasesForItem(item, observations.observations, takeFrequenciesIntoAccount,
                            null, null); //run it single threaded for now

        }

        final ExecutorService es;
        if (numThreads > 1) {
            es = Executors.newFixedThreadPool(numThreads);
        } else {
            es = null;
        }

        final boolean[] previousHidden = new boolean[this.slimGraph.getNumberOfVertices()];
        final Configuration previousStat = new Configuration();

        //Configuration is not null, but: it also was not initalized to anything
        determineCases(observations.observations, previousHidden, previousStat);
        //this pops back with all the cases; and having incremented the particular stats
        //(i.e previousStat)

        ArrayList<Future<?>> futureList = new ArrayList<Future<?>>();

        for (i = 0; i < this.allItemList.size(); i++) {
            final int item = i;

            Runnable run = new Runnable()
            {

                @Override
                public void run()
                {
                    //for 1 item, we get a WCL
                    //then we check it against all values of alpha, beta
                    //then we ge tthe score of the WCL
                    //there is a list of configurations, and we get
                    //
                    WeightedConfigurationList stats =
                            determineCasesForItem(item, observations.observations, takeFrequenciesIntoAccount,
                                    numThreads > 1 ? null : previousHidden, numThreads > 1 ? null : previousStat);
                    //if numthreads > 1 then null, else previousHidden

                    //J: the scoring function is critical; also this seems to be finding the best
                    //rate for alpha and beta.at a particular alpha, beta value
                    //this puts the score in for an item;
                    for (int a = 0; a < ReducedBoqa.this.ALPHA_GRID.length; a++) {
                        for (int b = 0; b < ReducedBoqa.this.BETA_GRID.length; b++) {
                            //This scor is a weighted sum!
                            scores[item][a][b] = stats.score(ReducedBoqa.this.ALPHA_GRID[a],
                                    ReducedBoqa.this.BETA_GRID[b]);
                            res.scores[item] = Util.logAdd(res.scores[item], scores[item][a][b]);
                        }
                    }


                    if (observations.observationStats != null) {
                        double fpr = observations.observationStats.falsePositiveRate();
                        if (fpr == 0) {
                            fpr = 0.0000001;
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

                        idealScores[item] = stats.score(fpr, fnr);
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

        for (i = 0; i < this.allItemList.size(); i++) {
            normalization = Util.logAdd(normalization, res.scores[i]);
            idealNormalization = Util.logAdd(idealNormalization, idealScores[i]);
        }

        for (i = 0; i < this.allItemList.size(); i++) {
            res.marginals[i] = Math.min(Math.exp(res.scores[i] - normalization), 1);
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
