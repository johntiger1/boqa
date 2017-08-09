/* Copyright (c) 2010-2012 Sebastian Bauer
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted (subject to the limitations in the
 * disclaimer below) provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the
 *   distribution.
 *
 * * Neither the name of Sebastian Bauer nor the names of its
 *   contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
 * GRANTED BY THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT
 * HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package sonumina.boqa.calculation;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ontologizer.association.Association;
import ontologizer.association.AssociationContainer;
import ontologizer.association.Gene2Associations;
import ontologizer.dotwriter.AbstractDotAttributesProvider;
import ontologizer.dotwriter.GODOTWriter;
import ontologizer.enumeration.GOTermEnumerator;
import ontologizer.enumeration.ItemEnumerator;
import ontologizer.go.Ontology;
import ontologizer.go.Term;
import ontologizer.go.TermID;
import ontologizer.set.PopulationSet;
import ontologizer.types.ByteString;
import sonumina.algorithms.Algorithms;
import sonumina.math.distribution.ApproximatedEmpiricalDistribution;
import sonumina.math.graph.SlimDirectedGraphView;

/**
 * This is core class implementing BOQA. Currently, it also implements other procedures based on semantic similarity
 * measures but this is planned to be refactored. In order to perform the algorithm, you need to setup a new BOQA object
 * a ontology and associations via setup(). One can then use assignMarginals() on the observations to obtain marginal
 * probabilities. Refer to the BOQATest class for a working example usage.
 *
 * @author Sebastian Bauer
 * @see setup, sonumina.boqa.tests.BOQATest
 */

//Default constructor
public class BOQA
{

    public boolean [] observations; //records if term i has been observed
    public boolean [] RegisteredObservaitons; // records whether the observation was true or false



    public BOQA()
    {

        for (int i = 1; i < 19; i++)
        {
            System.out.println("RO-BOQA"); //reducued online boqa

        }
    }

    //public int [];



    /** Our logger */
    private static Logger logger = LoggerFactory.getLogger(BOQA.class);

    private Ontology graph;

    private AssociationContainer assoc;

    /** Term enumerator */
    private GOTermEnumerator termEnumerator;

    /** Slim variant of the graph */
    //private SlimDirectedGraphView<Term> slimGraph;
    public SlimDirectedGraphView<Term> slimGraph;

    /** An array of all items */
    public ArrayList<ByteString> allItemList;

    /** Map items to their index */
    public HashMap<ByteString, Integer> item2Index;

    /** Links items to terms */
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

    /**
     * The factors of each combination.
     */
    public double[][] factors;

    /** Links items to directly associated terms */
    public int[][] items2DirectTerms;

    /**
     * Links items to the frequencies of corresponding directly associated terms. Frequencies are interpreted as
     * probabilities that the corresponding term is on.
     */
    public double[][] items2TermFrequencies; //really this should be an array of tuples (combine it with items2DirectTerms)


    /**
     * J: sorted freqeuncies (TODO look into making the array optimally sized, for example max # of annotations)
     * This contains the (ascending) order of the items2TermFrequencies, E.g., use item2TermFrequenciesOrder[0][2] to
     * determine the term that is associated to first item and has the third lowest frequency.
     */
    public int[][] item2TermFrequenciesOrder;

    /** Indicates whether an item have explicit frequencies */
    public boolean[] itemHasFrequencies;

    /** Contains all the ancestors of the terms */
    public int[][] term2Ancestors;

    /** Contains the parents of the terms */
    public int[][] term2Parents;

    /** Contains the children of the term */
    public int[][] term2Children;

    /** Contains the descendants of the (i.e., children, grand-children, etc.) */
    public int[][] term2Descendants;

    /** Contains the order of the terms */
    public int[] termsInTopologicalOrder;

    /** Contains the topological rank of the term */
    public int[] termsToplogicalRank;

    /** Contains the IC of the terms */
    public double[] terms2IC;

    /** Contains the term with maximum common ancestor of two terms */
    private int micaMatrix[][];

    /** Contains the jaccard index */
    private double jaccardMatrix[][];

    /** Contains the query cache, needs to be synched when accessed */
    private QuerySets queryCache;

    /** Used to parse frequency information */
    public static Pattern frequencyPattern = Pattern.compile("(\\d+)\\.?(\\d*)\\s*%");

    public static Pattern frequencyFractionPattern = Pattern.compile("(\\d+)/(\\d+)");

    /* Settings for generation of random data */
    // private final double ALPHA = 0.002; // 0.01
    private double ALPHA = 0.002;

    private double BETA = 0.10; // 0.1

    /* Settings for inference */
    private double ALPHA_GRID[] = new double[] { 1e-10, 0.0005, 0.001, 0.005, 0.01 };

    private double BETA_GRID[] = new double[] { 1e-10, 0.005, 0.01, 0.05, 0.1, 0.2, 0.4, 0.8, 0.9 };

    // private static boolean CONSIDER_FREQUENCIES_ONLY = false;
    // private final static String [] evidenceCodes = null;
    // private final static int SIZE_OF_SCORE_DISTRIBUTION = 250000;
    // public static int maxTerms = -1;

    private boolean CONSIDER_FREQUENCIES_ONLY = true;

    private final String[] evidenceCodes = null;// new String[]{"PCS","ICE"};

    private int SIZE_OF_SCORE_DISTRIBUTION = 250000;

    private final int NUMBER_OF_BINS_IN_APPROXIMATED_SCORE_DISTRIBUTION = 10000;

    private int maxTerms = -1; /* Defines the maximal number of terms a query can have */

    private int maxFrequencyTerms = 10; /* Maximal number of frequency terms (k in the paper) */

    /** False positives can be explained via inheritance */
    private static int VARIANT_INHERITANCE_POSITIVES = 1 << 0; //left shift by 0 = same

    /** False negatives can be explained via inheritance */
    private static int VARIANT_INHERITANCE_NEGATIVES = 1 << 1; // left shift 1 by 1 bit == 0?!!?
    //false, you are thinking of a right shift. A left shit will make it 2

    /** Model respects frequencies */
    private static int VARIANT_RESPECT_FREQUENCIES = 1 << 2;

    /** Defines the model as a combination of above flags */
    private int MODEL_VARIANT = VARIANT_RESPECT_FREQUENCIES | VARIANT_INHERITANCE_NEGATIVES;// |

    // VARIANT_INHERITANCE_POSITIVES;

    /** If set to true, empty observation are allowed */
    private boolean ALLOW_EMPTY_OBSERVATIONS = false;

    /** Precalculate the jaccard matrix */
    private boolean PRECALCULATE_JACCARD = false;

    /** Use cached MaxIC terms. Speeds up Resnik */
    private boolean PRECALCULATE_MAXICS = true;

    /** Use precalculated max items. Speeds up Resnik */
    private boolean PRECALCULATE_ITEM_MAXS = true;

    /** Cache the queries */
    private final boolean CACHE_RANDOM_QUERIES = true;

    /** Forbid illegal queries */
    private final boolean FORBID_ILLEGAL_QUERIES = true;

    /** Cache the score distribution during calculation */
    private boolean CACHE_SCORE_DISTRIBUTION = true;

    /** Precalculate score distribution. Always implies CACHE_SCORE_DISTRIBUTION. */
    private boolean PRECALCULATE_SCORE_DISTRIBUTION = true;

    /** Tries to load the score distribution */
    private boolean TRY_LOADING_SCORE_DISTRIBUTION = true;

    /** Identifies whether score distribution should be stored */
    private boolean STORE_SCORE_DISTRIBUTION = true;

    /** Defines the maximal query size for the cached distribution */
    private int MAX_QUERY_SIZE_FOR_CACHED_DISTRIBUTION = 20;

    /* Some configuration stuff */

    /**
     * Sets the number of terms that can be selected to be on during the simulation.
     *
     * @param maxTerms
     */
    public void setSimulationMaxTerms(int maxTerms)
    {
        this.maxTerms = maxTerms;
    }

    /**
     * Returns the number of terms that can be selected to be on during the simulation.
     *
     * @return
     */
    public int getSimulationMaxTerms()
    {
        return this.maxTerms;
    }

    /**
     * Set alpha value used for generateObservations() and used for the ideal FABN scoring.
     *
     * @param alpha
     */
    public void setSimulationAlpha(double alpha)
    {
        this.ALPHA = alpha;
    }

    /**
     * Returns the simulation alpha.
     *
     * @return
     */
    public double getSimulationAlpha()
    {
        return this.ALPHA;
    }

    /**
     * Set beta value used for generateObservations() and used for the ideal FABN scoring.
     *
     * @param alpha //todo
     */
    public void setSimulationBeta(double beta)
    {
        this.BETA = beta;
    }

    /**
     * Returns the simulation beta.
     *
     * @return
     */
    public double getSimulationBeta()
    {
        return this.BETA;
    }

    /**
     * Sets, whether only items with frequencies should be considered.
     *
     * @param frequencies
     */
    public void setConsiderFrequenciesOnly(boolean frequencies)
    {
        this.CONSIDER_FREQUENCIES_ONLY = frequencies;
    }

    /**
     * Returns, whether only items with frequencies should be considered.
     *
     * @return
     */
    public boolean getConsiderFrequenciesOnly()
    {
        return this.CONSIDER_FREQUENCIES_ONLY;
    }

    /**
     * Sets the maximum query size of for a cached distribution.
     *
     * @param size
     */
    public void setMaxQuerySizeForCachedDistribution(int size)
    {
        this.MAX_QUERY_SIZE_FOR_CACHED_DISTRIBUTION = size;
    }

    /**
     * Precalculate score distribution.
     *
     * @param precalc
     */
    public void setPrecalculateScoreDistribution(boolean precalc)
    {
        this.PRECALCULATE_SCORE_DISTRIBUTION = precalc;
    }

    /**
     * Sets whether jaccard similarity shall be precalculated.
     *
     * @param precalc
     */
    public void setPrecalculateJaccard(boolean precalc)
    {
        this.PRECALCULATE_JACCARD = precalc;
    }

    /**
     * Sets, whether maxICs should be precalculated.
     *
     * @param precalc
     */
    public void setPrecalculateMaxICs(boolean precalc)
    {
        this.PRECALCULATE_MAXICS = precalc;
    }

    /**
     * Set whether we cache the score distribution.
     *
     * @param cache
     */
    public void setCacheScoreDistribution(boolean cache)
    {
        this.CACHE_SCORE_DISTRIBUTION = cache;
    }

    /**
     * Set whether we store the score distribution.
     *
     * @param store
     */
    public void setStoreScoreDistriubtion(boolean store)
    {
        this.STORE_SCORE_DISTRIBUTION = store;
    }

    /**
     * Sets whether the matrix that contains the max ic term of two given terms shall be precalculated.
     *
     * @param precalc
     */
    public void setPrecalculateItemMaxs(boolean precalc)
    {
        this.PRECALCULATE_ITEM_MAXS = precalc;
    }

    /**
     * Sets whether score distribution should be loaded.
     *
     * @param loading
     * @return
     */
    public void setTryLoadingScoreDistribution(boolean loading)
    {
        this.TRY_LOADING_SCORE_DISTRIBUTION = loading;
    }

    /**
     * Sets the size of the score distribution.
     *
     * @param size
     */
    public void setSizeOfScoreDistribution(int size)
    {
        this.SIZE_OF_SCORE_DISTRIBUTION = size;
    }

    /**
     * Returns the size of the score distribution.
     *
     * @return
     */
    public int getSizeOfScoreDistribution()
    {
        return this.SIZE_OF_SCORE_DISTRIBUTION;
    }

    /**
     * Returns the number of terms considered in for frequency analysis.
     *
     * @return
     */
    public int getMaxFrequencyTerms()
    {
        return this.maxFrequencyTerms;
    }

    /**
     * Sets the number of terms considered in for frequency analysis.
     *
     * @param newMaxFrequencyTerms
     */
    public void setMaxFrequencyTerms(int newMaxFrequencyTerms)
    {
        this.maxFrequencyTerms = newMaxFrequencyTerms;
    }

    /**
     * Returns whether false negatives are propagated in a top-down fashion.
     *
     * @return
     */
    public boolean areFalseNegativesPropagated()
    {
        boolean FNProp = (this.MODEL_VARIANT & VARIANT_INHERITANCE_POSITIVES) != 0;
        return (this.MODEL_VARIANT & VARIANT_INHERITANCE_POSITIVES) != 0;
    }

    /**
     * Returns whether false positives are propagated in a bottom-up fashion.
     *
     * @return
     */
    public boolean areFalsePositivesPropagated()
    {
        boolean FPProp = (this.MODEL_VARIANT & VARIANT_INHERITANCE_NEGATIVES) != 0;
        return (this.MODEL_VARIANT & VARIANT_INHERITANCE_NEGATIVES) != 0; //for all intents and purposes,
        //these are constants (never changed and can be reduced to flags)
    }

    /**
     * Returns whether all false stuff is propagated.
     *
     * @return
     */
    public boolean allFalsesArePropagated()
    {
        return areFalseNegativesPropagated() && areFalsePositivesPropagated();
    }

    /**
     * Returns whether frequencies should be respected.
     *
     * @return
     */
    public boolean respectFrequencies()
    {
        return (this.MODEL_VARIANT & VARIANT_RESPECT_FREQUENCIES) != 0;
    }

    /**
     * Returns the case for the given node, given the hidden and observed states.
     *
     * @param node
     * @param hidden
     * @param observed
     * @return
     */
    //We know the propagation is at play here, so a True kid implies a true Actual

    //Anno prop only occurs in the Hidden and Observed layers
    //this seems to treat it like a TERM (by using the term2Children lookup table
    //STEP THE BACK: this only makes sense (i.e. NodeCase is ONLY defined for) observed nodes
    //in fact, observe: term2children will be quite large, even when only 1 disease
    private ReducedConfiguration.NodeCase getNodeCase(int node, boolean[] hidden, boolean[] observed)
    {

        //these two do false ___ propagation in the observed layer
        //the last one does no propagation (since recall propagation remains sequestered within
        //the layer)
        //also note: if we EVER either observe child in FPpro or NOT observe parent in FNpro,
        //then we are GUARANTEED to reach an exit case. Negate the statements, and we find:
        //
        //if we dont observe parent in FPpro => and-parents is 0, and HENCE, it is the case from the paper..?
        //yet this one does not take into account the hidden layer..

        //there seems to be some weird conflation of states (1;0) and probability (alpha, beta etc.)
        //in this code. It makes sense that If i false positive a node in the observed layer, this will propagate
        //upwards through the observed layer.

        //However, this doesn't really MEAN anything in terms of probability!
        //or does it? if i was a false positive, does this say anything about my parents or the hidden node?
        //perhaps we fall into the and-parent case only because of the annotation propagation rule
        //ALSO NOTE: we do NOT consider the 0,1 cases...

        if (areFalsePositivesPropagated()) {
            /* Here, we consider that false positives are inherited.
             * J: hence if my child is on, then no questions asked, i am on
              * (that is, when I am also on, I consider the case to be because of INHERIT_TRUE
              * and that correspondingly sets my calculations)*/
            for (int i = 0; i < this.term2Children[node].length; i++) {
                int chld = this.term2Children[node][i]; //TODO use a lambda here, or some python syn
                if (observed[chld]) { //theere is no and-parents here!
                    //in fact, we are looking at the children, and certainly not ALL THE children!

                    if (observed[node]) { //in the observed layer, parent points to child, so
                        //if child is on, then parenbt is on [contrary to how BN semantrics normally
                        //work, where really only parent affects child]
                        return ReducedConfiguration.NodeCase.INHERIT_TRUE; //the causality seems backwards

                        //this is explaining WHY the node is on--it is on because of ONE of its kids!
                        //it (the parent) is only on because its child was als oon
                        //however, this could have been handled immediately in a different graph
                        //preprocessing step

                        //we assume that it became true via inheritance
                    } else {
                        /* NaN */
                        logger
                            .error(
                                "A child of a node is on although the parent is not: Impossible configuration encountered!");
                        return ReducedConfiguration.NodeCase.FAULT;
                    }
                }

                //note that nothing happens if the child was not observed
            }
        }
        //Essentially, BOTH boil down to propagation of child states to parent states. NOTE: at this
        //point we do not necessarily have the probability assignments...
        //this also goes back to the "explaining away" discussed in the textbooks etc.
        if (areFalseNegativesPropagated()) { //note that this detemrinsitic propagation occurs in BOTH the hidden
            //and observed layers

            /* Here, we consider that false negatives are inherited */
            for (int i = 0; i < this.term2Parents[node].length; i++) {
                int parent = this.term2Parents[node][i];
                if (!observed[parent]) {
                    if (!observed[node]) { //i.e. that we get our falses from a parent!!
                        return ReducedConfiguration.NodeCase.INHERIT_FALSE;
                    } else {
                        /* NaN */
                        logger
                            .error(
                                "A parent of a node is off although the child is not: Impossible configuration encountered!");
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
            //1-t-1 correspondence, I'll allow it... The index could/should
            // work for both of them
            //also: here, since FN and FP propagation checks either failed or completed...
            //however, we need only show that one part of this doesn't work:
            //let's say NEITHER FP or FN are propagated. Then, this should still work:
            //might be some terrible latch logic where they just left it out even though
            //not intuive at all

            //Ok, so neither FN/FP are propagated. In this case,

            //hidden & observed,
            // (since these are "slices" of the BN, unoike
            //previously where we were iterating over an array of size ALL of the ITEMS)

            //when they are NOT inherited, what happens? a node CANNOT be made into a false
            //by another node--hence it IS conclusively false positve/negative

            //consider JUST the and-parents. they will ALWAYS be on, because their children are on?
            //essentially, we are arguing that the and of a node's parents will not be 0
            // when there is no propagation
            //but how would it have arisen "regularly"? i.e. how can parent's nodes be 0 if a child is on? since
            //it must annotate all the things anyways.
            //if the child is off, then the paper says that the parents will be off with 1 certainty
            //(i am not sure about the Implication, but it simple states that P(Q_i = 0| H_I = 1 and Q_parents = 0) = 0
            //so actually the inference goes the OTHER way (we can conclude q_i state given the parents), not the other
            //way around
            //so they are saying that parents cannot be 0 EVEN if the child is 0!!
            //the observed cases make sense, since naturally, they will cause their parents to be anded
            //however, the (hidden, Not observed) and (not hidden, not observed) require explanation
            //i.e. the false and true negatives. There is no guarantee that false negatives and true negatives
            //that parents are kept on

            //these correspond to the parents-and = 1 in the paper
            if (observed[node]) { //this synchronizes with the results from the paper (assuming above)
                return ReducedConfiguration.NodeCase.TRUE_POSITIVE;
            } else {
                return ReducedConfiguration.NodeCase.FALSE_NEGATIVE; //this is P=1 in the textbook
            }


        } else { // if hidden is off
            /* Term is truly off */
            if (!observed[node]) {
                return ReducedConfiguration.NodeCase.TRUE_NEGATIVE;
            } else {
                return ReducedConfiguration.NodeCase.FALSE_POSITIVE; //alpha (from the paper)
            }
        }
    }

    /**
     * Determines the cases of the observed states given the hidden states. Accumulates them in states.
     *
     * @param observedTerms
     * @param hidden
     * @param stats
     */
    private void determineCases(boolean[] observedTerms, boolean[] hidden, Configuration stats)
    {
        //total number of nodes
        int numTerms = this.slimGraph.getNumberOfVertices();

        for (int i = 0; i < numTerms; i++) {
            Configuration.NodeCase c = getNodeCase(i, hidden, observedTerms);
            stats.increment(c); //increment the case that c is in
        }
    }

    private long timeDuration;

    /**
     * Determines the case of the given items and the given observations.
     *
     * @param item the disease (not the term!)
     * @param observed
     * @param takeFrequenciesIntoAccount select, if frequencies should be taken into account.
     * @param hiddenStorage is the storage used to store the hidden states. It must correspond to the states of the
     *            previous item (item -1). If this is the first item, it must be 0.
     * @return
     */
    private WeightedConfigurationList determineCasesForItem(int item, boolean[] observed,
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

        /* Tracks the hidden state configuration that matches the observed state best */
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

        if (!takeFrequenciesIntoAccount) {
            /* New */
            int[] diffOn = this.diffOnTerms[item];
            int[] diffOff = this.diffOffTerms[item];

            /* Decrement config stats of the nodes we are going to change */
            //J: we will change them back (update) later!
            //presumably this is JUST like items2terms, except it takes into account some existing frequency stuff

            for (int element : diffOn) { //on the 1st iteration: wont this literally decrement everything tho?
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

        } else {
            /* Initialize stats */
            if (previousHidden != null) {
                for (int i = 0; i < hidden.length; i++) {
                    hidden[i] = false;
                }
            }
            stats.clear();
            determineCases(observed, hidden, stats); // not recursive!

            /*
             * Loop over all tracked configurations that may appear due to the given item being active
             * J: this is the complext, O(2^k) case mentioned in the paper
             */
            for (int c = 0; c < this.diffOnTermsFreqs[item].length; c++) {
                int[] diffOn = this.diffOnTermsFreqs[item][c];
                int[] diffOff = this.diffOffTermsFreqs[item][c];

                /* Decrement config stats of the nodes we are going to change */
                for (int element : diffOn) {
                    stats.decrement(getNodeCase(element, hidden, observed));
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

                /* Determine cases and store */
                statsList.add(stats.clone(), this.factors[item][c]);
            }
        }

        this.timeDuration += System.nanoTime() - now;
        logger.debug(this.timeDuration / (1000 * 1000) + " " + statsList.size());

        return statsList;
    }

    /**
     * Returns the log probability that the given term has the observed state given the hidden states. If one of its
     * more specific terms (descendants in this case) are on then the probability that the observed term is on is one.
     * Otherwise the probability depends on the false-positive/false-negative rate. J: note that right
     * now it only seems to deal with the alpha/beta case (false). It is saying that only those that take alpha
     * and beta into account in the observed layer will call this. Also, it actually also has inherit false
     * and true
     *
     * @param termIndex
     * @param alpha
     * @param beta
     * @param hidden
     * @param observed
     * @return
     */
    public double scoreNode(int termIndex, double alpha, double beta, boolean[] hidden, boolean[] observed)
    {

        double score = 0.0;

        Configuration.NodeCase c = getNodeCase(termIndex, hidden, observed);

        switch (c) {
            case FALSE_NEGATIVE:
                score = Math.log(beta);
                break;
            case FALSE_POSITIVE:
                score = Math.log(alpha);
                break;
            case TRUE_POSITIVE:
                score = Math.log(1 - beta);
                break;
            case TRUE_NEGATIVE:
                score = Math.log(1 - alpha);
                break;
            case INHERIT_FALSE:
                score = Math.log(1); //this is 0 (i.e. we assume it was fine)/did not contribute
                //to the product; normally, it COULD otherwise be a false neg/pos
                //note: the inheritors "are fully explained" by the other nodes

                break;
            case INHERIT_TRUE:
                score = Math.log(1);
                break;
            default:
                break;
        }
        return score;
    }

    /**
     * Score a hidden configuration given the observations.
     *
     * @param observedTerms
     * @param stats
     * @param score
     * @param hidden
     * @return
     */
    @SuppressWarnings("unused")
    private double scoreHidden(boolean[] observedTerms, double alpha, double beta, boolean[] hidden)
    {
        Configuration stats = new Configuration();
        determineCases(observedTerms, hidden, stats);
        double newScore = stats.getScore(alpha, beta);
        return newScore;
    }

    /**
     * Calculates the score, when the given item is activated.
     *
     * @param item which is supposed to be active.
     * @param observedTerms
     * @param stats, some statistics about false positives etc. @TODO update by removing stats, which got moved inside
     * @param takeFrequenciesIntoAccount
     * @return
     */
    public double score(int item, double alpha, double beta, boolean[] observedTerms,
        boolean takeFrequenciesIntoAccount)
    {
        WeightedConfigurationList stats =
            determineCasesForItem(item, observedTerms, takeFrequenciesIntoAccount, null, null);
        return stats.score(alpha, beta);
    }

    /**
     * Returns the result of a logical or operation of the parents state.
     *
     * @param v
     * @param states
     * @return
     */
    public boolean orParents(int v, boolean[] states)
    {
        int[] parents = this.term2Parents[v];
        for (int parent : parents) {
            if (states[parent]) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns the result of a logical and operation of the parents state.
     *
     * @param v
     * @param states
     * @return
     */
    public boolean andParents(int v, boolean[] states)
    {
        for (int i = 0; i<10000; i++) {
            System.out.println("i am here though");
        }
        int[] parents = this.term2Parents[v];
        for (int i = 0; i < parents.length; i++) {
            if (!states[parents[i]]) {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns the result of a logical and operation of the children state.
     *
     * @param v
     * @param states
     * @return
     */
    public boolean andChildren(int v, boolean[] states)
    {
        int[] children = this.term2Children[v];
        for (int i = 0; i < children.length; i++) {
            if (!states[children[i]]) {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns the result of a logical or operation of the children state.
     *
     * @param v
     * @param states
     * @return
     */
    public boolean orChildren(int v, boolean[] states)
    {
        int[] children = this.term2Children[v];
        for (int element : children) {
            if (states[element]) {
                return true;
            }
        }
        return false;
    }

    /**
     * Generates observation according to the model parameter for the given __ITEM___item.
     * J: this is most likely the TESTING code?
     * @param item
     * @return
     */
    public Observations generateObservations(int item, Random rnd)
    {
        int retry = 0;

        Observations o = null;

        do {
            int i;
            int[] falsePositive = new int[this.slimGraph.getNumberOfVertices()];
            int numFalsePositive = 0;
            int[] falseNegative = new int[this.slimGraph.getNumberOfVertices()];
            int numFalseNegative = 0;
            int numMissedInHidden = 0;

            int numPositive = 0;
            int numHidden = 0;

            boolean[] observations = new boolean[this.slimGraph.getNumberOfVertices()];
            boolean[] hidden = new boolean[this.slimGraph.getNumberOfVertices()];

            boolean CONSIDER_ONLY_DIRECT_ASSOCIATIONS = true;

            if (CONSIDER_ONLY_DIRECT_ASSOCIATIONS) { //Note: item = disease
                logger.debug("Item {} has {} annotations", item, this.items2DirectTerms[item].length);
                for (i = 0; i < this.items2DirectTerms[item].length; i++) {
                    //our "i" is actually all the children of items (all the terms that the disease
                    //annotates

                    boolean state = true;

                    //potentially it might not be on (state)
                    if (respectFrequencies()) {
                        state = rnd.nextDouble() < this.items2TermFrequencies[item][i];

                        logger.debug(this.items2DirectTerms[item][i] + "("
                            + this.items2TermFrequencies[item][i] + ")="
                            + state);
                    }

                    if (state) {
                        //this is our classic ancestor propagation (we just set it to
                        hidden[this.items2DirectTerms[item][i]] = state;
                        observations[this.items2DirectTerms[item][i]] = state;

                        activateAncestors(this.items2DirectTerms[item][i], hidden);
                        activateAncestors(this.items2DirectTerms[item][i], observations);

                        numPositive++;
                    } else {
                        numMissedInHidden++;
                    }
                }

            } else {
                //THIS actually seems to be only conbsidering direct assoc...
                //(since activateAncestors will activate ALL ANCESTORS,
                //note ancestors contains ALL ancestors (not just immediate parents)
                //TODO speed up here
                for (i = 0; i < this.items2Terms[item].length; i++) {
                    hidden[this.items2Terms[item][i]] = true;
                    observations[this.items2Terms[item][i]] = true;
                    numPositive++;
                }
            }

            /* Fill in false and true positives */
            //randomly set some nodes to be counted as FP/FN
            for (i = 0; i < observations.length; i++) {
                double r = rnd.nextDouble();
                if (observations[i]) {
                    if (r < this.BETA) {
                        falseNegative[numFalseNegative++] = i;
                        // System.out.println("false negative " + i);
                    }
                } else { //set on to off and vicea versa
                    if (r < this.ALPHA) {
                        falsePositive[numFalsePositive++] = i;
                        // System.out.println("false positive " + i);
                    }
                }
            }

            /* apply false negatives */
            if (areFalseNegativesPropagated()) {
                /*
                 * false negative, but also make all descendants negative. They are considered as inherited in this case
                 */
                for (i = 0; i < numFalseNegative; i++) {
                    observations[falseNegative[i]] = false; //set that paritcular node to false
                    //note that it doesn't do anuything about the PARENTS
                    deactivateDecendants(falseNegative[i], observations);
                    //they got turned off, since the and-parents will be negative.,
                }
            } else {
                /* false negative */
                for (i = 0; i < numFalseNegative; i++) {
                    observations[falseNegative[i]] = false;
                }

                /* fix for true path rule */
                //this has the effect of undoing some of the false negatives
                //(i.e. setting observations[falsenegative[i] back to TRUE)
                //however, it must be for the true path rule (i can follow any path back up)
                //to be true
                //hence only leaf nodes can become false negatives

                //Note annotation propagation and True path rule are synonyms
                //also: TODO there can be a speed up, since probably quite a bit of
                //repeated work (i.e. I already turned it on, no need to go through the entire
                // tree repeatedly)
                for (i = 0; i < observations.length; i++) {
                    if (observations[i]) {
                        activateAncestors(i, observations);
                    }
                }

                //this fixes the case where I am on but my parent is not
            }

            /* apply false positives */
            if (areFalsePositivesPropagated()) {
                /* fix for true path rule */
                for (i = 0; i < numFalsePositive; i++) {
                    observations[falsePositive[i]] = true;
                    activateAncestors(falsePositive[i], observations); //they got turned on
                    //by the node beiong positive
                }
            } else {
                /* False positive */
                for (i = 0; i < numFalsePositive; i++) {
                    observations[falsePositive[i]] = true;
                }

                /* fix for the true path rule (reverse case) */
                for (i = 0; i < observations.length; i++) {
                    if (!observations[i]) { //i am off, so my descendants should be off???
                        //yes. the reason is that I need all my parents to be on in order for me to be on
                        //hence it really is the reverse case (q_parents = 0 me = 1) => 0 chance
                        //alternatively stated, if i am off,

                        //however, while they are all quite similar cases, we must check that
                        // the total meaning is correct: P=>Q vs Q=>P
                        //here, this seems to be saying that if I am off, then my descendants must be off
                        //which is NOT the same as saying if my descendants are off then I am off
                        //however, it might be bi-implicative since it specified a conditional dist, with
                        //no directionality
                        //also, we always require the concept of a "node" we are checking
                        //if i am off then for my descendants, Q_pa-and(i) will be 0, hence there is NO
                        //situation where they can be ON
                        deactivateDecendants(i, observations);
                    }
                }
            }
            //actually maxterms is -1 here (i.e. there is limit on how many terms an item can be annotated
            //to). FALSE actually what it is saying is how many terms can comprise a query at once.
            if (this.maxTerms != -1) {
                IntArray sparse = new IntArray(observations);
                int[] mostSpecific = mostSpecificTerms(sparse.get()); //get an array of only
                //terms that were True (A[i] contains the index into the original array)
                if (mostSpecific.length > this.maxTerms) {
                    int[] newTerms = new int[this.maxTerms];

                    /* Now randomly choose maxTerms and place them in new Terms */
                    for (int j = 0; j < this.maxTerms; j++) {
                        int r = rnd.nextInt(mostSpecific.length - j);
                        newTerms[j] = mostSpecific[r];
                        mostSpecific[r] = mostSpecific[mostSpecific.length - j - 1]; /*
                                                                                      * Move last selectable term into
                                                                                      * the place of the chosen one
                                                                                      */
                    }
                    for (int j = 0; j < observations.length; j++) {
                        observations[j] = false;
                    }
                    for (int t : newTerms) {
                        observations[t] = true;
                        activateAncestors(t, observations);
                    }
                }
            }

            for (i = 0; i < hidden.length; i++) {
                if (hidden[i]) {
                    numHidden++; //i.e. whether or not hidden at i is true (by anno prop rule etc.)
                }
            }

            if (logger.isDebugEnabled()) {
                logger.debug("Number of terms that were missed in hidden: " + numMissedInHidden);
                logger.debug("Number of hidden positives:" + numPositive);
                logger.debug("Number of hidden negatives: " + numHidden);
            }

            numPositive = 0;
            numFalseNegative = 0;
            numFalsePositive = 0;
            for (i = 0; i < observations.length; i++) {
                if (observations[i]) {
                    if (!hidden[i]) {
                        numFalsePositive++;
                    }
                    numPositive++;
                } else {
                    if (hidden[i]) {
                        numFalseNegative++;
                    }
                }
            }

            logger.debug("Number of observed positives: {}", numPositive);
            logger.debug("Raw number of false positives: {}", numFalsePositive);
            logger.debug("Raw number of false negatives {}", numFalseNegative);

            if (numPositive == 0 && !this.ALLOW_EMPTY_OBSERVATIONS) {
                /* Queries with no query make no sense */
                retry++;
                continue;
            }

            Configuration stats = new Configuration();
            determineCases(observations, hidden, stats);

            if (logger.isDebugEnabled()) {
                logger.debug("Number of modelled false postives {} (alpha={}%)",
                    stats.getCases(Configuration.NodeCase.FALSE_POSITIVE), stats.falsePositiveRate());
                logger.debug("Number of modelled false negatives {}  (beta={}%)",
                    stats.getCases(Configuration.NodeCase.FALSE_NEGATIVE), stats.falseNegativeRate());
            }

            o = new Observations();
            o.item = item;
            o.observations = observations;
            o.observationStats = stats; //this doesnt seem to be completely necessary, since
            //we overwite o anyways..
            //however, one puyrpose is to link it back and then be able to return it
        } while (!this.ALLOW_EMPTY_OBSERVATIONS && retry++ < 50);
        return o;
    }

    /**
     * Deactivate the ancestors of the given node.
     *
     * @param i
     * @param observations
     */
    public void deactivateAncestors(int i, boolean[] observations)
    {
        for (int j = 0; j < this.term2Ancestors[i].length; j++) {
            observations[this.term2Ancestors[i][j]] = false;
        }
    }

    /**
     * Activates the ancestors of the given node.
     *
     * @param i
     * @param observations
     */
    public void activateAncestors(int i, boolean[] observations)
    {
        for (int j = 0; j < this.term2Ancestors[i].length; j++) {
            observations[this.term2Ancestors[i][j]] = true;
        }
    }

    /**
     * J: Deactivates the descendants of the given node.
     *
     * @param i
     * @param observations
     */
    private void deactivateDecendants(int i, boolean[] observations)
    {
        for (int j = 0; j < this.term2Descendants[i].length; j++) {
            observations[this.term2Descendants[i][j]] = false;
        }
    }

    /**
     * Calculates a "fingerprint" for the current data. Note that the fingerprint is not necessary unqiue but it should
     * be sufficient for the purpose.
     *
     * @return
     */
    private int fingerprint()
    {
        int fp = 0x3333;
        for (int i = 0; i < this.allItemList.size(); i++) {
            fp += this.allItemList.get(i).hashCode();
        }
        for (int i = 0; i < this.slimGraph.getNumberOfVertices(); i++) {
            fp += this.slimGraph.getVertex(i).getID().id;
            fp += this.slimGraph.getVertex(i).getName().hashCode();
        }
        fp += new Random(this.SIZE_OF_SCORE_DISTRIBUTION).nextInt();
        fp += new Random(this.MAX_QUERY_SIZE_FOR_CACHED_DISTRIBUTION).nextInt();
        return fp;
    }

    /**
     * Setups the BOQA for the given ontology and associations.
     *
     * @param ontology
     * @param associations
     */
    public void setup(Ontology ontology, AssociationContainer associations)
    {
        this.assoc = associations; //this is the main work
        this.graph = ontology;

        // graph.findRedundantISARelations();

        if (this.micaMatrix != null) {
            logger.error("setup() called a 2nd time.");
            this.micaMatrix = null;
        }

        //which diseases to consider
        HashSet<ByteString> itemsToBeConsidered = new HashSet<ByteString>(associations.getAllAnnotatedGenes());
        provideGlobals(itemsToBeConsidered);

        /*
         * If we want to consider items with frequencies only, we like to shrink the item list to contain only the
         * relevant items.
         */
        if (this.CONSIDER_FREQUENCIES_ONLY) {
            int oldSize = this.allItemList.size();

            itemsToBeConsidered = new HashSet<ByteString>();
            for (int i = 0; i < this.allItemList.size(); i++) {
                if (this.itemHasFrequencies[i]) {
                    itemsToBeConsidered.add(this.allItemList.get(i));
                }
            }
            if (itemsToBeConsidered.size() == 0) {
                throw new RuntimeException("No items left after frequency filtering");
            }
            provideGlobals(itemsToBeConsidered);

            System.out.println("There were " + oldSize + " items but we consider only " + this.allItemList.size()
                + " of them with frequencies.");
            System.out.println("Considering " + this.slimGraph.getNumberOfVertices() + " terms");
        }

        /** Here we precalculate the jaccard similiartiy of two given terms in a dense matrix */
        if (this.PRECALCULATE_JACCARD) {
            logger.info("Calculating Jaccard");
            double[][] newJaccardMatrix = new double[this.slimGraph.getNumberOfVertices()][];
            for (int i = 0; i < this.slimGraph.getNumberOfVertices(); i++) {
                newJaccardMatrix[i] = new double[this.slimGraph.getNumberOfVertices() - i - 1];
                for (int j = i + 1; j < this.slimGraph.getNumberOfVertices(); j++) {
                    newJaccardMatrix[i][j - i - 1] = jaccard(i, j);
                }
            }
            this.jaccardMatrix = newJaccardMatrix;
            logger.info("Calculated Jaccard");
        }

        /** Here we precalculate the maxICs of two given terms in a dense matrix */
        if (this.PRECALCULATE_MAXICS) {
            logger.info("Calculating max ICs");
            int[][] newMaxICMatrix = new int[this.slimGraph.getNumberOfVertices()][];
            for (int i = 0; i < this.slimGraph.getNumberOfVertices(); i++) {
                newMaxICMatrix[i] = new int[this.slimGraph.getNumberOfVertices() - i - 1];
                for (int j = i + 1; j < this.slimGraph.getNumberOfVertices(); j++) {
                    newMaxICMatrix[i][j - i - 1] = commonAncestorWithMaxIC(i, j);
                }
            }
            this.micaMatrix = newMaxICMatrix;

            logger.info("Calculated max ICs");
        }

        /**
         * Here we precalculate for each item the term which contributes as maximum ic term to the resnick calculation
         */
        if (this.PRECALCULATE_ITEM_MAXS) {
            logger.info("Calculating item maxs");
            this.resnikTermSim.maxScoreForItem =
                new double[this.allItemList.size()][this.slimGraph.getNumberOfVertices()];
            this.linTermSim.maxScoreForItem = new double[this.allItemList.size()][this.slimGraph.getNumberOfVertices()];
            this.jcTermSim.maxScoreForItem = new double[this.allItemList.size()][this.slimGraph.getNumberOfVertices()];

            for (int item = 0; item < this.allItemList.size(); item++) {
                /* The fixed set */
                int[] t2 = this.items2DirectTerms[item];

                /* The set representing a single query term */
                int[] t1 = new int[1];

                for (int to = 0; to < this.slimGraph.getNumberOfVertices(); to++) {
                    t1[0] = to;
                    this.resnikTermSim.maxScoreForItem[item][to] = scoreMaxAvg(t1, t2, this.resnikTermSim);
                    this.linTermSim.maxScoreForItem[item][to] = scoreMaxAvg(t1, t2, this.linTermSim);
                    this.jcTermSim.maxScoreForItem[item][to] = scoreMaxAvg(t1, t2, this.jcTermSim);
                }
            }

            logger.info("Calculated item maxs");
        }

        this.resnikTermSim.setupDistribution();
        this.linTermSim.setupDistribution();
        this.jcTermSim.setupDistribution();

        /* Choose appropriate values */
        double numOfTerms = getSlimGraph().getNumberOfVertices();
        //overwrite the other values!
        this.ALPHA_GRID =
            new double[] { 1e-10, 1 / numOfTerms, 2 / numOfTerms, 3 / numOfTerms, 4 / numOfTerms, 5 / numOfTerms,
            6 / numOfTerms };
        this.BETA_GRID = new double[] { 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95 };
    }

    /**
     * Writes a dot suitable for tikz.
     *
     * @param out
     * @param hpoTerms
     */
    public void writeDOTExample(File out, HashSet<TermID> hpoTerms)
    {
        /*
         * Basically, this defines a new command \maxbox whose text width as given by the second argument is not wider
         * than the first argument. The text which is then displayed in the box is used from the third argument.
         */
        String preamble = "d2tfigpreamble=\"\\ifthenelse{\\isundefined{\\myboxlen}}{\\newlength{\\myboxlen}}{}" +
            "\\newcommand*{\\maxbox}[3]{\\settowidth{\\myboxlen}{#2}" +
            "\\ifdim#1<\\myboxlen" +
            "\\parbox{#1}{\\centering#3}" +
            "\\else" +
            "\\parbox{\\myboxlen}{\\centering#3}" +
            "\\fi}\"";

        try {
            GODOTWriter.writeDOT(this.graph.getInducedGraph(this.termEnumerator.getAllAnnotatedTermsAsList()), out,
                null,
                hpoTerms, new AbstractDotAttributesProvider()
                {
                    @Override
                    public String getDotNodeAttributes(TermID id)
                    {
                        String termName;
                        Term term = BOQA.this.graph.getTerm(id);
                        if (BOQA.this.graph.isRootTerm(id)) {
                            termName = "Human Phenotype";
                        } else {
                            termName = term.getName();
                        }
                        String name = "\\emph{" + termName + "}";
                        int termIdx = BOQA.this.slimGraph.getVertexIndex(BOQA.this.graph.getTerm(id));
                        int numberOfItems = BOQA.this.termEnumerator.getAnnotatedGenes(id).totalAnnotatedCount();

                        String label =
                            "\\small" + name + "\\\\\\ " + numberOfItems + " \\\\\\ IC="
                                + String.format("%.4f", BOQA.this.terms2IC[termIdx]);
                        return "margin=\"0\" shape=\"box\"" + " label=\"\\maxbox{4.5cm}{" + name + "}{" + label
                            + "}\" " +
                            "style=\"rounded corners,top color=white,bottom color=black!10,draw=black!50,very thick\"";
                    }
                }, "nodesep=0.2; ranksep=0.1;" + preamble, false, false, null);
        } catch (IllegalArgumentException ex) {
            logger.error("Failed to write graphics due to: {}", ex.getLocalizedMessage());
        }
    }

    /**
     * Write score distribution.
     *
     * @throws IOException
     */
    public void writeScoreDistribution(File f, int item) throws IOException
    {
        int[] shuffledTerms = new int[this.slimGraph.getNumberOfVertices()];

        /* Initialize shuffling */
        for (item = 0; item < shuffledTerms.length; item++) {
            shuffledTerms[item] = item;
        }

        try (FileWriter out = new FileWriter(f)) {

            Random rnd = new Random();

            for (int j = 0; j < this.SIZE_OF_SCORE_DISTRIBUTION; j++) {
                int q = 10;
                int[] randomizedTerms = new int[q];

                chooseTerms(rnd, q, randomizedTerms, shuffledTerms);
                double randomScore = resScoreVsItem(randomizedTerms, item);
                out.write(randomScore + " \n");
            }
            out.flush();
        }

        logger.info("Score distribution for item {} with {} annotations written",
            this.allItemList.get(item), this.items2DirectTerms[item].length);
    }

    /**
     * @return
     */
    public static int getNumProcessors()
    {
        int numProcessors = Runtime.getRuntime().availableProcessors();
        return numProcessors;
    }

    /**
     * Calculates the set difference of a minus b.
     *
     * @param a
     * @param b
     * @return
     */
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

    /**
     * Helper function to create sub array of length elements from the given array.
     * this copies over the firrst n elements (i.e. subarray from 0 to length)
     * @param array
     * @param length
     * @return
     * TODO move this to a new class etc.
     */
    private static int[] subArray(int[] array, int length)
    {
        int[] a = new int[length];
        for (int i = 0; i < length; i++) {
            a[i] = array[i];
        }
        return a;
    }

    /**
     * A simple class maintaining ints.
     *
     * @author Sebastian Bauer
     */
    public static class IntArray
    {
        private int[] array;

        private int length;

        public IntArray(int maxLength)
        {
            this.array = new int[maxLength];
        }

        public IntArray(boolean[] dense)
        {
            int c = 0;
            for (boolean element : dense) {
                if (element) {
                    c++;
                }
            }

            this.array = new int[c]; //on the first pass, get the number of elements on in the dense array
            //on the next pass, actually put in all them selements
            c = 0; //my claim is that c WILL be the size of the array at the end, and hence we cannot just add
            //to it
            //only put the
            for (int i = 0; i < dense.length; i++) {
                if (dense[i]) {
                    this.array[c++] = i; //essentially this forms a coalescence operation
                    //(Brings all the true ones to a new array; we are guaranteed to just have an array of
                    //trues at the end)
                }
            }
            this.length = c;
        }

        //@TODO This is a deprecated method and will be removed
        public void add(int e)
        {
            this.array[this.length++] = e;
        } //should this not be larger?
        //we might overfill our array! -- i.e. arrayout of bounds exception?!
        //not necessarily, (consider just calling it using the length constructor, then it just fills
        //the array as regular


        public int[] get()//returns the first n elements (i.e. the entire array)
        {
            return subArray(this.array, this.length);
        }
    }

    /**
     * Provides some global variables, given the global graph, the global associations and the items.
     *
     * @param allItemsToBeConsidered
     */
    @SuppressWarnings("unused")
    private void provideGlobals(Set<ByteString> allItemsToBeConsidered)
    {
        int i;

        /* list all evidence codes */
        HashMap<ByteString, Integer> evidences = new HashMap<ByteString, Integer>();
        for (Gene2Associations g2a : this.assoc) {
            for (Association a : g2a) {
                //for each association in g2a
                if (a.getEvidence() != null) { //this is always null
                    /* Worst implementation ever! */
                    Integer evidence = evidences.get(a.getEvidence());
                    if (evidence == null) {
                        evidence = 1;
                    } else {
                        evidence++;
                    }

                    evidences.put(a.getEvidence(), evidence);
                }
            }
        }

        if (logger.isInfoEnabled()) {
            logger.info(allItemsToBeConsidered.size() + " items shall be considered");
            StringBuilder builder = new StringBuilder("Available evidences: ");
            for (Entry<ByteString, Integer> ev : evidences.entrySet()) {
                builder.append(ev.getKey().toString() + "->" + ev.getValue() + ",");
            }
            logger.info(builder.toString());
        }

        if (this.evidenceCodes != null) {
            evidences.clear();
            for (String ev : this.evidenceCodes) {
                evidences.put(new ByteString(ev), 1);
            }

            if (logger.isInfoEnabled()) {
                StringBuilder builder = new StringBuilder("Requested evidences: ");
                for (ByteString ev : evidences.keySet()) {
                    builder.append(ev.toString());
                }
                logger.info(builder.toString());
            }
        } else {
            /* Means take everything */
            evidences = null;
        }

        PopulationSet allItems = new PopulationSet("all");
        allItems.addGenes(allItemsToBeConsidered);
        this.termEnumerator =  //since a.getEvidence is always false
            allItems.enumerateGOTerms(this.graph, this.assoc, evidences != null ? evidences.keySet() : null);
        //a
        ItemEnumerator itemEnumerator = ItemEnumerator.createFromTermEnumerator(this.termEnumerator);

        /* Term stuff */
        Ontology inducedGraph = this.graph.getInducedGraph(this.termEnumerator.getAllAnnotatedTermsAsList());
        this.slimGraph = inducedGraph.getSlimGraphView();

        this.term2Parents = this.slimGraph.vertexParents; //recall this was an array of arrays
        this.term2Children = this.slimGraph.vertexChildren;
        this.term2Ancestors = this.slimGraph.vertexAncestors;
        this.term2Descendants = this.slimGraph.vertexDescendants; //presumably, this is a vertex induced subgraph
        this.termsInTopologicalOrder = this.slimGraph.getVertexIndices(inducedGraph.getTermsInTopologicalOrder());

        if (this.termsInTopologicalOrder.length != this.slimGraph.getNumberOfVertices()) {
            throw new RuntimeException("The ontology graph contains cycles.");
        }
        this.termsToplogicalRank = new int[this.termsInTopologicalOrder.length];
        for (i = 0; i < this.termsInTopologicalOrder.length; i++) {
            this.termsToplogicalRank[this.termsInTopologicalOrder[i]] = i; //this is just copying them over
        }

        /* Item stuff */
        this.allItemList = new ArrayList<ByteString>();
        this.item2Index = new HashMap<ByteString, Integer>();
        i = 0;
        for (ByteString item : itemEnumerator) {
            this.allItemList.add(item);
            this.item2Index.put(item, i);
            i++;
        }

        logger.info(i + " items passed criterias (supplied evidence codes)");

        /* Fill item matrix */
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

        /* Fill direct item matrix */
        this.items2DirectTerms = new int[this.allItemList.size()][];
        i = 0;
        for (ByteString item : itemEnumerator) {
            int j = 0;

            ArrayList<TermID> tids = itemEnumerator.getTermsDirectlyAnnotatedToTheItem(item);

            if (false) {
                /* Perform sanity check */
                for (TermID s : tids) {
                    for (TermID d : tids) {
                        if (this.graph.existsPath(s, d) || this.graph.existsPath(d, s)) {
                            System.out.println("Item \"" + item + "\" is annotated to " + s.toString() + " and "
                                + d.toString());
                        }
                    }
                }
            }

            // S: System.out.println(item.toString());
            this.items2DirectTerms[i] = new int[tids.size()];

            for (TermID tid : tids) {
                this.items2DirectTerms[i][j++] = this.slimGraph.getVertexIndex(this.graph.getTerm(tid));
            }
            i++;
        }

        /* Fill in frequencies for directly annotated terms. Also sort them */
        this.items2TermFrequencies = new double[this.allItemList.size()][];
        this.itemHasFrequencies = new boolean[this.allItemList.size()];
        this.item2TermFrequenciesOrder = new int[this.allItemList.size()][];
        for (i = 0; i < this.items2DirectTerms.length; i++) {
            /**
             * A term and the corresponding frequency. We use this for sorting.
             *
             * @author Sebastian Bauer
             */
            class Freq implements Comparable<Freq>
            {
                public int termIdx;

                public double freq;

                @Override
                public int compareTo(Freq o)
                {
                    if (this.freq > o.freq) {
                        return 1;
                    }
                    if (this.freq < o.freq) {
                        return -1;
                    }
                    return 0;
                }
            }

            this.items2TermFrequencies[i] = new double[this.items2DirectTerms[i].length];
            this.item2TermFrequenciesOrder[i] = new int[this.items2DirectTerms[i].length];
            Freq[] freqs = new Freq[this.items2DirectTerms[i].length];

            ByteString item = this.allItemList.get(i);
            Gene2Associations as = this.assoc.get(item);

            // Disabled
            // if (as.getAssociations().size() != items2DirectTerms[i].length)
            // throw new IllegalArgumentException("Number of associations differs (" + as.getAssociations().size() +
            // ") from the number of directly annotated terms (" + items2DirectTerms[i].length + ").");

            for (int j = 0; j < this.items2DirectTerms[i].length; j++) {
                boolean hasExlipictFrequency = false;

                /* Default frequency */
                double f = 1.0;

                TermID tid = this.slimGraph.getVertex(this.items2DirectTerms[i][j]).getID();

                /* Find frequency. We now have a O(n^3) algo. Will be optimized later */
                for (Association a : as) {
                    if (a.getTermID().equals(tid) && a.getAspect() != null) {
                        f = getFrequencyFromString(a.getAspect().toString());
                        if (f < 1.0) {
                            hasExlipictFrequency = true;
                        }
                        /* We assume that the term appears only once */
                        break;
                    }
                }

                this.items2TermFrequencies[i][j] = f;
                freqs[j] = new Freq();
                freqs[j].termIdx = j;// items2DirectTerms[i][j];
                freqs[j].freq = f;

                if (hasExlipictFrequency) {
                    this.itemHasFrequencies[i] = true;
                }
            }

            /* Now sort and remember the indices */
            Arrays.sort(freqs);
            for (int j = 0; j < this.items2DirectTerms[i].length; j++) {
                this.item2TermFrequenciesOrder[i][j] = freqs[j].termIdx;
            }
        }

        createDiffVectors();

        /* Calculate IC */
        this.terms2IC = new double[this.slimGraph.getNumberOfVertices()];
        for (i = 0; i < this.slimGraph.getNumberOfVertices(); i++) {
            Term t = this.slimGraph.getVertex(i);
            this.terms2IC[i] =
                -Math
                    .log(((double) this.termEnumerator.getAnnotatedGenes(t.getID()).totalAnnotatedCount()
                        / this.allItemList
                            .size()));
        }

        ArrayList<Integer> itemIndices = new ArrayList<Integer>();
        for (int o = 0; o < this.allItemList.size(); o++) {
            itemIndices.add(o);
        }

        if (false) {
            System.out.println("Start TSP");
            long start = System.nanoTime();
            Algorithms.approximatedTSP(itemIndices, itemIndices.get(0),
                new Algorithms.IVertexDistance<Integer>()
                {
                    @Override
                    public double distance(Integer ai, Integer bi)
                    {
                        int[] at = BOQA.this.items2Terms[ai.intValue()];
                        int[] bt = BOQA.this.items2Terms[bi.intValue()];
                        return Algorithms.hammingDistanceSparse(at, bt);
                    }
                });
            System.out.println("End (" + ((System.nanoTime() - start) / 1000 / 1000) + "ms)");
        }
    }

    /**
     * Create the diff annotation vectors.
     */
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
        logger.info(sum + " differences detected (" + (double) sum / this.allItemList.size() + " per item)");

        this.diffOnTermsFreqs = new int[this.allItemList.size()][][];
        this.diffOffTermsFreqs = new int[this.allItemList.size()][][];
        this.factors = new double[this.allItemList.size()][];
        for (int item = 0; item < this.allItemList.size(); item++) {
            int numTerms = this.items2TermFrequencies[item].length;
            int numTermsWithExplicitFrequencies = 0;
            int numConfigs = 0;

            /*
             * Determine the number of terms that have non-1.0 frequency. We restrict them to the top 6 (the less
             * probable) due to complexity issues and hope that this a good enough approximation.
             */
            for (i = 0; i < numTerms && i < this.maxFrequencyTerms; i++) {
                if (this.items2TermFrequencies[item][this.item2TermFrequenciesOrder[item][i]] >= 1.0) {
                    break;
                }
                numTermsWithExplicitFrequencies++;
            }

            /* We try each possible activity/inactivity combination of terms with explicit frequencies */
            SubsetGenerator sg = new SubsetGenerator(numTermsWithExplicitFrequencies, numTermsWithExplicitFrequencies);
            SubsetGenerator.Subset s;

            /* First, determine the number of configs (could calculate binomial coefficient of course) */
            while ((s = sg.next()) != null) {
                numConfigs++;                 //since this is an iterator, this is the only way to get size
            }

            this.diffOnTermsFreqs[item] = new int[numConfigs][];
            this.diffOffTermsFreqs[item] = new int[numConfigs][];
            this.factors[item] = new double[numConfigs];

            /* Contains the settings of the previous run */
            IntArray prevArray = new IntArray(this.slimGraph.getNumberOfVertices());

            int config = 0;

            while ((s = sg.next()) != null) {
                boolean[] hidden = new boolean[this.slimGraph.getNumberOfVertices()];
                boolean[] taken = new boolean[numTermsWithExplicitFrequencies];

                double factor = 0.0;

                /* First, activate variable terms according to the current selection */
                for (i = 0; i < s.r; i++) {
                    int ti = this.item2TermFrequenciesOrder[item][s.j[i]]; /*
                                                                            * index of term within the all directly
                                                                            * associated indices
                                                                            */
                    int h = this.items2DirectTerms[item][ti]; /* global index of term */
                    hidden[h] = true;
                    activateAncestors(h, hidden);
                    factor += Math.log(this.items2TermFrequencies[item][ti]);
                    taken[s.j[i]] = true;
                }

                /* Needs also respect the inactive terms in the factor */
                for (i = 0; i < numTermsWithExplicitFrequencies; i++) {
                    if (!taken[i]) {
                        factor +=
                            Math.log(1 - this.items2TermFrequencies[item][this.item2TermFrequenciesOrder[item][i]]);
                    }
                }

                /* Second, activate mandatory terms */
                for (i = numTermsWithExplicitFrequencies; i < numTerms; i++) {
                    int ti = this.item2TermFrequenciesOrder[item][i];
                    int h = this.items2DirectTerms[item][ti]; /* global index of term */
                    hidden[h] = true;
                    activateAncestors(h, hidden);
                    /* Factor is always 0 */
                }

                /* Now make a sparse representation */
                IntArray newArray = new IntArray(hidden);

                /* And record the difference */
                this.diffOnTermsFreqs[item][config] = setDiff(newArray.get(), prevArray.get());
                this.diffOffTermsFreqs[item][config] = setDiff(prevArray.get(), newArray.get());
                this.factors[item][config] = factor;

                prevArray = newArray;
                config++;
            }
        }
    }

    /**
     * Returns the number of items that are annotated to term i.
     *
     * @param i
     * @return
     */
    public int getNumberOfItemsAnnotatedToTerm(int i)
    {
        Term t = this.slimGraph.getVertex(i);
        return this.termEnumerator.getAnnotatedGenes(t.getID()).totalAnnotatedCount();
    }

    /**
     * Return the IC of term i
     *
     * @param i
     * @return
     */
    public double getIC(int i)
    {
        return this.terms2IC[i];
    }

    /**
     * Converts the frequency string to a double value.
     *
     * @param freq
     * @return
     */
    private double getFrequencyFromString(String freq)
    {
        double f = 1.0;

        if (freq == null || freq.length() == 0) {
            return 1.0;
        }

        Matcher matcher = frequencyPattern.matcher(freq);
        if (matcher.matches()) {
            String fractionalPart = matcher.group(2);
            if (fractionalPart == null || fractionalPart.length() == 0) {
                fractionalPart = "0";
            }

            f =
                Double.parseDouble(matcher.group(1)) + Double.parseDouble(fractionalPart)
                    / Math.pow(10, fractionalPart.length());
            f /= 100.0;
        } else {
            matcher = frequencyFractionPattern.matcher(freq);
            if (matcher.matches()) {
                f = Double.parseDouble(matcher.group(1)) / Double.parseDouble(matcher.group(2));
            } else {
                if (freq.equalsIgnoreCase("very rare")) {
                    f = 0.01;
                } else if (freq.equalsIgnoreCase("rare")) {
                    f = 0.05;
                } else if (freq.equalsIgnoreCase("occasional")) {
                    f = 0.075;
                } else if (freq.equalsIgnoreCase("frequent")) {
                    f = 0.33;
                } else if (freq.equalsIgnoreCase("typical")) {
                    f = 0.50;
                } else if (freq.equalsIgnoreCase("common")) {
                    f = 0.75;
                } else if (freq.equalsIgnoreCase("hallmark")) {
                    f = 0.90;
                } else if (freq.equalsIgnoreCase("obligate")) {
                    f = 1;
                } else {
                    logger.info("Unknown frequency identifier: {}", freq);
                }
            }
        }
        return f;
    }

    /**
     * This is a container for the results of the class.
     *
     * @author Sebastian Bauer
     */
    static public class Result
    {
        /** Contains the marginal probability for each item */
        private double[] marginals;

        /** Contains the marginal probability for each item */
        private double[] marginalsIdeal;

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

    /**
     * Provides the marginals for the observations. (non threaded--builder paradigm?)
     *
     * @param observations
     * @param takeFrequenciesIntoAccount
     * @return
     */
    public Result assignMarginals(Observations observations, boolean takeFrequenciesIntoAccount)
    {
        return assignMarginals(observations, takeFrequenciesIntoAccount, 1);
    }

    /**
     * Provides the marginals for the observations.
     *
     * @param observations
     * @param takeFrequenciesIntoAccount
     * @param numThreads defines the number of threads to be used for the calculation.
     * @return
     */
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
            res.stats[i] = new Configuration(); //each node gets its own cofnig?
            //In fact each NODE in the BN gets its own config? (it's confusing since stats
            //is also a variable inside configuration...)
            //however the stats inside config is more "true" to what it means (the 4 buckets)
        }
        for (i = 0; i < res.scores.length; i++) {
            res.scores[i] = Double.NEGATIVE_INFINITY; //-Math.
            //Math.log(0); //this is etting it to negative infinity??
            //for example, one speedup is to not compute log(0) each time
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
        final Configuration previousStat = new Configuration();

        //Configuration is not null, but: it also was not initalized to anything
        determineCases(observations.observations, previousHidden, previousStat);
        //this pops back with all the cases; and having incremented the particular stats
        //(i.e previousStat)

        ArrayList<Future<?>> futureList = new ArrayList<Future<?>>();

        for (i = 0; i < this.allItemList.size(); i++) {
            final int item = i;

            /* Construct the runnable suitable for the calculation for a single item */
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

                    for (WeightedConfiguration wc : stats) //this is a list of weighted configs,
                        // which itself is simply a tuple of factor and config
                    {
                        System.out.println("Itme " + item);
                        System.out.println(wc.stat);

                    }
                    //J: the scoring function is critical; also this seems to be finding the best
                    //rate for alpha and beta.at a particular alpha, beta value
                    //this puts the score in for an item;
                    for (int a = 0; a < BOQA.this.ALPHA_GRID.length; a++) {
                        for (int b = 0; b < BOQA.this.BETA_GRID.length; b++) {
                            //This scor is a weighted sum!
                            scores[item][a][b] = stats.score(BOQA.this.ALPHA_GRID[a],
                                    BOQA.this.BETA_GRID[b]);
                            res.scores[item] = Util.logAdd(res.scores[item], scores[item][a][b]);
                        }
                    }

                    /* This is used only for benchmarks, where we know the true configuration */
                    if (observations.observationStats != null) {
                        /* Calculate ideal scores */
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

        /*
         * There is a possibility that ideal marginal is not as good as the marginal for the unknown parameter
         * situation, i.e., if the initial signal got such disrupted that another item is more likely. This may produce
         * strange plots. Therefore, we take the parameter estimated marginals as the ideal one if they match the
         * reality better.
         */
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

    static long time;

    static long lastTime;

    /**
     * Return a common ancestor of t1 and t2 that have max ic.
     *
     * @param t1
     * @param t2
     * @return
     */
    private int commonAncestorWithMaxIC(int t1, int t2)
    {
        if (this.micaMatrix != null) {
            if (t1 < t2) {
                return this.micaMatrix[t1][t2 - t1 - 1];
            } else if (t2 < t1) {
                return this.micaMatrix[t2][t1 - t2 - 1];
            } else {
                return t1;
            }
        }

        /* A rather slow implementation */ // TODO: 16/05/17
        int[] ancestorsA;
        int[] ancestorsB;

        if (this.term2Ancestors[t1].length > this.term2Ancestors[t2].length) {
            ancestorsA = this.term2Ancestors[t1];
            ancestorsB = this.term2Ancestors[t2];
        } else {
            ancestorsA = this.term2Ancestors[t1];
            ancestorsB = this.term2Ancestors[t2];
        }

        int bestTerm = -1;
        double bestIC = Double.NEGATIVE_INFINITY;

        for (int term : ancestorsA) {
            for (int element : ancestorsB) {
                if (term == element) {
                    double ic = this.terms2IC[term];

                    if (ic > bestIC) {
                        bestIC = ic;
                        bestTerm = term;
                    }
                    break;
                }
            }
        }

        if (bestTerm == -1) {
            throw new RuntimeException("No best term found, which is strange.");
        }

        return bestTerm;
    }

    /**
     * Returns the jaccard index of the given two terms.
     *
     * @param t1
     * @param t2
     * @return
     */
    public double jaccard(int t1, int t2)
    {
        if (t1 == t2) {
            return 1;
        }

        if (this.jaccardMatrix != null) {
            if (t1 < t2) {
                return this.jaccardMatrix[t1][t2 - t1 - 1];
            } else {
                return this.jaccardMatrix[t2][t1 - t2 - 1];
            }
        }

        Term tt1 = this.slimGraph.getVertex(t1);
        Term tt2 = this.slimGraph.getVertex(t2);
        HashSet<ByteString> tt1a =
            new HashSet<ByteString>(this.termEnumerator.getAnnotatedGenes(tt1.getID()).totalAnnotated);
        HashSet<ByteString> tt2a =
            new HashSet<ByteString>(this.termEnumerator.getAnnotatedGenes(tt2.getID()).totalAnnotated);
        HashSet<ByteString> union = new HashSet<ByteString>(tt1a);
        union.addAll(tt2a);

        tt1a.retainAll(tt2a);

        return (double) tt1a.size() / union.size();
    }

    /**
     * Returns a minimal length array of terms of which the induced graph is the same as of the given terms. These are
     * the leaf terms.
     *
     * @param terms
     * @return
     */
    public int[] mostSpecificTerms(int[] terms)
    {
        ArrayList<TermID> termList = new ArrayList<TermID>(terms.length);
        for (int term : terms) {
            termList.add(this.slimGraph.getVertex(term).getID());
        }

        Ontology termGraph = this.graph.getInducedGraph(termList);

        ArrayList<Term> leafTermList = termGraph.getLeafTerms();

        int[] specifcTerms = new int[leafTermList.size()];
        int i = 0;

        for (Term t : termGraph.getLeafTerms()) {
            specifcTerms[i++] = this.slimGraph.getVertexIndex(t);
        }

        return specifcTerms;
    }

    /**
     * Gets a sparse representation of the most specific terms in the observation map.
     *
     * @param observations
     * @return
     */
    private int[] getMostSpecificTermsSparse(boolean[] observations)
    {
        int numObservedTerms = 0;
        for (boolean observation : observations) {
            if (observation) {
                numObservedTerms++;
            }
        }

        int[] observedTerms = new int[numObservedTerms];
        for (int i = 0, j = 0; i < observations.length; i++) {
            if (observations[i]) {
                observedTerms[j++] = i;
            }
        }

        return mostSpecificTerms(observedTerms);
    }

    /**
     * Defines a function to determine the term similarity.
     *
     * @author Sebastian Bauer
     */
    public static interface ITermSim
    {
        public double termSim(int t1, int t2);

        /**
         * @return the name of the method.
         */
        public String name();
    }

    /**
     * Class implementing a term similarity measure.
     *
     * @author Sebastian Bauer
     */
    public abstract class AbstractTermSim implements ITermSim
    {
        /** Contains for each item the maximal score for the given term */
        public double[][] maxScoreForItem;

        /** Stores the score distribution */
        private ApproximatedEmpiricalDistributions scoreDistributions;

        /** Lock for the score distribution */
        private ReentrantReadWriteLock scoreDistributionLock = new ReentrantReadWriteLock();

        /**
         * Returns the score distribution for the given item for the given query size. If the score distribution has not
         * been created yet, create it using the supplied queries.
         *
         * @param querySize
         * @param item
         * @param queries
         * @return
         */
        private ApproximatedEmpiricalDistribution getScoreDistribution(int querySize, int item, int[][] queries)
        {
            this.scoreDistributionLock.readLock().lock();
            ApproximatedEmpiricalDistribution d =
                this.scoreDistributions.getDistribution(item * (BOQA.this.MAX_QUERY_SIZE_FOR_CACHED_DISTRIBUTION + 1)
                    + querySize);
            this.scoreDistributionLock.readLock().unlock();

            if (d == null) {
                /* Determine score distribution */
                double[] scores = new double[BOQA.this.SIZE_OF_SCORE_DISTRIBUTION];
                double maxScore = Double.NEGATIVE_INFINITY;

                for (int j = 0; j < BOQA.this.SIZE_OF_SCORE_DISTRIBUTION; j++) {
                    scores[j] = scoreVsItem(queries[j], item, this);
                    if (scores[j] > maxScore) {
                        maxScore = scores[j];
                    }
                }

                ApproximatedEmpiricalDistribution d2 =
                    new ApproximatedEmpiricalDistribution(scores,
                        BOQA.this.NUMBER_OF_BINS_IN_APPROXIMATED_SCORE_DISTRIBUTION);

                this.scoreDistributionLock.writeLock().lock();
                d =
                    this.scoreDistributions.getDistribution(item
                        * (BOQA.this.MAX_QUERY_SIZE_FOR_CACHED_DISTRIBUTION + 1) + querySize);
                if (d == null) {
                    this.scoreDistributions.setDistribution(item
                        * (BOQA.this.MAX_QUERY_SIZE_FOR_CACHED_DISTRIBUTION + 1) + querySize,
                        d2);
                }
                this.scoreDistributionLock.writeLock().unlock();
            }

            return d;
        }

        /**
         * Sets up the score distribution. At the moment, this must be called before maxScoreForItem is setup.
         */
        public void setupDistribution()
        {
            /** Instantiates the query cache */
            if (BOQA.this.CACHE_RANDOM_QUERIES) {
                boolean distributionLoaded = false;
                String scoreDistributionsName =
                    "scoreDistributions-" + name() + "-" + BOQA.this.allItemList.size() + "-"
                        + BOQA.this.CONSIDER_FREQUENCIES_ONLY + "-"
                        + BOQA.this.SIZE_OF_SCORE_DISTRIBUTION + ".gz";

                BOQA.this.queryCache = new QuerySets(BOQA.this.MAX_QUERY_SIZE_FOR_CACHED_DISTRIBUTION + 1);

                if ((BOQA.this.CACHE_SCORE_DISTRIBUTION || BOQA.this.PRECALCULATE_SCORE_DISTRIBUTION)
                    && BOQA.this.TRY_LOADING_SCORE_DISTRIBUTION) {
                    File inFile = new File(scoreDistributionsName);
                    try (InputStream underlyingStream = new GZIPInputStream(new FileInputStream(inFile));
                        ObjectInputStream ois = new ObjectInputStream(underlyingStream)) {

                        int fingerprint = ois.readInt();
                        if (fingerprint == fingerprint()) {
                            this.scoreDistributions = (ApproximatedEmpiricalDistributions) ois.readObject();
                            distributionLoaded = true;
                            logger.info("Score distribution loaded from \"{}\"", inFile.getAbsolutePath());
                        }
                    } catch (FileNotFoundException e) {
                    } catch (IOException e) {
                        e.printStackTrace();
                    } catch (ClassNotFoundException e) {
                        e.printStackTrace();
                    }
                }

                if (BOQA.this.PRECALCULATE_SCORE_DISTRIBUTION) {
                    if (!distributionLoaded) {
                        this.scoreDistributions =
                            new ApproximatedEmpiricalDistributions(BOQA.this.allItemList.size()
                                * (BOQA.this.MAX_QUERY_SIZE_FOR_CACHED_DISTRIBUTION + 1));
                    }

                    logger.info("Precaculating score distribution for {}", name());

                    Random rnd = new Random(9);
                    ExecutorService es = null;
                    ArrayList<Future<?>> futureList = new ArrayList<Future<?>>();

                    if (getNumProcessors() > 1) {
                        es = Executors.newFixedThreadPool(getNumProcessors());
                    } else {
                        es = null;
                    }

                    for (int i = 0; i < BOQA.this.allItemList.size(); i++) {
                        final long seed = rnd.nextLong();
                        final int item = i;

                        Runnable run = new Runnable()
                        {
                            @Override
                            public void run()
                            {
                                Random rnd = new Random(seed);

                                for (int qs = 1; qs <= BOQA.this.MAX_QUERY_SIZE_FOR_CACHED_DISTRIBUTION; qs++) {
                                    int[][] queries = getRandomizedQueries(rnd, qs);
                                    getScoreDistribution(qs, item, queries);
                                }
                            }
                        };

                        if (es != null) {
                            futureList.add(es.submit(run));
                        } else {
                            run.run();
                        }
                    }

                    /* Cleanup */
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
                            throw new RuntimeException(e);
                        }
                    }

                    logger.info("Score distribution has been precalculated");

                    if (BOQA.this.STORE_SCORE_DISTRIBUTION && !distributionLoaded) {
                        File outFile = new File(scoreDistributionsName);
                        try (OutputStream underlyingStream = new GZIPOutputStream(new FileOutputStream(outFile));
                            ObjectOutputStream oos = new ObjectOutputStream(underlyingStream)) {
                            /*
                             * The fingerprint shall ensure that the score distribution and ontology/associations are
                             * compatible
                             */
                            oos.writeInt(fingerprint());

                            /* Finally, Write store distribution */
                            oos.writeObject(this.scoreDistributions);

                            logger.info("Score distribution written to \"{}\"", outFile.getAbsolutePath());
                        } catch (IOException e) {
                            logger.warn("Failed to write score distribution: {}", e.getMessage(), e);
                        }
                    }
                }
            }

        }
    }

    /**
     * Term similarity measure according to Resnik
     */
    public final AbstractTermSim resnikTermSim = new AbstractTermSim()
    {
        @Override
        public double termSim(int t1, int t2)
        {
            return BOQA.this.terms2IC[commonAncestorWithMaxIC(t1, t2)];
        }

        @Override
        public String name()
        {
            return "resnik";
        }
    };

    /**
     * Term similarity measure according to Lin. Note that the similarity of terms with information content of 0 is
     * defined as 1 here.
     */
    private final AbstractTermSim linTermSim = new AbstractTermSim()
    {
        @Override
        public double termSim(int t1, int t2)
        {
            double nominator = 2 * BOQA.this.terms2IC[commonAncestorWithMaxIC(t1, t2)];
            double denominator = BOQA.this.terms2IC[t1] + BOQA.this.terms2IC[t2];
            if (nominator <= 0.0 && denominator <= 0.0) {
                return 1;
            }
            return nominator / denominator;
        }

        @Override
        public String name()
        {
            return "lin";
        };
    };

    /**
     * Term similarity measure according to Jiang and Conrath
     */
    private final AbstractTermSim jcTermSim = new AbstractTermSim()
    {
        @Override
        public double termSim(int t1, int t2)
        {
            return (1)
                / (1 + BOQA.this.terms2IC[t1] + BOQA.this.terms2IC[t2] - 2 * BOQA.this.terms2IC[commonAncestorWithMaxIC(
                    t1, t2)]);
        }

        @Override
        public String name()
        {
            return "jc";
        }
    };

    /**
     * Score two list of terms according to max-avg-of-best method using the given term similarity measure.
     *
     * @param tl1
     * @param tl2
     * @param termSim
     * @return
     */
    private double scoreMaxAvg(int[] tl1, int[] tl2, ITermSim termSim)
    {
        double totalScore = 0;
        for (int t1 : tl1) {
            double maxScore = Double.NEGATIVE_INFINITY;

            for (int t2 : tl2) {
                double score = termSim.termSim(t1, t2);
                if (score > maxScore) {
                    maxScore = score;
                }
            }

            totalScore += maxScore;
        }
        totalScore /= tl1.length;
        return totalScore;
    }

    /**
     * Score two list of terms according to resnik-max-avg-of-best method.
     *
     * @param tl1
     * @param tl2
     * @return
     */
    public double resScoreMaxAvg(int[] tl1, int[] tl2)
    {
        return scoreMaxAvg(tl1, tl2, this.resnikTermSim);
    }

    /**
     * Score two list of terms according to lin-max-avg-of-best method.
     *
     * @param tl1
     * @param tl2
     * @return
     */
    public double linScoreMaxAvg(int[] tl1, int[] tl2)
    {
        return scoreMaxAvg(tl1, tl2, this.linTermSim);
    }

    /**
     * Score two list of terms according to jc-max-avg-of-best method.
     *
     * @param tl1
     * @param tl2
     * @return
     */
    public double jcScoreMaxAvg(int[] tl1, int[] tl2)
    {
        return scoreMaxAvg(tl1, tl2, this.jcTermSim);
    }

    /**
     * Sim score avg one list of a term vs an item.
     *
     * @param tl1
     * @param item
     * @return
     */
    private double scoreMaxAvgVsItem(int[] tl1, int item, AbstractTermSim termSim)
    {
        if (termSim.maxScoreForItem != null) {
            double score = 0;
            for (int t1 : tl1) {
                score += termSim.maxScoreForItem[item][t1];
            }
            score /= tl1.length;
            return score;
        }

        return scoreMaxAvg(tl1, this.items2DirectTerms[item], termSim);
    }

    /**
     * Sim score avg one list of a term vs an item according to Resnik.
     *
     * @param tl1
     * @param item
     * @return
     */
    public double resScoreMaxAvgVsItem(int[] tl1, int item)
    {
        return scoreMaxAvgVsItem(tl1, item, this.resnikTermSim);
    }

    /**
     * Sim score avg one list of a term vs an item according to Lin.
     *
     * @param tl1
     * @param item
     * @return
     */
    public double linScoreMaxAvgVsItem(int[] tl1, int item)
    {
        return scoreMaxAvgVsItem(tl1, item, this.linTermSim);
    }

    /**
     * Sim score avg one list of a term vs an item according to Jiang and Conrath.
     *
     * @param tl1
     * @param item
     * @return
     */
    public double jcScoreMaxAvgVsItem(int[] tl1, int item)
    {
        return scoreMaxAvgVsItem(tl1, item, this.jcTermSim);
    }

    /**
     * Score one list of terms vs. an item using the default method and using the supplied term similarity measure.
     *
     * @param tl1
     * @param item
     * @param termSim
     * @return
     */
    private double scoreVsItem(int[] tl1, int item, AbstractTermSim termSim)
    {
        return scoreMaxAvgVsItem(tl1, item, termSim);
    }

    /**
     * Score one list of terms vs an item using the default method..
     *
     * @param tl1
     * @param item
     * @return
     */
    public double resScoreVsItem(int[] tl1, int item)
    {
        return scoreVsItem(tl1, item, this.resnikTermSim);
    }

    /**
     * Creates an array suitable for shuffling.
     *
     * @return
     */
    public int[] newShuffledTerms()
    {
        int[] shuffledTerms = new int[this.slimGraph.getNumberOfVertices()];

        /* Initialize shuffling */
        for (int i = 0; i < shuffledTerms.length; i++) {
            shuffledTerms[i] = i;
        }

        return shuffledTerms;
    }

    /**
     * Determined the pvalue of the given score for the given item.
     *
     * @param rnd defines the random number to be used
     * @param observedTerms
     * @param randomizedTerms
     * @param querySize
     * @param res
     * @param item
     * @param score
     * @param termSim
     * @return
     */
    private int simPValue(Random rnd, int[] observedTerms,
        int[] randomizedTerms, int querySize, Result res, int item,
        double score, AbstractTermSim termSim)
    {
        /* Turn it into a p value by considering the distribution */
        if (this.CACHE_RANDOM_QUERIES) {
            if (querySize > this.MAX_QUERY_SIZE_FOR_CACHED_DISTRIBUTION) {
                querySize = this.MAX_QUERY_SIZE_FOR_CACHED_DISTRIBUTION;
            }

            int[][] queries = getRandomizedQueries(rnd, querySize);

            if (this.CACHE_SCORE_DISTRIBUTION || this.PRECALCULATE_SCORE_DISTRIBUTION) {
                ApproximatedEmpiricalDistribution d = termSim.getScoreDistribution(querySize, item, queries);
                res.marginals[item] = 1 - (d.cdf(score, false) - d.prob(score));
            } else {
                int count = 0;

                for (int j = 0; j < this.SIZE_OF_SCORE_DISTRIBUTION; j++) {
                    double randomScore = scoreVsItem(queries[j], item, termSim);
                    if (randomScore >= score) {
                        count++;
                    }
                }

                res.marginals[item] = count / (double) this.SIZE_OF_SCORE_DISTRIBUTION;
            }
        } else {
            int count = 0;
            int[] shuffledTerms = newShuffledTerms();

            for (int j = 0; j < this.SIZE_OF_SCORE_DISTRIBUTION; j++) {
                chooseTerms(rnd, observedTerms.length, randomizedTerms, shuffledTerms);
                double randomScore = scoreVsItem(randomizedTerms, item, termSim);
                if (randomScore >= score) {
                    count++;
                }
            }
            res.marginals[item] = count / (double) this.SIZE_OF_SCORE_DISTRIBUTION;
        }
        return querySize;
    }

    /**
     * Makes the calculation according to a sim score avg max. We handle the observations as an item and compare it to
     * all other items. Also calculates the significance (stored in the marginal attribute).
     *
     * @param observations
     * @param pval to be set to true if significance should be determined
     * @param rnd the random source
     * @return
     */
    public Result simScore(boolean[] observations, boolean pval, AbstractTermSim termSim, Random rnd)
    {
        int[] observedTerms = getMostSpecificTermsSparse(observations);
        int[] randomizedTerms = new int[observedTerms.length];

        int querySize = observedTerms.length;

        Result res = new Result();
        res.scores = new double[this.allItemList.size()];
        res.marginals = new double[this.allItemList.size()];

        long startTime = System.currentTimeMillis();
        long lastTime = startTime;

        for (int i = 0; i < this.allItemList.size(); i++) {
            long time = System.currentTimeMillis();

            if (time - lastTime > 5000) {
                System.out
                    .println(termSim.name() + ": " + (time - startTime) + "ms " + i / (double) this.allItemList.size());
                lastTime = time;
            }

            /* Determine and remember the plain score */
            double score = scoreMaxAvgVsItem(observedTerms, i, termSim);
            res.scores[i] = score;

            querySize = simPValue(rnd, observedTerms, randomizedTerms, querySize, res, i, score, termSim);

        }

        return res;
    }

    /**
     * Makes the calculation according to Resnik avg max. We handle the observations as an item and compare it to all
     * other items. Also calculates the significance (stored in the marginal attribute).
     *
     * @param observations the input observations.
     * @param pval to be set to true if significance should be determined
     * @param rnd the random source
     * @return
     */
    public Result resnikScore(boolean[] observations, boolean pval, Random rnd)
    {
        return simScore(observations, pval, this.resnikTermSim, rnd);
    }

    /**
     * Makes the calculation according to Lin avg max. We handle the observations as an item and compare it to all other
     * items. Also calculates the significance (stored in the marginal attribute).
     *
     * @param observations the input observations.
     * @param pval to be set to true if significance should be determined
     * @param rnd the random source
     * @return
     */
    public Result linScore(boolean[] observations, boolean pval, Random rnd)
    {
        return simScore(observations, pval, this.linTermSim, rnd);
    }

    /**
     * Makes the calculation according to JC avg max. We handle the observations as an item and compare it to all other
     * items. Also calculates the significance (stored in the marginal attribute).
     *
     * @param observations the input observations.
     * @param pval to be set to true if significance should be determined
     * @param rnd the random source
     * @return
     */
    public Result jcScore(boolean[] observations, boolean pval, Random rnd)
    {
        return simScore(observations, pval, this.jcTermSim, rnd);
    }

    /**
     * Returns the term similarity according to Mathur and Dinakarpadnian.
     *
     * @param t1 term 1
     * @param t2 term 2
     * @return
     */
    public double mbTermSim(int t1, int t2)
    {
        return jaccard(t1, t2) * (this.terms2IC[t1] + this.terms2IC[t2]) / 2;
    }

    /**
     * Returns the msim according to Mathur and Dinakarpadnian, i.e., the maximum simimlarity between t1 and all of tl2.
     *
     * @param t1
     * @param tl2
     * @return
     */
    public double msim(int t1, int tl2[])
    {
        double s = 0.0;
        for (int element : tl2) {
            double snew = mbTermSim(t1, element);
            if (snew > 0.0) {
                s = snew;
            }
        }
        return s;
    }

    /**
     * Returns the unsymetric mbsim according to Mathur and Dinakarpadnian.
     *
     * @param tl1
     * @param tl2
     * @return
     */
    public double mbsimUnsym(int tl1[], int tl2[])
    {
        double s = 0.0;
        for (int element : tl1) {
            s += msim(element, tl2);
        }

        return s / tl1.length;
    }

    /**
     * Returns (symetric) mbsim according to Mathur and Dinakarpadnian.
     *
     * @param tl1
     * @param tl2
     * @return
     */
    public double mbsim(int tl1[], int tl2[])
    {
        return (mbsimUnsym(tl1, tl2) + mbsimUnsym(tl2, tl1)) / 2;
    }

    /**
     * Makes the calculation according to Mathur and Dinakarpadnian. We handle the observations as an item and compare
     * it to all other items.
     *
     * @param observations the input observations.
     * @param pval to be set to true if significance should be determined
     * @param rnd the random source
     * @return
     */
    public Result mbScore(boolean[] observations)
    {
        int[] observedTerms = getMostSpecificTermsSparse(observations);

        Result res = new Result();
        res.scores = new double[this.allItemList.size()];
        res.marginals = new double[this.allItemList.size()];

        long startTime = System.currentTimeMillis();
        long lastTime = startTime;

        for (int i = 0; i < this.allItemList.size(); i++) {
            long time = System.currentTimeMillis();

            if (time - lastTime > 5000) {
                System.out.println("mbScore: " + (time - startTime) + "ms " + i / (double) this.allItemList.size());
                lastTime = time;
            }

            /* Determine and remember the plain score */
            double score = mbsim(observedTerms, this.items2DirectTerms[i]);
            res.scores[i] = score;
        }
        return res;
    }

    /** Lock for randomized queries */
    private ReentrantReadWriteLock queriesLock = new ReentrantReadWriteLock();

    /**
     * Returns an array containing randomized term query. In the returned array, the first index distinguishes each
     * random query, and the second index distinguishes the terms.
     *
     * @param rnd source of random.
     * @param querySize defines the size of the query.
     * @return
     */
    private int[][] getRandomizedQueries(Random rnd, int querySize)
    {
        this.queriesLock.readLock().lock();
        int[][] queries = this.queryCache.getQueries(querySize);
        this.queriesLock.readLock().unlock();

        if (queries == null) {
            this.queriesLock.writeLock().lock();
            queries = this.queryCache.getQueries(querySize);
            if (queries == null) {
                int[] shuffledTerms = newShuffledTerms();

                queries = new int[this.SIZE_OF_SCORE_DISTRIBUTION][querySize];
                for (int j = 0; j < this.SIZE_OF_SCORE_DISTRIBUTION; j++) {
                    chooseTerms(rnd, querySize, queries[j], shuffledTerms);
                }

                this.queryCache.setQueries(querySize, queries);
            }
            this.queriesLock.writeLock().unlock();
        }
        return queries;
    }

    /**
     * Select size number of terms that are stored in chosen.
     *
     * @param rnd
     * @param size
     * @param chosen
     * @param storage
     */
    public void chooseTerms(Random rnd, int size, int[] chosen, int[] storage)
    {
        if (this.FORBID_ILLEGAL_QUERIES) {
            boolean valid;
            int tries = 0;

            do {
                choose(rnd, size, chosen, storage);
                valid = true;

                outer: for (int i = 0; i < size; i++) {
                    for (int j = 0; j < size; j++) {
                        if (i == j) {
                            continue;
                        }

                        /* If a chosen term is descendant of another one, we reject the query. */
                        if (this.slimGraph.isDescendant(chosen[i], chosen[j])) {
                            valid = false;
                            break outer;
                        }
                    }
                }
                tries++;
            } while (!valid && tries < 512);
        } else {
            choose(rnd, size, chosen, storage);
        }
    }

    /**
     * Chooses size randomly selected values from storage. Storage is manipulated by this call. Selected values are
     * stored in chosen.
     *
     * @param rnd
     * @param size number of elements that are chosen
     * @param chosen where the chosen values are deposited.
     * @param storage defines the elements from which to choose
     */
    public static void choose(Random rnd, int size, int[] chosen, int[] storage)
    {
        /*
         * Choose terms randomly as the size of observed terms. We avoid drawing the same term but alter shuffledTerms
         * such that it can be used again in the next iteration. Note that this duplicates code from the above.
         */
        for (int k = 0; k < size; k++) {
            int chosenIndex = rnd.nextInt(storage.length - k);
            int chosenTerm = storage[chosenIndex];

            /* Place last term at the position of the chosen term */
            storage[chosenIndex] = storage[storage.length - k - 1];

            /* Place chosen term at the last position */
            storage[storage.length - k - 1] = chosenTerm;

            chosen[k] = chosenTerm;
        }
    }

    /**
     * Returns the mica of term, i.e., the a common ancestor of the given terms whose information content is maximal.
     *
     * @param t1
     * @param t2
     * @return
     */
    public int getCommonAncestorWithMaxIC(int t1, int t2)
    {
        return commonAncestorWithMaxIC(t1, t2);
    }

    /**
     * Returns the current slim graph.
     *
     * @return
     */
    public SlimDirectedGraphView<Term> getSlimGraph()
    {
        return this.slimGraph;
    }

    /**
     * Returns the index of the given term.
     *
     * @param t
     * @return
     */
    public int getTermIndex(Term t)
    {
        return this.slimGraph.getVertexIndex(t);
    }

    /**
     * Returns the ontology.
     *
     * @return
     */
    public Ontology getOntology()
    {
        return this.graph;
    }

    /**
     * Returns the association container.
     *
     * @return
     */
    public AssociationContainer getAssociations()
    {
        return this.assoc;
    }

    /**
     * Returns the terms that are directly annotated to the given item.
     *
     * @param itemId
     * @return
     */
    public int[] getTermsDirectlyAnnotatedTo(int itemId)
    {
        return this.items2DirectTerms[itemId];
    }

    /**
     * Returns the frequencies of terms directly annotated to the given item. The order of the entries match the order
     * of getTermsDirectlyAnnotatedTo().
     *
     * @param itemId
     * @return
     */
    public double[] getFrequenciesOfTermsDirectlyAnnotatedTo(int itemId)
    {
        return this.items2TermFrequencies[itemId];
    }

    /**
     * Returns the parents of a given term.
     *
     * @param t
     * @return
     */
    public int[] getParents(int t)
    {
        return this.term2Parents[t];
    }

    /**
     * Returns whether the item has associated frequencies.
     *
     * @param item
     * @return
     */
    public boolean hasItemFrequencies(int item)
    {
        return this.itemHasFrequencies[item];
    }

    /**
     * Returns the evidence codes that should be respected. May be null in case all evidence codes are respected.
     *
     * @return
     */
    public String[] getEvidenceCodes()
    {
        return this.evidenceCodes;
    }

    /**
     * Returns the ic of the given term.
     *
     * @param t
     * @return
     */
    public double getTermIC(int t)
    {
        return this.terms2IC[t];
    }

    /**
     * Returns the number of items.
     *
     * @return
     */
    public int getNumberOfItems()
    {
        return this.allItemList.size();
    }
}
