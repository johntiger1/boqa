package calculation_john_pack;

import ontologizer.association.AssociationContainer;
import ontologizer.enumeration.GOTermEnumerator;
import ontologizer.enumeration.ItemEnumerator;
import ontologizer.go.Ontology;
import ontologizer.go.Term;
import ontologizer.go.TermID;
import ontologizer.set.PopulationSet;
import ontologizer.types.ByteString;
import sonumina.math.graph.SlimDirectedGraphView;

import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

/**
 * Created by johnchen on 18/05/17.
 */
public class ReducedBoqa {

    //Represents an array of disease distributions.
    //Each row represents a disease distribution, under a specific phenotype set to TRUE.
    //hence the dimensions of the array should be: # of phenotypes * # of diseases
    public double [][] multiDiseaseDistributions;

    public static int iteration = 0;
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

    public int getNumberOfUnobservedPhenotypes()
    {
        int count=0;

        //go thru the arraylist
        //
        //for (int i = 0; i < )
        //actually, this may be maintained in o
        //count = this.o.observations.length - this.o.real_observations.size();
        for (boolean b : this.o.observations)
        {
            count += b ? 0: 1;
        }

        return count;

    }


    public Ontology graph;
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


        SlimDirectedGraphView<Term> slimGraphView = this.graph.getSlimGraphView();
        this.term2Parents = slimGraphView.vertexParents; //recall this was an array of arrays
        this.term2Children = slimGraphView.vertexChildren;
        this.term2Ancestors = slimGraphView.vertexAncestors;
        this.term2Descendants = slimGraphView.vertexDescendants; //presumably, this is a vertex induced subgraph

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

        //
        multiDiseaseDistributions = new double[this.slimGraph.getNumberOfVertices()][this.allItemList.size()];
        System.out.println("num rows" + this.slimGraph.getNumberOfVertices());
        System.out.println("num cols" + this.allItemList.size());
        dfs();
    }

    //Returns the difference between the observed phenotypes, and the unobserved ones.
    //Note: this can be done in cosntant time, if it is done alongside the others


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
                //one reeason the following line works is because containsKey is a necessary condition for real_observations.get(chld)
                //the opposite case is if o.real_observations.get() is false, which is not handled anyways--this is not a False Positive!
                if (o.observations[chld] && o.real_observations[chld]) { //this explains false positive propagation
                    //translation: if both are observed has no real analogue anymore
                    //actually, if we actually do do true propagation, then
                    //we can just check that both are in
                    //however, this implies that we should have "truly observed"
                    //since, it could just be on solely because one of its kids was accidentally on


                    //if both are unobserved, it DOES have an analogue though! Both are not in
                    //the registered observations

                    if (o.observations[node] && (o.real_observations[node])) { //in the observed layer, parent points to child, so
                        //if child is on, then parenbt is on [contrary to how BN semantrics normally
                        //work, where really only parent affects child]
                        return ReducedConfiguration.NodeCase.INHERIT_TRUE; //the causality seems backwards
                        //it (the parent) is only on because its child was als oon
                        //however, this could have been handled immediately in a different graph
                        //preprocessing step

                        //we assume that it became true via inheritance
                    } else if (o.observations[node]  && (!o.real_observations[node])) {
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

                //note that before in the real_obs case, null represents FUN, while false represents FON
                //but consider a 1 0 (observed false) paired with a unobserved true...
                //is just checking that real_obs is false enough--yes, because it is default state...
                //unobs neg state leading to obs neg state though?
                //esnure matching on states...
                //i.e. nobs to unobs or obs to obs
                if (!o.observations[parent] || !o.real_observations[parent])  {
                    if (!o.observations[node] || !o.real_observations[node]) {
                        return ReducedConfiguration.NodeCase.INHERIT_FALSE; //wasn't actually false but inherited it..
                    } else {
                        /* NaN (the only other case is where it isnt null and it WAS observed!)*/
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
            if (o.observations[node] && o.real_observations[node]) {
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
                if (!o.observations[node])
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

            if (!o.observations[node] || !o.real_observations[node]) {

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

    //this "initializes the array" to a known default state, one where the target disease has NO hiddens.
    //Hence why the veyr first run of determineCasesForItem succeeds
    //Note also diffOff is also good here, since we literally have NO hiddens on so nothing to turn off
    //we most likely need something like this (initializing the STATS configuration to some known state)



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

    //returns the set difference of A and B
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
                c[cc++] = element; //if not in B, put the element in c and increment the c counter
            }
        } //copy everything over
        int[] nc = new int[cc];
        for (int i = 0; i < cc; i++) {
            nc[i] = c[i];
        }

        //retruns a newly allocated array
        return nc;
    }

    private static int[] setDiffSpecial(int[] a, int[] b, boolean [] skip)
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

            //only put it in if we are NOT skipping
            if (!skip[element]) {
                if (!inB) {
                    c[cc++] = element; //if not in B, put the element in c and increment the c counter
                }
            }
        } //copy everything over
        int[] nc = new int[cc];
        for (int i = 0; i < cc; i++) {
            nc[i] = c[i];
        }

        //retruns a newly allocated array
        return nc;
    }

    private void createPhenoVectors(boolean []skip) {

        //System.out.println("starting createPhenoVecs");
        //long start = System.nanoTime();
        int numberOfUnobservedPhenotypes = getNumberOfUnobservedPhenotypes();
        phenoOn = new int[numberOfUnobservedPhenotypes][];
        phenoOff = new int[numberOfUnobservedPhenotypes][];

        //this MEANS: we MUST turn off all of the phenotypes at the start!
                //this makes sense, since the one "behind" it is the empty set; hence all must be taken!


        //all terms annotated to the first item are diffOnTerms for that item as well
        int i ;
        int j = 0;
        int prev_j = 0;
        int sum =0;
        int starting = 0;

        //potentially skip many phenotypes
        while (o.observations[topo_sorted[j]])
        {
            j++;
            //TODO unlikely edge case: what if we observe everything?!
        }

        //assert j is total#-unobserved or something (not necessarily)
        //j is not really a descriptive statisitc, rather it only refers to the first, or the min
        prev_j = j;
        this.phenoOn[0] = this.term2Ancestors[topo_sorted[j++]]; /* For the first step, all terms must be activated */
        this.phenoOff[0] = new int[0];   //this is an empty array
        for (i = 1; i < numberOfUnobservedPhenotypes; i++) {
            //if unobserved ONLY

            //find the next unobserved phenotype, in the order given by topo_osrt


            while (o.observations[topo_sorted[j]] )
            {
                /*
                when would it ever fail:
                    if we run the outer loop too many times for example (note that j cannot decrement)
                 */
                if (j>=topo_sorted.length)
                {
                    //uhoh
                }
                j++;


            }


            //j must be from the previous incidence...

            int prevOnTerms[] = this.term2Ancestors[topo_sorted[prev_j]];      //items2terms[0] on first iteration
            int newOnTerms[] = this.term2Ancestors[topo_sorted[j]];
            prev_j = j;
            j++;
            this.phenoOn[i] = setDiffSpecial(newOnTerms, prevOnTerms, skip);
            this.phenoOff[i] = setDiffSpecial(prevOnTerms, newOnTerms, skip);
            sum += this.phenoOn[i].length + this.phenoOff[i].length;

        }

        System.out.println("number of deltas is " + sum);


        // System.out.println("done createPhenoVecs. Took" + (System.nanoTime()-start));



    }

    int [][][] phenoOnMT;
    int [][][] phenoOffMt;

    //Creates the phenovectors, partitioning them using the "linear" method.
    //

    //should use this.num_threads
    private void createPhenoVectors_multithread(int num_threads, boolean [] skip) {

//        System.out.println("starting createPhenoVecs");
//        long start = System.nanoTime();
        int numberOfUnobservedPhenotypes = getNumberOfUnobservedPhenotypes();
        int i ;
        int j = 0;
        int prev_j = 0;
        int sum =0;
        int starting = 0;
        prev_j = j;
        phenoOnMT = new int[num_threads][][];
        phenoOffMt = new int[num_threads][][];

        //assign inner sizes (task delegations) to each thread
        //z is the thread number

        int tasks_per_thread = numberOfUnobservedPhenotypes/num_threads;
        int remainder = numberOfUnobservedPhenotypes % num_threads;
        for (int z = 1; z < num_threads; z++)
        {
            phenoOnMT[z] = new int [tasks_per_thread][];
            phenoOffMt[z] = new int [tasks_per_thread][];
        }

        //assign any extra left over work to the main thread
        phenoOnMT[0] = new int [tasks_per_thread + remainder][];
        phenoOffMt[0] = new int[tasks_per_thread+ remainder][];


        int thread_num;
        /*
                #threads ,#unobserved phenotypes/#threads, deltas between then
        alternatively: we can do one of them, and then just pass the indices to different threads
        the only thing to care for is that every #phen/#threads, we set everything on specially
        --dangerous coupling and OBO possible

        the other obvious thing is to: delegate it to the calling thread
        Given that it is getting a slice of them, wipe it, then
        set all of the correct phenotypes ON--how long this will take?

        poor since reducing all the nodes is fine, but recomputing it for WHAT disease you are on
        and so forth is a little tricky



        dunkirk phantom of the opera similarity!

         */


        for ( thread_num =0; thread_num<num_threads;thread_num++ )
        {
            while (o.observations[topo_sorted[j]])
            {
                j++;
                //TODO unlikely edge case: what if we observe everything?!
            }

            //assert j is total#-unobserved or something (not necessarily)
            //j is not really a descriptive statisitc, rather it only refers to the first, or the min

            this.phenoOnMT[thread_num][0] = this.term2Ancestors[topo_sorted[j++]]; /* For the first step, all terms must be activated */
            //assert it gets a fresh copy (no phenotypes set) each time
            this.phenoOffMt[thread_num][0] = new int[0];

            for (i = 1; i < this.phenoOnMT[thread_num].length; i++) {
                //if unobserved ONLY

                //find the next unobserved phenotype, in the order given by topo_osrt


                while (o.observations[topo_sorted[j]] )
                {
                /*
                when would it ever fail:
                    if we run the outer loop too many times for example (note that j cannot decrement)
                 */
                    if (j>=topo_sorted.length)
                    {
                        //uhoh
                        System.out.println("wwoos");
                    }
                    j++;




                }


                //j must be from the previous incidence...

                int prevOnTerms[] = this.term2Ancestors[topo_sorted[prev_j]];      //items2terms[0] on first iteration
                int newOnTerms[] = this.term2Ancestors[topo_sorted[j]];
                prev_j = j;
                j++;

                //either need a special setdiff OR we can remove (prune) phenotypes that
                //are already observed

                //or we could just use more set differences!
                //subtract the set of phenotypes that already have been observed
//                phenoOnMT[thread_num][i] = setDiff(newOnTerms, prevOnTerms);
//                phenoOffMt[thread_num][i] = setDiff(prevOnTerms, newOnTerms);
                phenoOnMT[thread_num][i] = setDiffSpecial(newOnTerms, prevOnTerms,skip);
                phenoOffMt[thread_num][i] = setDiffSpecial(prevOnTerms, newOnTerms,skip);
                sum += phenoOnMT[thread_num][i].length + phenoOffMt[thread_num][i].length;

            }

        }

        //this MEANS: we MUST turn off all of the phenotypes at the start!
        //this makes sense, since the one "behind" it is the empty set; hence all must be taken!


        //all terms annotated to the first item are diffOnTerms for that item as well


        //potentially skip many phenotypes



        System.out.println("number of deltas is " + sum);


//         System.out.println("done createPhenoVecs (MT). Took" + (System.nanoTime()-start));
    }


    private void createDiffVectors()
    {

        System.out.println("starting createDiffVecs");
        long start = System.nanoTime();


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
            //jagged arrays autodone!
            this.diffOnTerms[i] = setDiff(newOnTerms, prevOnTerms);
            this.diffOffTerms[i] = setDiff(prevOnTerms, newOnTerms);

            sum += this.diffOnTerms[i].length + this.diffOffTerms[i].length;
        }

        System.out.println("done createDiffVecs. Took" + (System.nanoTime()-start));


    }

    //provides TWO arrays such that phenoOn specifies the nodes which must
    //be switched on vs. the previous phenotype, while diffOff specifies which
    //nodes must be switched OFf vs. the previous phenotype
    //uses a HashSset based approach
    int [][] phenoOn; //A[i] specifies how to get from A[i] to A[i+1]
    int [][] phenoOff;
    HashSet [] pOn;
    HashSet [] pOff;
    private void computePhenoDifferentials()
    {
        //between n numbers there are only n-1 spaces
        pOn = new HashSet[getNumberOfUnobservedPhenotypes()];
        pOff = new HashSet[getNumberOfUnobservedPhenotypes()];

        //each hashset only stores the delta from one pheno to the next->10k at the most
        HashSet hOn = new HashSet();
        HashSet hOff = new HashSet();

        //since we must compute n different duplicates, we do in fact need
        //double for loop


        //if duplicates are encountered skip
        //else, add it to the OTHER hashset

        //after done BOTH for loops we have computed the diff between
        //what i am actually computing now is the term with its next term differences
        //
        int z;
        for (int i =0; i < this.slimGraph.getNumberOfVertices(); i++) {

            if (!this.o.observations[i]) {

                //Add in all the ancestors of the pheno (includes itself)
                //Since we only store the differentials between them, we cannot use
                //phenoOn and phenoOff directly
                //(without applying it many times)
                for (int k = 0; k < term2Ancestors[i].length; k++) {
                    hOn.add(term2Ancestors[i][k]);
                }

                //goes over the term2ancestors of the next term
                //this should compute the delta between n and 0
                //essentially, i+1%x describes a function that everywhere is i, except at the very boundary, in the range [0, n-1]

                //LOOKS AT THE NEXT (wrapping around if necessary) term's things
                //j iterates over the terms of the next phenotype
                for (int j = 0; j < term2Ancestors[z=(i + 1) % pOn.length].length; j++) {


                    if (!hOn.contains(term2Ancestors[z][j])) {
                        hOff.add(term2Ancestors[z][j]);
                        //these are definitely in B
                        //the ones in A - the duplicates are the ones exclusive to just A
                    } else //previous hashset and this one both contain the term
                    {
                        hOn.remove(term2Ancestors[z][j]);
                        //now at this point, h will only have items which are exclusive ot itself
                        //then:
                    }

                    //we need a case of A-B => this implies what to TURN OFF
                    //(actually, they have a neat trick: )
                    //it just is to run it in the reverse way!
                }

                //what must be turned on from i to i+1
                pOn[i] = (HashSet) hOff.clone();
                //what must be turned off from i to i+1
                pOff[i] = (HashSet) hOn.clone();

                hOn.clear();
                hOff.clear();


            }
        }


    }

    //Takes in a stats, and implictly needs the diffOn to be created
    //it shall update the stats with the result.

    //however, the new setting on will not score anything, instead it will
    //
    //actually it may be impossible to decouple these
    //it isn't. what we must do however is pass out the
    //actually: the diffOnTerms is ALREADY OK
    //i.e. it is already sufficient to specify the deltas/differentials
    //between one diseaase to the next
    //hence just replace the code below with a for loop that loops thru the
    //phenotypes
    //rename this: just return a multiarray of all the different disease probabiloitiews
    //under the different phenotypes

    //This will need to be a multifaceted, all-encompassing function!
    //(i.e. the way it is written, it is difficult to encapsulate/decouple)


    private void actuateDiseaseDifferentials( int item, Observations o,
                                             boolean takeFrequenciesIntoAccount, boolean[] previousHidden, ReducedConfiguration previousStats)
    {

        int numTerms = this.slimGraph.getNumberOfVertices();

        if (previousHidden == null && previousStats != null) {
            throw new IllegalArgumentException();
        }
        if (previousHidden != null && previousStats == null) {
            throw new IllegalArgumentException();
        }

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

        int[] diffOn = this.diffOnTerms[item];
        int[] diffOff = this.diffOffTerms[item];

        useDeltaToSetToCurrentDisease(hidden, diffOn, diffOff, stats, o);

        //assert: pOff and pOn should have the same length (they represent the # of transitions between unobserved phenotypes)


        int k;
        //the last one specifies how to get back to n-1
        //System.out.println("starting looping thru all phens");
        long start = System.nanoTime();


        //iterate up to the last one (take special action on the VERY last)
        //its ok since it depends only on the # of unobserved phenotpes, which is phenoOn.length
        for (int j =0; j<phenoOn.length; j++) {
            //dangerous dependency: now we have a monotone list of numbers, but we cannot retrieve their original indices
            //System.out.println("done");
            //start = System.nanoTime();

            //This will be problematic if we are no longer sequential..
            //k refers to the current. Note that it is possible we cannot find it!
            //in this case, we should TERMINATE

            //POTENTIALLY quite slow--due to various array "fetch" operations!
            //This is dangeroulsy coupled to phenoOn (especially the length)
            //Essentially, we are asserting the following invariant:
            //the number of unobserved phenotypes is equal to the number of observations false in the array
            //however, consider whaat happens if k >j
            //this seems to be overloaded, somewhat...
            //it all boils down to the fact that phenOn [-1] and pOn[-1 ] was never convincingly defined
            //also: we must iterate over ALL phenotypes, i.e. we must search thru the ENTIRE topo_sorted array
            //since we may observe at any point and or time
            //at the same time we must keep it synchronized with phenoOn.length as well
            //it would make a lot of sense to maintain an sorted arraylist of unobserved phenotypes somewhere..
            /*
            We look through the topo_sorted array for uncheck phenotypes.
            In HIGH theory: we should have no transitions computed for phenotypes that are already oberved
            This preserves the correctness of the decrementing and so forth

            K< topo_sorted.length!! is the key thing
            new reduction: we only need topo_sorted, and just ensure that it is truncated appropriately
            what the below is trying to do is to scan over all phenotyes again!
             */

            //we DONT need this line at all, if everything is ACTUALLY set up correctly
            //however, the lnegth depends on the
            //i.e. we will always ensure that k=j, since obs[topo_sorted[j]] is not true

            //this line works if
            //for(k=j;k < o.observations[topo_sorted[k]];k++);
            //also inefficient
            //there are a couple of problems here
            //we skip when necessary and so forth...


            //go thru all the phenotypes to search for an unobserved one
//            while (k < phenoOn.length && o.observations[topo_sorted[k]])
//            {
//
//
//
//            }
            //this finds the first index past j where observations is NOT true
            //ensure that the array is not all observed!


            //array op--this is very rapid.
            //System.out.println("done looping through the 'find' op took " + (System.nanoTime()-start));

            //these js will be used in pOn
            decrementStaleNodes(o, hidden, stats, j);
            updateObservationsBasedOnPhenotype(o, j);


            //Note that while the function expects an item (hidden node) to change, here we are actually changing
            //the observed node!

            //true: we should be able to pass this off to a function, specifying the action occurring.
            //pass in the iterator, and the action to take n stats, and so forth

            incrementUpdatedNodes(o, hidden, stats, j);


            //fill in the COLUMN of the matrix!
            //this will need to be topo_sorted[k] now...
            multiDiseaseDistributions[topo_sorted[j]][item] = stats.getScore(this.ALPHA_GRID[0],initial_beta, experimental_beta);
        }

        //reset to the original/baseline phenotype


        resetToConsistentStateForDiseaseDiff(o, stats, diffOn, diffOff, hidden);

        //System.out.println("done looping through all the phenos. Took" + (System.nanoTime()-start));
    }

    int curr_thread = 0;

    //move the for loop (other the indices) into here
    private void actuateDiseaseDifferentials_MT( int start_disease, int end_disease,
                                                 int start_pheno, int end_pheno,
                                                 Observations o,
                                              boolean takeFrequenciesIntoAccount,
                                                 boolean[] previousHidden,
                                                 ReducedConfiguration previousStats,
                                                 int thread_id)
    {


        int numTerms = this.slimGraph.getNumberOfVertices();

        if (previousHidden == null && previousStats != null) {
            throw new IllegalArgumentException();
        }
        if (previousHidden != null && previousStats == null) {
            throw new IllegalArgumentException();
        }

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
        int[] diffOn ;
        int[] diffOff ;
        int pheno_counter = start_pheno;

//        System.out.printf(" I am %d, i am responsibl" +
//                        "e for phenos %d to %d and diseases %d to %d\n",thread_id, start_pheno
//        , end_pheno, start_disease, end_disease);
//        System.out.println("i am " + thread_id + " i have stats at " + stats.toString());
        //        System.out.println("predifferences" + (end_pheno-start_pheno) );
//        System.out.printf("Going from this disease %d to this: %d", start_disease, end_disease);
//        System.out.printf("Going from this pheno %d to this: %d", start_pheno, end_pheno);
        int counter = 0;
        for (int item = start_disease; item < end_disease; item++)
        {

            diffOn = this.diffOnTerms[item];
            diffOff = this.diffOffTerms[item];

            useDeltaToSetToCurrentDisease(hidden, diffOn, diffOff, stats, o);
            long start = System.nanoTime();
//            System.out.println("differences" + (end_pheno-start_pheno) );
//            System.out.println("length" + phenoOnMT[thread_id].length);
//            assert end_pheno-start_pheno == phenoOnMT[thread_id].length;
            //for some, j should be multiplied by thread_id+1 *j
            //it's the solution to the modulus quesiton
            //thread_id*j + j
            //first thread does follow this, but all other threads: OFFSET + thread_id*j

            for (int j =0; j<phenoOnMT[thread_id].length; j++) {

                if (item==0 && j==0){
                System.out.println("I am thread" + thread_id + " on iteration " + iteration +
                        "my pON[0] is"
                + Arrays.toString(phenoOnMT[thread_id][0]));}

                int freeze_disease = 230;
                int freeze_pheno = 203;
                if (item == freeze_disease && j == freeze_pheno)
                {
                    System.out.println("Before decrement I am thread" + thread_id + "stats to be swcored is" + stats);
                }
                counter+= decrementStaleNodes_MT(o, hidden, stats, phenoOnMT[thread_id][j],phenoOffMt[thread_id][j]);
                if (item == freeze_disease && j == freeze_pheno)
                {
                    System.out.println("After decrement I am thread" + thread_id + "stats to be swcored is" + stats);
                }
                updateObservationsBasedOnPhenotype_MT(o, phenoOnMT[thread_id][j], phenoOffMt[thread_id][j]);

                incrementUpdatedNodes_MT(o, hidden, stats, phenoOnMT[thread_id][j], phenoOffMt[thread_id][j]);
                    //this means phenocounter is writing its value down the col to too many elements
                if (item == freeze_disease && j == freeze_pheno)
                {
                    System.out.println("After increment" + thread_id + "stats to be swcored is" + stats);
                }


                multiDiseaseDistributions[topo_sorted[pheno_counter++]][item] = stats.getScore(this.ALPHA_GRID[0],initial_beta, experimental_beta);

                //assert start_pheno < end_pheno
                    //implicit dependence on phenoOnMT.length being right


            }
            pheno_counter = start_pheno;
            resetToConsistentStateForDiseaseDiff(o, stats, diffOn, diffOff, hidden);
            //should require an incrmeent though?!
        }

//        System.out.println("i am finishing " + thread_id + " i have stats at " + stats
//        + "also this is my decs" + counter);

//        System.out.println("done my job" + thread_id);



    }

    int [] topo_sorted;
    /*
      Performs an after-the-fact dfs of a "graph"
     for each of the descendants of the root term, get their descendants and so forth
    this produces a DFS
     */
    public void dfs()
    {

//        System.out.println("This is the root term" + graph.getRootTerm());
//        System.out.println("These are the children");
//        long start = System.nanoTime();

        ArrayList<Term> topo_sort = graph.getTermsInTopologicalOrder();
        //now: we need the topo sort
        //convert it to the indices.
        topo_sorted = new int[topo_sort.size()];

        int i = 0;
        //start = System.nanoTime();
        for (Term t: topo_sort)
        {
            //assume this is the same as a[i], then i++
            //topo_sorted[i++] = slimGraph.getVertexIndex(t);
            //use i, then increment it

            topo_sorted[i] = slimGraph.getVertexIndex(t);
            i++;
        }
        //System.out.println("filling new array took " + (start-System.nanoTime()));

        //if we compute via this array, then we will have GOOD things!



    }
    private void incrementUpdatedNodes_MT(Observations o, boolean[] hidden, ReducedConfiguration stats,
                                       int[]phenoOnSlice,int[]phenoOffSlice) {

        for (int x : phenoOnSlice)
        {
            stats.increment(getNodeCase(x, hidden, o)); //decrements an arbitrary one of the same class!
        }

        for (int x : phenoOffSlice)
        {
            stats.increment(getNodeCase(x, hidden, o)); //decrements an arbitrary one of the same class!
        }
    }

    private void incrementUpdatedNodes(Observations o, boolean[] hidden, ReducedConfiguration stats, int j) {
        for (int x : phenoOn[j])
        {
            stats.increment(getNodeCase(x, hidden, o)); //decrements an arbitrary one of the same class!
        }

        for (int x : phenoOff[j])
        {
            stats.increment(getNodeCase(x, hidden, o)); //decrements an arbitrary one of the same class!
        }



//        for (Object val : pOn[j])
//        {
//            stats.increment(getNodeCase((int) val, hidden, o));
//        }
//
//        //re-evalaute the ones that were changed
//        for (Object val : pOff[j])
//        {
//            stats.increment(getNodeCase((int) val, hidden, o));
//        }
    }

    private void updateObservationsBasedOnPhenotype(Observations o, int j) {
        //ojbects that implement the iterator can sue the for each
//        for (Object val : pOff[j])
//        {
//            o.real_observations[((int) val)] = false;
//        }
//
//        for (Object val : pOn[j])
//        {
//            //TODO think about how we can infer duplicates/etc. using this system
//            //(since the record obs already stores the stae of the previous observation)
//            //we do not actually need the indexing terms2ancestors[i][j], since that is already completely
//            //reflected in the state of the o
//            o.real_observations[((int) val)] = true;
//
//        }

        for (int x : phenoOff[j])
        {
            //this is the same as the array[array[2]] = true semantics
            o.real_observations[x] = false; //decrements an arbitrary one of the same class!
        }



        for (int x : phenoOn[j])
        {
            o.real_observations[x] = true; //decrements an arbitrary one of the same class!
        }
    }

    private void updateObservationsBasedOnPhenotype_MT(Observations o,int [] phenoOnSlice, int [] phenoOffSlice
                                                    ) {
        for (int x : phenoOffSlice)
        {
            //this is the same as the array[array[2]] = true semantics
            o.real_observations[x] = false; //decrements an arbitrary one of the same class!
        }
        for (int x : phenoOnSlice)
        {
            o.real_observations[x] = true; //decrements an arbitrary one of the same class!
        }
    }


    //Decrement the nodes that are "stale" from the previous to our current
    //essentially, changes everything except the nodes which remain exactly the same
    private void decrementStaleNodes(Observations o, boolean[] hidden, ReducedConfiguration stats, int j) {
//        Iterator i = pOn[j].iterator();
//        while (i.hasNext()) {
//            stats.decrement(getNodeCase((int) (i.next()), hidden, o));
//        }

        //        for (i = pOff[j].iterator(); i.hasNext();)
//        {
//            stats.decrement(getNodeCase((int) i.next(), hidden, o));
//        }
        for (int x : phenoOn[j])
        {
            stats.decrement(getNodeCase(x, hidden, o)); //decrements an arbitrary one of the same class!
        }

        for (int x : phenoOff[j])
        {
            stats.decrement(getNodeCase(x, hidden, o));
        }


    }

    private int decrementStaleNodes_MT(Observations o, boolean[] hidden, ReducedConfiguration stats,
                                        int [] phenoOnSlice, int [] phenoOffSlice) {

        int counter = 0;
        for (int x : phenoOnSlice)
        {
            stats.decrement(getNodeCase(x, hidden, o)); //decrements an arbitrary one of the same class!
            counter++;
        }

        for (int x : phenoOffSlice)
        {
            stats.decrement(getNodeCase(x, hidden, o));
            counter++;
        } ///should also incrment here??
//        System.out.println("dec: " + counter);
        return counter;
    }

    //Resets the d-p combo into a "consistent" state, such that the NEXT diffON
    //can be applied.
    private void resetToConsistentStateForDiseaseDiff(Observations o, ReducedConfiguration stats, int[] diffOn, int[] diffOff
    , boolean []hidden) {

        List <Integer> temp = new ArrayList<>();
        //undo all the changes we have done via phenoOn and phenoOff
        for (int i = 0; i < o.observations.length; i++)
        {
            //unobserved must be false. this is OK since when we infer a node, we set it to True (observed) as well
            if (!o.observations[i] && o.real_observations[i])
            {
                o.real_observations[i] = false;
                temp.add(i);
            }

        }
        //this is problematic
        //it MIGHT be that this is : "resetting it to the same disease each time"
        //else this is a semi-convoluted way of resetting
        //try doing it from the real obs that were reset
//        for (int element : diffOn) {
//            stats.decrement(getNodeCase(element, hidden, o));
//        }
//        for (int element : diffOff) {
//            stats.decrement(getNodeCase(element, hidden, o)); //lookup the hidden and observed too
//        }

        for (int x : temp)
        {
            stats.decrement(getNodeCase(x, hidden, o));
        }

        //TODO it would make more sense to do the decrements here (since it reoresents things that have changed)


//        for (int x: phenoOff[phenoOff.length-1])
//        {
//            o.real_observations[x] = false;
//        }
//
//        for (int x: phenoOn[phenoOn.length-1])
//        {
//            o.real_observations[x] = true;
//        }

//        for (Object val : pOff[pOff.length-1])
//        {
//            o.real_observations[((int) val)] = false;
//        }
//
//        for (Object val : pOn[pOn.length-1])
//        {
//            o.real_observations[((int) val)] = true;
//        }
    }


    private void useDeltaToSetToCurrentDisease( boolean[] hidden,  int[] diffOn, int[] diffOff,
                                                ReducedConfiguration stats, Observations obs) {


            /* Change nodes states */ //why is hidden[0] always on, esp. when we set observed to be on only
        for (int i = 0; i < diffOn.length; i++) {
            hidden[diffOn[i]] = true; //each element in the diffOn array, will change
            /*
            //all the changes that will be made to i for the current context will already have been made
            in particular, i cannot appear below in diffOff!
             */
            stats.increment(getNodeCase(diffOn[i], hidden, obs));

            //(regarde, diffon is fixed for the item)
        }
        for (int i = 0; i < diffOff.length; i++) {
            hidden[diffOff[i]] = false; //what we need to switch off to get to disease i
            stats.increment(getNodeCase(diffOff[i], hidden, obs));
        }



    }

    private WeightedConfigurationList determineCasesForItem(int item, Observations o,
                                                            boolean takeFrequenciesIntoAccount, boolean[] previousHidden, ReducedConfiguration previousStats)
    {
        //TODO just pass in this field instead of computing it each time
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
            /* New */ //note that at 0, many thing smust be turned on (since it is the first one)
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

    public Result assignMarginals( Observations observations, final boolean takeFrequenciesIntoAccount,
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

        //presumably these are problematic. They should be updated, to whatever they were/should be
        //They were FINE in the previous iteration, because in determine cases we did have an initializing case.
        //i.e. diffOn[0] specified the entire disease to put on,
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
//                    for (WeightedConfiguration wc : stats) //this is a list of weighted configs, which itself is simply a list of
//                    {
//                        System.out.println("Itme " + item);
//                        System.out.println(wc.stat);
//
//
//
//                    }
                    //since multithreading with the differentials would be too difficult,
                    //we only sequentially set flags in the stats, (when we only have one process)

                    //PERHAPS: we should assign stats to some of the res.
                    //In general, res.stats is not updated again, sadly

                    //System.out.println(stats.toString());
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

                    //System.out.println(observations.observationStats != null);
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

    int num_threads = 1;

    //ideally it should be always an array, with it being size 1 * n for a single thing
    int [][]phenoOff_array=null;
    //Returns an array of Result objects
    //will now be quite an expensive operation
    public void assignMultiShotMarginals( Observations observations, final boolean takeFrequenciesIntoAccount,
                                   final int numThreads)
    {
        this.num_threads = numThreads;

        System.out.println("starting compute pheno diff");
        long start = System.nanoTime();
        //computePhenoDifferentials();
//        createPhenoVectors(observations.observations);//must be recomputed on every call
        createPhenoVectors_multithread(numThreads, observations.observations);

        System.out.println("done compute pheno diff. Took" + (System.nanoTime()-start));

        int i;

        final Result res = new Result();
        int num_diseases = this.allItemList.size();
        res.scores = new double[num_diseases]; //its an array size of the entire BN
        res.marginals = new double[num_diseases];
        res.marginalsIdeal = new double[num_diseases];
        res.stats = new ReducedConfiguration[num_diseases]; // an ARRAY of ReducedConfigurations..
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
        final double[][][] scores = new double[num_diseases][this.ALPHA_GRID.length][this.BETA_GRID.length];
        final double[] idealScores = new double[num_diseases]; //this is likely a "target" whom we wish to approximate
        //as much as possible
        //interestingly we also have ideal marginals, which I am not sure how the two differ by

        boolean [] all_completed = new boolean[num_diseases];

        final ExecutorService es;
        es = Executors.newFixedThreadPool(num_threads);
            //use newCached for performance
            /*
            actually i claim that fixed is better, since we broke it up into that way
             */

        final boolean[] previousHidden = new boolean[this.slimGraph.getNumberOfVertices()];
        final ReducedConfiguration previousStat = new ReducedConfiguration();

        //Initialize the previousStat config(no diseases, but with current observations).
        // This step is necessary for the createDiffVectors to make sense.
        determineCases(observations, previousHidden, previousStat);

        int curr_start_indexP = 0;

        int curr_end_indexP= 0;
        int num_diseases_per_thread = num_diseases/num_threads;
//        System.out.println("hey this is " + num_diseases_per_thread );
        int remainder = num_diseases%num_threads;
        //actualy it can take both more diseases AND more pehnotypes
        //note that we cannot use the sme info from the phenoOn, as we may have different cardinality P and D
        int curr_start_indexD = 0;
        int curr_end_indexD= num_diseases;
        for (i = 0; i < num_threads; i++) {

            curr_start_indexP = curr_end_indexP;
//                System.out.println(curr_end_indexP);
            curr_end_indexP+=phenoOffMt[i].length; //uses the NEXT one

//            if (i != 0) {curr_start_indexD = curr_end_indexD;}
//            curr_end_indexD += num_diseases_per_thread;
            //this specifies how many phenotypes we are going over

            //later an es.execute() call
            class MatrixWorker implements Runnable{
                public int id;
                public Observations o;
                public int start_item;
                public int end_item;
                public int pheno_start;
                public int pheno_end;
                public int num_elements;
                public boolean []hidden;
                public ReducedConfiguration rc;
                public MatrixWorker(int id, Observations o, int start_item, int end_item, int pheno_start,
                                    int pheno_end, boolean[]hidden, ReducedConfiguration rc)
                {
                    this.id = id;
                    this.o=o;
                    this.start_item=start_item;
                    this.end_item=end_item;
                    this.pheno_start=pheno_start;
                    this.pheno_end = pheno_end;
                    this.hidden = Arrays.copyOf(hidden, hidden.length);
                    //also need to copy previousStat and so forth
                    this.rc=rc.clone(); //recall observations has a config inside it!
                }
                public void run()
                {
                    actuateDiseaseDifferentials_MT(start_item,end_item,pheno_start, pheno_end,o, takeFrequenciesIntoAccount,
                            numThreads > 1 ? this.hidden : previousHidden, numThreads > 1 ? this.rc: previousStat, id);



                }

            }
            //let 8 processes do the job of actuateDiseaseDiff
            //using an executor service would probably be best here (thread destruction/creation)

            //just submit the tasks for running here
            //
            //es.execute

//                System.out.println(o.observationStats);
                MatrixWorker m = new MatrixWorker(i, new Observations(o), curr_start_indexD,
                        curr_end_indexD, curr_start_indexP, curr_end_indexP, previousHidden, previousStat);
//            System.out.println("pheno start is " + curr_start_indexP);
//            System.out.println("pheno ebd is " + curr_end_indexP);
                es.execute(m);



        }

            es.shutdown();

            try {
                es.awaitTermination(10, TimeUnit.SECONDS);

//                for (int a = 0; a < all_completed.length; a++)
//                {
//                    if (!all_completed[a])
//                    {
//                        System.out.println("OH NO!");
//                    }
//                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        // System.out.println(idealNormalization + " " + normalization);
        // if (exitNow)
        // System.exit(10);
    }


