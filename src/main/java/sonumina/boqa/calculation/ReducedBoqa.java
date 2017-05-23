package sonumina.boqa.calculation;

import ontologizer.association.AssociationContainer;
import ontologizer.go.Ontology;

import java.util.ArrayList;
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
    public void setup(Ontology ontology, AssociationContainer assocs)
    {
        this.assoc = assocs;
        this.graph = ontology;
        //without graph inducing--how long will the computation take?
        //why is inducement a valid thing?
        //is it because of the only one valid configuration for a disease type thing ? => 0's in the Hidden
        //still mess it up since they could be 1 in the observed (False Positive). Hence
        //not entirely sure of rational to break it up
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
}
