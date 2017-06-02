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

import java.util.HashMap;

/**
 * Class to count the different cases.
 *
 * @author Sebastian Bauer
 */
final public class ReducedConfiguration  implements Cloneable
{
    public static enum NodeCase
    {
        FAULT,
        TRUE_OBSERVED_POSITIVE,
        FALSE_POSITIVE,
        TRUE_NEGATIVE,
        FALSE_UNOBSERVED_NEGATIVE,
        INHERIT_TRUE,
        INHERIT_FALSE,
        //TRUE_OBSERVED_POSITIVE, //IN FACT, we can just get rid of true positive with this
        FALSE_OBSERVED_NEGATIVE
    }
//
//    public HashMap<Integer,Boolean> registered_observations;
//
//    public ReducedConfiguration(HashMap<Integer,Boolean> reg_obs){
//        this.registered_observations = reg_obs;
//
//    }

    //internal record of the states of all the nodes @TODO optimize to just link it directly
    //to the ontology?
    final private int[] stats = new int[NodeCase.values().length]; ///this is literally just 8 or so
    //values wide!!!
    //J I finally get it: these are the BUCKETS from before, that we use to multiply everything
    final public void increment(NodeCase c)
    {
        this.stats[c.ordinal()]++;
    }

    final public void decrement(NodeCase c)
    {
        this.stats[c.ordinal()]--;
    }

    @Override
    public String toString()
    {
        String str = "";
        for (int i = 0; i < this.stats.length; i++) {
            str += " " + NodeCase.values()[i].name() + ": " + this.stats[i] + "\n";
        }

        return str;
    }

    /**
     * Get the number of observed cases for the given case.
     *
     * @param c
     * @return
     */
    final public int getCases(NodeCase c)
    {
        return this.stats[c.ordinal()];
    }

    /**
     * Returns the total number of cases that were tracked.
     *
     * @return
     */
    final public int getTotalCases()
    {
        int c = 0;
        for (int stat : this.stats) {
            c += stat;
        }
        return c;
    }

    /**
     * Returns the false positive rate.
     *
     * @return
     */
    final public double falsePositiveRate()
    {
        return getCases(ReducedConfiguration .NodeCase.FALSE_POSITIVE)
                / (double) (getCases(ReducedConfiguration .NodeCase.FALSE_POSITIVE)
                + getCases(ReducedConfiguration .NodeCase.TRUE_NEGATIVE));
    }

    /**
     * Return false negative rate.
     *
     * @return
     */
    final public double falseNegativeRate()
    {
        return getCases(ReducedConfiguration .NodeCase.FALSE_UNOBSERVED_NEGATIVE)
                / (double) (getCases(ReducedConfiguration .NodeCase.FALSE_UNOBSERVED_NEGATIVE)
                + getCases(ReducedConfiguration .NodeCase.TRUE_OBSERVED_POSITIVE));
    }

    /**
     * Returns the log score of the summarized ReducedConfiguration .
     * J: This is the famous 4 term expoenntiation, multiplication equation
     * @param alpha
     * @param naive_beta
     * @return
     */
    final public double getScore(double alpha, double naive_beta, double experimental_beta)
    {
        return Math.log(naive_beta) * getCases(NodeCase.FALSE_UNOBSERVED_NEGATIVE) +
                Math.log(experimental_beta) * getCases(NodeCase.FALSE_OBSERVED_NEGATIVE) +
                //True positives can only occur via being observed
                // pretty hard wired into the stats...
                //pretty dificult to inverse the operation (recover the false neg and TRUE false
                // neg from before--costly7
                //we dont even have their identities inside the stats! so it will actually be impossible
                //to check agains the existing observations.

                Math.log(alpha) * getCases(NodeCase.FALSE_POSITIVE) +
                Math.log(1 - experimental_beta) * getCases(NodeCase.TRUE_OBSERVED_POSITIVE) + //True positives can only occur via being observed
                Math.log(1 - alpha) * getCases(NodeCase.TRUE_NEGATIVE) +

                Math.log(1) * getCases(NodeCase.INHERIT_FALSE) + /* 0 */
                Math.log(1) * getCases(NodeCase.INHERIT_TRUE); /* 0 */
    }

    /**
     * Adds the given stat to this one.
     *
     * @param toAdd
     */
    final public void add(ReducedConfiguration  toAdd)
    {
        for (int i = 0; i < this.stats.length; i++) {
            this.stats[i] += toAdd.stats[i];
        }
    }

    /**
     * Clear the stats.
     */
    final public void clear()
    {
        for (int i = 0; i < this.stats.length; i++) {
            this.stats[i] = 0;
        }
    }

    @Override
    final public ReducedConfiguration  clone()
    {
        ReducedConfiguration  c = new ReducedConfiguration ();
        for (int i = 0; i < this.stats.length; i++) {
            c.stats[i] = this.stats[i];
        }
        return c;
    }

    public boolean equals(ReducedConfiguration  obj)
    {
        for (int i = 0; i < obj.stats.length; i++) {
            if (obj.stats[i] != this.stats[i]) {
                return false;
            }
        }
        return true;
    }
}
