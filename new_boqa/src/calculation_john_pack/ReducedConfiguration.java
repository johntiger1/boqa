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

package calculation_john_pack;

/**
 * Class to count the different cases.
 *
 * @author Sebastian Bauer
 */
final public class ReducedConfiguration implements Cloneable
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

    final public int[] stats = new int[ReducedConfiguration.NodeCase.values().length]; ///this is literally just 8 or so

//    @Override
//    public boolean equals(Object o)
//    {
//        if (!(o instanceof ReducedConfiguration))
//        {
//
//        }
//
//        return true;
//    }

    public int getCases(ReducedConfiguration.NodeCase c)
    {
        return this.stats[c.ordinal()]; //what if we used an array?!
    }

    public int getTotalNodes()
    {
        int total = 0;
        for (int count: stats)
        {
            total +=count;
        }

        return total;
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
        return Math.log(naive_beta) * getCases(ReducedConfiguration.NodeCase.FALSE_UNOBSERVED_NEGATIVE) +
                Math.log(experimental_beta) * getCases(NodeCase.FALSE_OBSERVED_NEGATIVE) +
                //True positives can only occur via being observed
                // pretty hard wired into the stats...
                //pretty dificult to inverse the operation (recover the false neg and TRUE false
                // neg from before--costly7
                //we dont even have their identities inside the stats! so it will actually be impossible
                //to check agains the existing observations.

                Math.log(alpha) * getCases(NodeCase.FALSE_POSITIVE) +
                Math.log(1 - experimental_beta) * getCases(NodeCase.TRUE_OBSERVED_POSITIVE) + //True positives can only occur via being observed
                //note that adding to one is satisfied: 1-experimental beta+e_beta =1
                //the misleading thing is that there is NO case where an obs can be true without having been examined via beta (e_beta)
                Math.log(1 - alpha) * getCases(NodeCase.TRUE_NEGATIVE);

//                +
//                Math.log(1) * getCases(NodeCase.INHERIT_FALSE) + /* 0 */
//                Math.log(1) * getCases(NodeCase.INHERIT_TRUE); /* 0 */
    }

    //Takes very little time ~395ns
    public void testSpeed(double alpha, double naive_beta, double experimental_beta)
    {
        long start = System.nanoTime();
        double sum = Math.log(1) * getCases(NodeCase.INHERIT_FALSE) + /* 0 */
                Math.log(1) * getCases(NodeCase.INHERIT_TRUE);
        //System.out.println("done looping through the 'find' op took " + (System.nanoTime()-start));
        System.out.println("done calculating two log0s, took" + (System.nanoTime()-start));

    }


    public void decrement(ReducedConfiguration.NodeCase c)
    {
        this.stats[c.ordinal()]--;
    }

    public void increment(ReducedConfiguration.NodeCase c)
    {
        this.stats[c.ordinal()]++;
    }

    @Override
    public String toString()
    {
        String str = "";
        for (int i = 0; i < this.stats.length; i++) {
            str += " " + ReducedConfiguration.NodeCase.values()[i].name() + ": " + this.stats[i] + "\n";
        }

        return str;
    }
    @Override
    public ReducedConfiguration clone()
    {
        ReducedConfiguration c = new ReducedConfiguration();
        for (int i = 0; i < this.stats.length; i++) {
            c.stats[i] = this.stats[i];
        }
        return c;
    }

    final public double falsePositiveRate()
    {
        return getCases(ReducedConfiguration.NodeCase.FALSE_POSITIVE)
                / (double) (getCases(ReducedConfiguration.NodeCase.FALSE_POSITIVE)
                + getCases(ReducedConfiguration.NodeCase.TRUE_NEGATIVE));
    }

    /**
     * Return false negative rate.
     *
     * @return
     */
    final public double falseNegativeRate()
    {
        int num_false_obs_neg = getCases(NodeCase.FALSE_OBSERVED_NEGATIVE);
        int num_false_unobs_neg = getCases(NodeCase.FALSE_UNOBSERVED_NEGATIVE);
        return (num_false_obs_neg + num_false_unobs_neg)
                / (double) (num_false_obs_neg + num_false_unobs_neg
                + getCases(NodeCase.TRUE_OBSERVED_POSITIVE));
    }

}
