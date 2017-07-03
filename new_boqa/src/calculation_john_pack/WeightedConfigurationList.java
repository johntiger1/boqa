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

import java.util.ArrayList;
import java.util.Iterator;

/**
 * This is a list of weighted stats.
 *
 * @author Sebastian Bauer
 */
public class WeightedConfigurationList implements Iterable<WeightedConfiguration>
{
    private ArrayList<WeightedConfiguration> tupelList = new ArrayList<WeightedConfiguration>(10); //TODO typo and hardcoded
    //interestingly this is actually what default araylist size is anyways
    ArrayList<String> jim = new ArrayList<>();
    //Why we should remove it:
    //preusumably this is the top 10 diseases

    @Override
    public Iterator<WeightedConfiguration> iterator()
    {
        return this.tupelList.iterator();
    }

    public void add(ReducedConfiguration stat, double factor)
    {
        WeightedConfiguration t = new WeightedConfiguration();
        t.stat = stat;
        t.factor = factor;
        this.tupelList.add(t);
    }

    public double score(double alpha, double e_beta, double n_beta)
    {
        double sumOfScores = Math.log(0);
        //System.out.println("my value is " + sumOfScores);

        for (WeightedConfiguration tupel : this.tupelList) {
            double score = tupel.stat.getScore(alpha, n_beta, e_beta) + tupel.factor; /*
                                                                             * Multiply score by factor, remember that
                                                                             * we are operating in log space
                                                                             */
            sumOfScores = Util.logAdd(sumOfScores, score);
            //System.out.println("my new value is " + sumOfScores); //somehow, it becomes NON NEGATIVE INFINITY AFTER
        }
        return sumOfScores;
    }

    public int size()
    {
        return this.tupelList.size();
    }

}
