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

import java.util.HashMap;

/**
 * Class maintaining observations and stats regarding to true positives.
 *
 * @author Sebastian Bauer
 */
public class Observations
{
    public int item;

    //this is the way it works:
    //default state is observations= [0,..,0]
    //whenever we make an observation, we flip the corresponding elt to 1, AND set the registeredobservation
    //the registeredobservation contains whether or not it was "truly" set to 0/1
    //we msut take great care to only check registered observations if the corresondonding observations is 1
    //(otherwise it has no meaning)
    //dangerous parallel array, just use an object!

    //perhaps just make it an object itself!

    //it makes most sense just to have an arraylist of <item: observation_value> pairings
    //In general, we can make these checks by asserting abscenjnece in real_observations
    @Deprecated
    public boolean[] observations;

    //public boolean [] RegisteredObservaitons; // records whether the observation was true or false
    public Observations()
    {
        //observations should be set when necessary (fixed size)
    }
    //public ArrayList<obs_value> real_observations= new ArrayList<>();
    public boolean[] real_observations;
    public ReducedConfiguration observationStats; //this is our familiar 4 term!

//    //TODO UNSAFE METHOD
//    public void recordObservation(int index, boolean value)
//    {
//        observations[index] = true;
//        RegisteredObservaitons[index] = value;
//
//    }
}
