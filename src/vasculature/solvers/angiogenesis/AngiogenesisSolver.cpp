/*

 Copyright (c) 2005-2015, University of Oxford.
 All rights reserved.

 University of Oxford means the Chancellor, Masters and Scholars of the
 University of Oxford, having an administrative office at Wellington
 Square, Oxford OX1 2JD, UK.

 This file is part of Chaste.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
 contributors may be used to endorse or promote products derived from this
 software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#include "UblasVectorInclude.hpp"
#include "VascularNode.hpp"
#include "AngiogenesisSolver.hpp"

template<unsigned DIM>
AngiogenesisSolver<DIM>::AngiogenesisSolver(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork) :
        mpNetwork(pNetwork),
        mGrowthVelocity(1.0),
        mTimeIncrement(1.0),
        mEndTime(100.0)
{

}

template<unsigned DIM>
AngiogenesisSolver<DIM>::~AngiogenesisSolver()
{

}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::Run()
{
    // Loop over the time (replace with simulation time)
    double current_time = 0.0;

    while(current_time < mEndTime)
    {
        current_time += mTimeIncrement;

        // Move any migrating nodes
        std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = mpNetwork->GetNodes();
        for(unsigned idx = 0; idx < nodes.size(); idx++)
        {
            if(nodes[idx]->IsMigrating() && nodes[idx]->GetNumberOfSegments()==1)
            {
                // Get the segment direction vector
                //c_vector<double,DIM> direction = nodes[idx]->GetLocationVector() -
            }
        }

    }

}

// Explicit instantiation
template class AngiogenesisSolver<2> ;
template class AngiogenesisSolver<3> ;

