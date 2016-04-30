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

#include "RandomNumberGenerator.hpp"
#include "CaVesselSegment.hpp"
#include "CaVessel.hpp"
#include "Owen2011SproutingRule.hpp"

template<unsigned DIM>
Owen2011SproutingRule<DIM>::Owen2011SproutingRule()
    : LatticeBasedSproutingRule<DIM>(),
      mHalfMaxVegf(0.1),
      mVegfField()
{

}

template <unsigned DIM>
boost::shared_ptr<Owen2011SproutingRule<DIM> > Owen2011SproutingRule<DIM>::Create()
{
    MAKE_PTR(Owen2011SproutingRule<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
Owen2011SproutingRule<DIM>::~Owen2011SproutingRule()
{

}

template<unsigned DIM>
void Owen2011SproutingRule<DIM>::SetHalfMaxVegf(double halfMaxVegf)
{
    mHalfMaxVegf = halfMaxVegf;
}

template<unsigned DIM>
std::vector<boost::shared_ptr<VascularNode<DIM> > > Owen2011SproutingRule<DIM>::GetSprouts(const std::vector<boost::shared_ptr<VascularNode<DIM> > >& rNodes)
{

    if(!this->mpGrid)
    {
        EXCEPTION("A regular grid is required for this type of sprouting rule.");
    }

    if(!this->mpVesselNetwork)
    {
        EXCEPTION("A vessel network is required for this type of sprouting rule.");
    }

    // Get the VEGF field
    this->mVegfField = this->mpSolver->GetPointSolution();

    // Set up the output sprouts vector
    std::vector<boost::shared_ptr<VascularNode<DIM> > > sprouts;

    // Loop over all nodes and randomly select sprouts
    for(unsigned idx = 0; idx < rNodes.size(); idx++)
    {
        if(rNodes[idx]->GetNumberOfSegments() != 2)
        {
            continue;
        }

        // Check we are not too close to the end of the vessel
        if(this->mVesselEndCutoff > 0.0)
        {
            if(rNodes[idx]->GetVesselSegment(0)->GetVessel()->GetClosestEndNodeDistance(rNodes[idx]->GetLocationVector())< this->mVesselEndCutoff)
            {
                continue;
            }
            if(rNodes[idx]->GetVesselSegment(1)->GetVessel()->GetClosestEndNodeDistance(rNodes[idx]->GetLocationVector())< this->mVesselEndCutoff)
            {
                continue;
            }
        }

        // Check we are not too close to an existing candidate
        if(this->mTipExclusionRadius>0.0)
        {
            bool too_close = false;
            for(unsigned jdx=0; jdx<sprouts.size(); jdx++)
            {
                if(rNodes[idx]->GetDistance(sprouts[jdx]->GetLocationVector()) < this->mTipExclusionRadius)
                {
                    too_close = true;
                }
            }
            if(too_close)
            {
                continue;
            }
        }

        // Get the grid index of the node
        unsigned grid_index = this->mpGrid->GetNearestGridIndex(rNodes[idx]->GetLocationVector());
        double vegf_conc = this->mVegfField[grid_index];
        double prob_tip_selection = this->mSproutingProbability*SimulationTime::Instance()->GetTimeStep()*vegf_conc/(vegf_conc + mHalfMaxVegf);

        if (RandomNumberGenerator::Instance()->ranf() < prob_tip_selection)
        {
            sprouts.push_back(rNodes[idx]);
        }
    }
    return sprouts;
}

// Explicit instantiation
template class Owen2011SproutingRule<2> ;
template class Owen2011SproutingRule<3> ;