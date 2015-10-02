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

#include "GeometryTools.hpp"
#include "OffLatticeTipAttractionGrowthDirectionModifier.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned DIM>
OffLatticeTipAttractionGrowthDirectionModifier<DIM>::OffLatticeTipAttractionGrowthDirectionModifier()
    : AbstractGrowthDirectionModifier<DIM>(),
      mpNetwork()
{

}

template <unsigned DIM>
boost::shared_ptr<OffLatticeTipAttractionGrowthDirectionModifier<DIM> > OffLatticeTipAttractionGrowthDirectionModifier<DIM>::Create()
{
    MAKE_PTR(OffLatticeTipAttractionGrowthDirectionModifier<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
OffLatticeTipAttractionGrowthDirectionModifier<DIM>::~OffLatticeTipAttractionGrowthDirectionModifier()
{

}

template<unsigned DIM>
void OffLatticeTipAttractionGrowthDirectionModifier<DIM>::SetNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

template<unsigned DIM>
c_vector<double, DIM> OffLatticeTipAttractionGrowthDirectionModifier<DIM>::GetGrowthDirection(c_vector<double, DIM> currentDirection,
                                                                                    boost::shared_ptr<VascularNode<DIM> > pNode)
{
    // Get the closest node in the search cone
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = mpNetwork->GetNodes();

    double min_distance = 1.e12;
    c_vector<double, DIM> min_direction = zero_vector<double>(DIM);
    for(unsigned idx=0; idx<nodes.size(); idx++)
    {
        if(IsPointInCone<3>(nodes[idx]->GetLocationVector(), pNode->GetLocationVector(), pNode->GetLocationVector() + currentDirection * 100.0, M_PI/3.0))
        {
            double distance = norm_2(pNode->GetLocationVector() - nodes[idx]->GetLocationVector());
            if(distance < min_distance)
            {
                min_distance = distance;
                min_direction = nodes[idx]->GetLocationVector() - pNode->GetLocationVector();
                min_direction /= norm_2(min_direction);
            }
        }
    }

    double strength;
    double crictical_distance = 100.0;
    if(min_distance >= crictical_distance)
    {
        strength = 0.0;
    }
    else
    {
        strength = 2.0 *  (1.0 - (min_distance * min_distance) / (crictical_distance * crictical_distance));
    }

    this->SetStrength(strength);
    return min_direction;
}

// Explicit instantiation
template class OffLatticeTipAttractionGrowthDirectionModifier<2> ;
template class OffLatticeTipAttractionGrowthDirectionModifier<3> ;
