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
#include "AbstractSproutingRule.hpp"

template<unsigned DIM>
AbstractSproutingRule<DIM>::AbstractSproutingRule()
    : mSproutingProbability(0.1),
      mNodes()
{

}

template<unsigned DIM>
AbstractSproutingRule<DIM>::~AbstractSproutingRule()
{

}

template<unsigned DIM>
void AbstractSproutingRule<DIM>::SetSproutingProbability(double probability)
{
    mSproutingProbability = probability;
}

template<unsigned DIM>
void AbstractSproutingRule<DIM>::SetNodes(std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes)
{
    mNodes = nodes;
}

template<unsigned DIM>
std::vector<bool> AbstractSproutingRule<DIM>::WillSprout()
{
    std::vector<bool> sprout_flags;
    for(unsigned idx = 0; idx < mNodes.size(); idx++)
    {
        double prob = RandomNumberGenerator::Instance()->ranf();
        bool will_sprout = false;
        // Only non-tip nodes can sprout
        if(mNodes[idx]->GetNumberOfSegments()==2 && prob < mSproutingProbability)
        {
            will_sprout = true;
        }
        sprout_flags.push_back(will_sprout);
    }
    return sprout_flags;
}

template<unsigned DIM>
std::vector<c_vector<double, DIM> > AbstractSproutingRule<DIM>::GetSproutDirection(std::vector<bool> sprout_indices)
{

    std::vector<bool> indices;
    if(sprout_indices.size()>0)
    {
        indices = sprout_indices;
    }
    else
    {
        indices = WillSprout();
    }

    std::vector<c_vector<double, DIM> > directions;
    for(unsigned idx = 0; idx < mNodes.size(); idx++)
    {
        if(indices[idx])
        {
            c_vector<double, DIM> sprout_direction;
            if(RandomNumberGenerator::Instance()->ranf()>=0.5)
            {
                sprout_direction = unit_vector<double>(DIM,0);
            }
            else
            {
                sprout_direction = -unit_vector<double>(DIM,0);
            }
            directions.push_back(sprout_direction);
        }
        else
        {
            directions.push_back(zero_vector<double>(DIM));
        }
    }
    return directions;
}

// Explicit instantiation
template class AbstractSproutingRule<2> ;
template class AbstractSproutingRule<3> ;
