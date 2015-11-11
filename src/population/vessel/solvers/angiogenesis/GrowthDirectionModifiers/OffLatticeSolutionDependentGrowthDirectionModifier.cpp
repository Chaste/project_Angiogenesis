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
#include "OffLatticeSolutionDependentGrowthDirectionModifier.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"

template<unsigned DIM>
OffLatticeSolutionDependentGrowthDirectionModifier<DIM>::OffLatticeSolutionDependentGrowthDirectionModifier()
    : AbstractGrowthDirectionModifier<DIM>(),
      mpSolver(),
      mProbeLength(5.0)
{

}

template<unsigned DIM>
OffLatticeSolutionDependentGrowthDirectionModifier<DIM>::~OffLatticeSolutionDependentGrowthDirectionModifier()
{

}

template <unsigned DIM>
boost::shared_ptr<OffLatticeSolutionDependentGrowthDirectionModifier<DIM> > OffLatticeSolutionDependentGrowthDirectionModifier<DIM>::Create()
{
    MAKE_PTR(OffLatticeSolutionDependentGrowthDirectionModifier<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
void OffLatticeSolutionDependentGrowthDirectionModifier<DIM>::SetSolver(boost::shared_ptr<AbstractHybridSolver<DIM> > pSolver)
{
    mpSolver = pSolver;
}

template<unsigned DIM>
c_vector<double, DIM> OffLatticeSolutionDependentGrowthDirectionModifier<DIM>::GetGrowthDirection(c_vector<double, DIM> currentDirection,
                                                                                                  boost::shared_ptr<VascularNode<DIM> > pNode)
{

    // If there is a PDE get the direction of highest solution gradient for the specified species
    if(mpSolver)
    {
        // Make points
        std::vector<c_vector<double, DIM> > locations;
        locations.push_back(pNode->GetLocationVector());
        locations.push_back(locations[0] + mProbeLength * unit_vector<double>(DIM,0));
        locations.push_back(locations[0] - mProbeLength * unit_vector<double>(DIM,0));
        locations.push_back(locations[0] + mProbeLength * unit_vector<double>(DIM,1));
        locations.push_back(locations[0] - mProbeLength * unit_vector<double>(DIM,1));
        if(DIM==3)
        {
            locations.push_back(locations[0] + mProbeLength * unit_vector<double>(DIM,2));
            locations.push_back(locations[0] - mProbeLength * unit_vector<double>(DIM,2));
        }

        // Get the solution
        std::vector<double> solutions = mpSolver->GetSolutionAtPoints(locations);

        // Get the gradients
        std::vector<double> gradients;
        for(unsigned idx=1; idx<solutions.size();idx++)
        {
            gradients.push_back((solutions[idx] - solutions[0]) / mProbeLength);
        }

        // Get the index of the max gradient
        double max_grad = 0.0;
        int index = -1;

        for(unsigned idx = 0; idx<gradients.size(); idx++)
        {
            if(gradients[idx]>max_grad)
            {
                max_grad = gradients[idx];
                index = idx;
            }
        }
        if(index == -1)
        {
            return zero_vector<double>(DIM);
        }
        else if(index == 0)
        {
            return unit_vector<double>(DIM,0);
        }
        else if(index == 1)
        {
            return -unit_vector<double>(DIM,0);
        }
        else if(index == 2)
        {
            return unit_vector<double>(DIM,1);
        }
        else if(index == 3)
        {
            return -unit_vector<double>(DIM,1);
        }
        else if(index == 4)
        {
            return unit_vector<double>(DIM,2);
        }
        else if(index == 5)
        {
            return -unit_vector<double>(DIM,2);
        }
        else
        {
            EXCEPTION("Unexpected index while obtaining solution dependent direction vectors.");
        }
    }
    else
    {
        return zero_vector<double>(DIM);
    }
}

// Explicit instantiation
template class OffLatticeSolutionDependentGrowthDirectionModifier<2> ;
template class OffLatticeSolutionDependentGrowthDirectionModifier<3> ;
