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
 * Redistributions in binary form must reproduce the abovea copyright notice,
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

#include <algorithm>
#include "LinearSteadyStateDiffusionReactionPde.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
LinearSteadyStateDiffusionReactionPde<ELEMENT_DIM, SPACE_DIM>::LinearSteadyStateDiffusionReactionPde() :
            AbstractDiscreteContinuumLinearEllipticPde<ELEMENT_DIM, ELEMENT_DIM>(),
            mLinearInUTerm(0.0 * unit::per_second),
            mDiscreteLinearSourceStrengths()
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
boost::shared_ptr<LinearSteadyStateDiffusionReactionPde<ELEMENT_DIM, SPACE_DIM> > LinearSteadyStateDiffusionReactionPde<ELEMENT_DIM, SPACE_DIM>::Create()
{
    MAKE_PTR(LinearSteadyStateDiffusionReactionPde<ELEMENT_DIM>, pSelf);
    return pSelf;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double LinearSteadyStateDiffusionReactionPde<ELEMENT_DIM, SPACE_DIM>::ComputeLinearInUCoeffInSourceTerm(const ChastePoint<SPACE_DIM>& rX, Element<ELEMENT_DIM, SPACE_DIM>* pElement)
{
    return mLinearInUTerm/unit::per_second;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
units::quantity<unit::rate> LinearSteadyStateDiffusionReactionPde<ELEMENT_DIM, SPACE_DIM>::ComputeLinearInUCoeffInSourceTerm(unsigned gridIndex)
{
    if(gridIndex >= mDiscreteLinearSourceStrengths.size())
    {
        EXCEPTION("Requested out of bound grid index in discrete sources. Maybe you forgot to update the source strengths.");
    }
    return mLinearInUTerm + mDiscreteLinearSourceStrengths[gridIndex];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSteadyStateDiffusionReactionPde<ELEMENT_DIM, SPACE_DIM>::SetContinuumLinearInUTerm(units::quantity<unit::rate> linearInUTerm)
{
    mLinearInUTerm = linearInUTerm;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSteadyStateDiffusionReactionPde<ELEMENT_DIM, SPACE_DIM>::UpdateDiscreteSourceStrengths()
{

    AbstractDiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::UpdateDiscreteSourceStrengths();
    if(this->mUseRegularGrid)
    {
        if(!this->mpRegularGrid)
        {
            EXCEPTION("A grid has not been set for the determination of source strengths.");
        }
        mDiscreteLinearSourceStrengths = std::vector<units::quantity<unit::rate> >(this->mpRegularGrid->GetNumberOfPoints(), 0.0);

        for(unsigned idx=0; idx<this->mDiscreteSources.size(); idx++)
        {
            this->mDiscreteSources[idx]->SetRegularGrid(mpRegularGrid);
            if(this->mDiscreteSources[idx]->IsLinearInSolution())
            {
                std::vector<units::quantity<unit::rate> > result = this->mDiscreteSources[idx]->GetLinearRegularGridValues();
                std::transform(mDiscreteLinearSourceStrengths.begin( ), mDiscreteLinearSourceStrengths.end( ),
                               result.begin( ), mDiscreteLinearSourceStrengths.begin( ),std::plus<units::quantity<unit::rate> >( ));
            }
        }
    }
    else
    {
        if(!this->mpMesh)
        {
            EXCEPTION("A mesh has not been set for the determination of source strengths.");
        }

    }
}

// Explicit instantiation
template class LinearSteadyStateDiffusionReactionPde<2>;
template class LinearSteadyStateDiffusionReactionPde<3>;