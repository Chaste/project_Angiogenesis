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
#include "DiscreteContinuumLinearEllipticPde.hpp"
#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::DiscreteContinuumLinearEllipticPde() :
            AbstractLinearEllipticPde<ELEMENT_DIM, ELEMENT_DIM>(),
            mDiffusionTensor(identity_matrix<double>(SPACE_DIM)),
            mDiffusivity(0.003),
            mConstantInUTerm(0.0),
            mLinearInUTerm(0.0),
            mVariableName("Default"),
            mDiscreteSources(),
            mpRegularGrid(),
            mpMesh(),
            mUseRegularGrid(true),
            mDiscreteConstantSourceStrengths(std::vector<double>(1,0.0)),
            mDiscreteLinearSourceStrengths(std::vector<double>(1,0.0))
{
    mDiffusionTensor *= mDiffusivity;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
boost::shared_ptr<DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM> > DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::Create()
{
    MAKE_PTR(DiscreteContinuumLinearEllipticPde<ELEMENT_DIM>, pSelf);
    return pSelf;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::AddDiscreteSource(boost::shared_ptr<DiscreteSource<SPACE_DIM> > pDiscreteSource)
{
    mDiscreteSources.push_back(pDiscreteSource);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::ComputeConstantInUSourceTerm(const ChastePoint<SPACE_DIM>& rX, Element<ELEMENT_DIM, SPACE_DIM>* pElement)
{
    return mConstantInUTerm;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::ComputeConstantInUSourceTerm(unsigned gridIndex)
{
    if(gridIndex >= mDiscreteLinearSourceStrengths.size())
    {
        EXCEPTION("Requested out of bound grid index in discrete sources. Maybe you forgot to update the source strengths.");
    }
    return mConstantInUTerm + mDiscreteConstantSourceStrengths[gridIndex];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double, SPACE_DIM, SPACE_DIM> DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>&)
{
    return mDiffusionTensor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::ComputeLinearInUCoeffInSourceTerm(const ChastePoint<SPACE_DIM>& rX, Element<ELEMENT_DIM, SPACE_DIM>* pElement)
{
    return mLinearInUTerm;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::ComputeLinearInUCoeffInSourceTerm(unsigned gridIndex)
{
    if(gridIndex >= mDiscreteLinearSourceStrengths.size())
    {
        EXCEPTION("Requested out of bound grid index in discrete sources. Maybe you forgot to update the source strengths.");
    }
    return mLinearInUTerm + mDiscreteLinearSourceStrengths[gridIndex];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::ComputeIsotropicDiffusionTerm()
{
    return mDiffusivity;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<boost::shared_ptr<DiscreteSource<SPACE_DIM> > > DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::GetDiscreteSources()
{
    return mDiscreteSources;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::string& DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::GetVariableName()
{
    return mVariableName;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::SetContinuumConstantInUTerm(double constantInUTerm)
{
    mConstantInUTerm = constantInUTerm;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::SetIsotropicDiffusionConstant(double diffusivity)
{
    mDiffusivity = diffusivity;
    mDiffusionTensor = identity_matrix<double>(SPACE_DIM)* mDiffusivity;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::SetContinuumLinearInUTerm(double linearInUTerm)
{
    mLinearInUTerm = linearInUTerm;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::SetRegularGrid(boost::shared_ptr<RegularGrid<ELEMENT_DIM, SPACE_DIM> > pRegularGrid)
{
    mpRegularGrid = pRegularGrid;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::SetMesh(boost::shared_ptr<TetrahedralMesh<ELEMENT_DIM, SPACE_DIM> > pMesh)
{
    mpMesh = pMesh;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::SetUseRegularGrid(bool useRegularGrid)
{
    mUseRegularGrid = useRegularGrid;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::SetVariableName(const std::string& rVariableName)
{
    mVariableName = rVariableName;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DiscreteContinuumLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::UpdateDiscreteSourceStrengths()
{
    if(mUseRegularGrid)
    {
        if(!mpRegularGrid)
        {
            EXCEPTION("A grid has not been set for the determination of source strengths.");
        }
        mDiscreteConstantSourceStrengths = std::vector<double>(mpRegularGrid->GetNumberOfPoints(), 0.0);
        mDiscreteLinearSourceStrengths = std::vector<double>(mpRegularGrid->GetNumberOfPoints(), 0.0);

        for(unsigned idx=0; idx<mDiscreteSources.size(); idx++)
        {
            mDiscreteSources[idx]->SetRegularGrid(mpRegularGrid);
            std::vector<double> result = mDiscreteSources[idx]->GetRegularGridValues();
            if(mDiscreteSources[idx]->IsLinearInSolution())
            {
                std::transform(mDiscreteLinearSourceStrengths.begin( ), mDiscreteLinearSourceStrengths.end( ),
                               result.begin( ), mDiscreteLinearSourceStrengths.begin( ),std::plus<double>( ));
            }
            else
            {
                std::transform(mDiscreteConstantSourceStrengths.begin( ), mDiscreteConstantSourceStrengths.end( ),
                               result.begin( ), mDiscreteConstantSourceStrengths.begin( ),std::plus<double>( ));
            }
        }
    }
    else
    {
        if(!mpMesh)
        {
            EXCEPTION("A mesh has not been set for the determination of source strengths.");
        }

        //mesh size
//        mDiscreteConstantSourceStrengths = std::vector<double>(mpRegularGrid->GetNumberOfPoints(), 0.0);
//        mDiscreteLinearSourceStrengths = std::vector<double>(mpRegularGrid->GetNumberOfPoints(), 0.0);

    }
}

// Explicit instantiation
template class DiscreteContinuumLinearEllipticPde<2>;
template class DiscreteContinuumLinearEllipticPde<3>;