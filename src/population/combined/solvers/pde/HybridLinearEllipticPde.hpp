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

#ifndef HYBRIDLINEARELLIPTICPDE_HPP_
#define HYBRIDLINEARELLIPTICPDE_HPP_

#include <string>
#include "ChastePoint.hpp"
#include "UblasIncludes.hpp"
#include "SmartPointers.hpp"
#include "UblasVectorInclude.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "SimpleCellPopulation.hpp"
#include "CaVascularNetwork.hpp"
#include "DiscreteSource.hpp"
#include "GeometryTools.hpp"

/*
 * Linear Elliptic PDE with discrete or averaged cell and
 * vessel source terms.
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class HybridLinearEllipticPde : public AbstractLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>
{
    c_matrix<double, SPACE_DIM, SPACE_DIM> mDiffusionTensor;
    double mDiffusivity;
    double mConstantInUTerm;
    double mLinearInUTerm;
    std::string mVariableName;
    boost::shared_ptr<SimpleCellPopulation<SPACE_DIM> > mpPopulation;
    boost::shared_ptr<CaVascularNetwork<SPACE_DIM> > mpNetwork;

    /* The discrete sources
    */
    std::vector<boost::shared_ptr<DiscreteSource<SPACE_DIM> > > mDiscreteSources;

public:

    HybridLinearEllipticPde() :
            mDiffusionTensor(identity_matrix<double>(SPACE_DIM)),
            mDiffusivity(1.e-3),
            mConstantInUTerm(0.0),
            mLinearInUTerm(0.0),
            mVariableName("Default"),
            mpPopulation(),
            mpNetwork(),
            mDiscreteSources()
    {
        mDiffusionTensor *= mDiffusivity;
    }

    static boost::shared_ptr<HybridLinearEllipticPde<ELEMENT_DIM, SPACE_DIM> > Create()
    {
        MAKE_PTR(HybridLinearEllipticPde<ELEMENT_DIM>, pSelf);
        return pSelf;
    }

    void AddDiscreteSource(boost::shared_ptr<DiscreteSource<SPACE_DIM> > pDiscreteSource)
    {
        mDiscreteSources.push_back(pDiscreteSource);
    }

    double ComputeConstantInUSourceTerm(const ChastePoint<SPACE_DIM>& rX, Element<ELEMENT_DIM, SPACE_DIM>* pElement)
    {
        return mConstantInUTerm;
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<SPACE_DIM>& rX, Element<ELEMENT_DIM, SPACE_DIM>* pElement)
    {
        return mLinearInUTerm;
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>&)
    {
        return mDiffusionTensor;
    }

    double GetDiffusionConstant()
    {
        return mDiffusivity;
    }

    double GetConstantInUTerm(c_vector<double, SPACE_DIM> location = zero_vector<double>(SPACE_DIM), double spacing = 0.0)
    {
        double consumption_term = 0.0;
        for(unsigned idx=0; idx<mDiscreteSources.size(); idx++)
        {
            if(mDiscreteSources[idx]->GetType()==SourceType::SOLUTION && !mDiscreteSources[idx]->IsLinearInSolution())
            {
                consumption_term =  -1.e-5 * mDiscreteSources[idx]->GetValue(location).second;
            }

            if(mDiscreteSources[idx]->GetType()==SourceType::MULTI_POINT && !mDiscreteSources[idx]->IsLinearInSolution())
            {
                std::pair<bool, double> result = mDiscreteSources[idx]->GetValue(location, spacing/2.0);
                if(result.first)
                {
                    consumption_term = result.second;
                }
            }
        }
        return mConstantInUTerm - consumption_term;
    }

    std::vector<boost::shared_ptr<DiscreteSource<SPACE_DIM> > > GetDiscreteSources()
    {
        return mDiscreteSources;
    }

    double GetLinearInUTerm(c_vector<double, SPACE_DIM> location  = zero_vector<double>(SPACE_DIM), double spacing = 0.0)
    {
        double consumption_term = 0.0;
        for(unsigned idx=0; idx<mDiscreteSources.size(); idx++)
        {
            if(mDiscreteSources[idx]->GetType()==SourceType::SOLUTION && mDiscreteSources[idx]->IsLinearInSolution())
            {
                consumption_term =  -1.e-5 * mDiscreteSources[idx]->GetValue(location).second;
            }

            if(mDiscreteSources[idx]->GetType()==SourceType::MULTI_POINT && mDiscreteSources[idx]->IsLinearInSolution())
            {
                std::pair<bool, double> result = mDiscreteSources[idx]->GetValue(location, spacing/2.0);
                if(result.first)
                {
                    consumption_term = result.second;
                }
            }
        }

        return mLinearInUTerm - consumption_term;
    }

    const std::string& GetVariableName()
    {
        return mVariableName;
    }

    void SetCellPopulation(boost::shared_ptr<SimpleCellPopulation<SPACE_DIM> > pPopulation)
    {
        mpPopulation = pPopulation;
    }

    void SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<SPACE_DIM> > pNetwork)
    {
        mpNetwork = pNetwork;
    }

    void SetConstantInUTerm(double constantInUTerm)
    {
        mConstantInUTerm = constantInUTerm;
    }

    void SetLinearInUTerm(double linearInUTerm)
    {
        mLinearInUTerm = linearInUTerm;
    }

    void SetDiffusionConstant(double diffusivity)
    {
        mDiffusivity = diffusivity;
        mDiffusionTensor = identity_matrix<double>(SPACE_DIM)* mDiffusivity;
    }

    void SetVariableName(const std::string& rVariableName)
    {
        mVariableName = rVariableName;
    }

};

#endif /*HYBRIDLINEARELLIPTICPDE_HPP_*/
