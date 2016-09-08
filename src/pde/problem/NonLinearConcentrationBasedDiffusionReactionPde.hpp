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

#ifndef NonLinearConcentrationBasedDiffusionReactionPde_HPP_
#define NonLinearConcentrationBasedDiffusionReactionPde_HPP_

#include <string>
#include "ChastePoint.hpp"
#include "UblasIncludes.hpp"
#include "SmartPointers.hpp"
#include "UblasVectorInclude.hpp"
#include "DiscreteSource.hpp"
#include "GeometryTools.hpp"
#include "RegularGrid.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractNonlinearEllipticPde.hpp"

/**
 * Non-Linear Elliptic PDE with both continuum and discrete source terms.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class NonLinearConcentrationBasedDiffusionReactionPde : public AbstractNonlinearEllipticPde<SPACE_DIM>
{
    /**
     * The diffusion tensor
     */
    c_matrix<double, SPACE_DIM, SPACE_DIM> mDiffusionTensor;

    /**
     * The diffusion constant for isotropic diffusion
     */
    double mDiffusivity;

    double mThreshold;

    /**
     * The continuum constant in U term, discrete terms are added to this.
     */
    double mConstantInUTerm;

    /**
     * The continuum linear in U term, discrete terms are added to this.
     */
    double mLinearInUTerm;

    /**
     * The name of the quantity in the pde
     */
    std::string mVariableName;

    /**
     * The collection of discrete sources for addition to the continuum terms
     */
    std::vector<boost::shared_ptr<DiscreteSource<SPACE_DIM> > > mDiscreteSources;

    /**
     * The grid for solvers using regular grids
     */
    boost::shared_ptr<RegularGrid<ELEMENT_DIM, SPACE_DIM> > mpRegularGrid;

    /**
     * The mesh for solvers using finite element meshes
     */
    boost::shared_ptr<TetrahedralMesh<ELEMENT_DIM, SPACE_DIM> > mpMesh;

    /**
     * Whether to use a regular grid or mesh for discrete source calculations
     */
    bool mUseRegularGrid;

    /**
     * The constant source strengths for each point on the grid or mesh
     */
    std::vector<double> mDiscreteConstantSourceStrengths;

    /**
     * The linear source strengths for each point on the grid or mesh
     */
    std::vector<double> mDiscreteLinearSourceStrengths;

public:

    /**
     * Constructor
     */
    NonLinearConcentrationBasedDiffusionReactionPde();

    /**
     * Factory Constructor
     * @return a pointer to an instance of the pde
     */
    static boost::shared_ptr<NonLinearConcentrationBasedDiffusionReactionPde<ELEMENT_DIM, SPACE_DIM> > Create();

    /**
     * Add a discrete source to the pde
     * @param pDiscreteSource a pointer the discrete source
     */
    void AddDiscreteSource(boost::shared_ptr<DiscreteSource<SPACE_DIM> > pDiscreteSource);


    void SetThreshold(double threshold)
    {
        mThreshold = threshold;
    }

    double GetThreshold()
    {
        return mThreshold;
    }

    double ComputeLinearSourceTerm(const ChastePoint<SPACE_DIM>& rX)
    {
        return 0.0;
    }

    double ComputeNonlinearSourceTerm(const ChastePoint<SPACE_DIM>& rX, double u)
    {
//        if (u>mThreshold)
//        {
//            return mConstantInUTerm;
//        }
//        else
//        {
//            return mConstantInUTerm * (u/mThreshold);
//        }
        if(u<0.0)
        {
            return 0.0;
        }
        else
        {
            return mConstantInUTerm * u / (mThreshold + u);
        }
    }

    double ComputeNonlinearSourceTermPrime(const ChastePoint<SPACE_DIM>& rX, double u)
    {
//        if (u>mThreshold)
//        {
//            return 0.0;
//        }
//        else
//        {
//            return mConstantInUTerm/mThreshold;
//        }
        if(u<0.0)
        {
            return mConstantInUTerm * mThreshold/((mThreshold+ 0.0)*(mThreshold + 0.0));
        }
        else
        {
            return mConstantInUTerm * mThreshold/((mThreshold + u)*(mThreshold + u));
        }

    }

    double ComputeNonlinearSourceTerm(unsigned gridIndex, double u)
    {
//        if (u>mThreshold)
//        {
//            return mConstantInUTerm;
//        }
//        else
//        {
//            return mConstantInUTerm * (u/mThreshold);
//        }

        if(u<0.0)
        {
            return mDiscreteConstantSourceStrengths[gridIndex];
        }
        else
        {
            return mConstantInUTerm * u / (mThreshold + u) + mDiscreteConstantSourceStrengths[gridIndex] + mDiscreteLinearSourceStrengths[gridIndex]*u;
        }
    }

    double ComputeNonlinearSourceTermPrime(unsigned gridIndex, double u)
    {
//        if (u>mThreshold)
//        {
//            return 0.0;
//        }
//        else
//        {
//            return mConstantInUTerm/mThreshold;
//        }
        if(u<0.0)
        {
            return mConstantInUTerm * mThreshold/((mThreshold + 0.0)*(mThreshold + 0.0)) + mDiscreteLinearSourceStrengths[gridIndex];
        }
        else
        {
            return mConstantInUTerm * mThreshold/((mThreshold + u)*(mThreshold + u)) + mDiscreteLinearSourceStrengths[gridIndex];
        }

    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTermPrime(const ChastePoint<SPACE_DIM>& rX, double u)
    {
        return zero_matrix<double>(SPACE_DIM);
    }
    /**
     * Overwritten method to return the constant in U contribution to the Chaste FE solver
     * @param rX grid location
     * @param pElement pointer to containing element
     * @return source strength
     */
    double ComputeConstantInUSourceTerm(const ChastePoint<SPACE_DIM>& rX, Element<ELEMENT_DIM, SPACE_DIM>* pElement);

    /**
     * Overwritten method to return the constant in U contribution to the regular grid solvers
     * @param gridIndex grid index
     * @return source strength
     */
    double ComputeConstantInUSourceTerm(unsigned gridIndex=0);

    /**
     * Overwritten method to return the diffusion term to the Chaste FE solver
     * @return the diffusion matrix
     */
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>&, double u);

    /**
     * Return the diffusion constant for isotropic diffusion
     * @return the diffusion constant
     */
    double ComputeIsotropicDiffusionTerm();

    /**
     * Overwritten method to return the linear in U contribution to the Chaste FE solver
     * @param rX grid location
     * @param pElement pointer to containing element
     * @return source strength
     */
    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<SPACE_DIM>& rX, Element<ELEMENT_DIM, SPACE_DIM>* pElement);

    /**
     * Overwritten method to return the linear in U contribution to the regular grid solvers
     * @param gridIndex grid index
     * @return source strength
     */
    double ComputeLinearInUCoeffInSourceTerm(unsigned gridIndex=0);

    /**
     * Return the collection of discrete sources
     * @return vector of pointers to the discrete sources
     */
    std::vector<boost::shared_ptr<DiscreteSource<SPACE_DIM> > > GetDiscreteSources();

    /**
     * Return the name of the pde variable
     * @return the name of the pde variable
     */
    const std::string& GetVariableName();

    /**
     * Set the continuum constant in U term
     * @param constantInUTerm the continuum constant in U term
     */
    void SetContinuumConstantInUTerm(double constantInUTerm);

    /**
     * Set the isotropic diffusion constant
     * @param diffusivity the isotropic diffusion constant
     */
    void SetIsotropicDiffusionConstant(double diffusivity);

    /**
     * Set the linear constant in U term
     * @param linearInUTerm the linear constant in U term
     */
    void SetContinuumLinearInUTerm(double linearInUTerm);

    /**
     * Set the regular grid
     * @param pRegularGrid the regular grid
     */
    void SetRegularGrid(boost::shared_ptr<RegularGrid<ELEMENT_DIM, SPACE_DIM> > pRegularGrid);

    /**
     * Set the finite element mesh
     * @param pMesh the finite element mesh
     */
    void SetMesh(boost::shared_ptr<TetrahedralMesh<ELEMENT_DIM, SPACE_DIM> > pMesh);

    /**
     * Set whether to use a regular grid
     * @param useRegularGrid whether to use a regular grid
     */
    void SetUseRegularGrid(bool useRegularGrid);

    /**
     * Set the name of the pde variable
     * @param rVariableName the name of the pde variable
     */
    void SetVariableName(const std::string& rVariableName);

    /**
     * Update the discrete source strengths
     */
    void UpdateDiscreteSourceStrengths();

};

#endif /*NonLinearConcentrationBasedDiffusionReactionPde_HPP_*/
