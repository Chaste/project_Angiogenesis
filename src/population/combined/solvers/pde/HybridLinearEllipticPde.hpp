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
#include "DiscreteSource.hpp"
#include "GeometryTools.hpp"
#include "RegularGrid.hpp"
#include "TetrahedralMesh.hpp"

/**
 * Linear Elliptic PDE with both continuum and discrete source terms.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class HybridLinearEllipticPde : public AbstractLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>
{
    /**
     * The diffusion tensor
     */
    c_matrix<double, SPACE_DIM, SPACE_DIM> mDiffusionTensor;

    /**
     * The diffusion constant for isotropic diffusion
     */
    double mDiffusivity;

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
    boost::shared_ptr<RegularGrid> mpRegularGrid;

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
    HybridLinearEllipticPde();

    /**
     * Factory Constructor
     * @return a pointer to an instance of the pde
     */
    static boost::shared_ptr<HybridLinearEllipticPde<ELEMENT_DIM, SPACE_DIM> > Create();

    /**
     * Add a discrete source to the pde
     * @param pDiscreteSource a pointer the discrete source
     */
    void AddDiscreteSource(boost::shared_ptr<DiscreteSource<SPACE_DIM> > pDiscreteSource);

    /**
     * Overwritten method to return the constant in U contribution to the Chaste FE solver
     * @param rX grid location
     * @param pElement pointer to containing element
     * @return source strength
     */
    double ComputeConstantInUSourceTerm(const ChastePoint<SPACE_DIM>& rX, Element<ELEMENT_DIM, SPACE_DIM>* pElement);

    /**
     * Overwritten method to return the diffusion term to the Chaste FE solver
     * @return the diffusion matrix
     */
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>&);

    /**
     * Overwritten method to return the linear in U contribution to the Chaste FE solver
     * @param rX grid location
     * @param pElement pointer to containing element
     * @return source strength
     */
    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<SPACE_DIM>& rX, Element<ELEMENT_DIM, SPACE_DIM>* pElement);

    /**
     * Return the constant in U contribution to the hybrid solvers
     * @param location the sample location
     * @param spacing the distance from the sample location demarking 'coincident' discrete points
     * @return the constant in U contribution
     */
    double GetConstantInUTerm(unsigned gridIndex=0);

    /**
     * Return the diffusion constant for isotropic diffusion
     * @return the diffusion constant
     */
    double GetDiffusionConstant();

    /**
     * Return the collection of discrete sources
     * @return vector of pointers to the discrete sources
     */
    std::vector<boost::shared_ptr<DiscreteSource<SPACE_DIM> > > GetDiscreteSources();

    /**
     * Return the linear in U contribution to the hybrid solvers
     * @param location the sample location
     * @param spacing the distance from the sample location demarking 'coincident' discrete points
     * @return the linear in U contribution
     */
    double GetLinearInUTerm(unsigned gridIndex=0);

    /**
     * Return the name of the pde variable
     * @return the name of the pde variable
     */
    const std::string& GetVariableName();

    /**
     * Set the continuum constant in U term
     * @param constantInUTerm the continuum constant in U term
     */
    void SetConstantInUTerm(double constantInUTerm);

    /**
     * Set the isotropic diffusion constant
     * @param diffusivity the isotropic diffusion constant
     */
    void SetDiffusionConstant(double diffusivity);

    /**
     * Set the linear constant in U term
     * @param linearInUTerm the linear constant in U term
     */
    void SetLinearInUTerm(double linearInUTerm);

    /**
     * Set the regular grid
     * @param pRegularGrid the regular grid
     */
    void SetRegularGrid(boost::shared_ptr<RegularGrid> pRegularGrid);

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

#endif /*HYBRIDLINEARELLIPTICPDE_HPP_*/
