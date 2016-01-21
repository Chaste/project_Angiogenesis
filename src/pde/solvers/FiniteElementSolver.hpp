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

#ifndef FINITEELEMENTSOLVER_HPP_
#define FINITEELEMENTSOLVER_HPP_

#include "SmartPointers.hpp"
#include "AbstractHybridSolver.hpp"
#include "HybridMesh.hpp"
#include "Part.hpp"
#include "HybridLinearEllipticPde.hpp"
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkUnstructuredGrid.h>

template<unsigned DIM>
class FiniteElementSolver : public AbstractHybridSolver<DIM>
{
    using AbstractHybridSolver<DIM>::Solve;

    std::vector<double> mFeSolution;
    vtkSmartPointer<vtkUnstructuredGrid> mFeVtkSolution;
    boost::shared_ptr<HybridMesh<DIM, DIM> > mpMesh;

    std::vector<double> mDomainBounds;
    double mDomainOutOfBoundsValue;

    boost::shared_ptr<CaVascularNetwork<DIM> > mpNetwork;
    double mNetworkBounds;
    bool mUseNewton;
    bool mUseLinearSolveForGuess;
    std::vector<double> mGuess;

public:

    FiniteElementSolver();

    ~FiniteElementSolver();

    /*
     * Construct a new instance of the class and return a shared pointer to it.
     */
    static boost::shared_ptr<FiniteElementSolver<DIM> > Create();

    std::vector<double> GetSolutionAtPoints(std::vector<c_vector<double, DIM> > samplePoints, bool useVtkSampling = true);

    std::vector<double> GetSolutionAtPointsUsingVtk(std::vector<c_vector<double, DIM> > samplePoints);

    std::vector<double> GetSolutionOnRegularGrid(boost::shared_ptr<RegularGrid<DIM, DIM> > pGrid, bool useVtkSampling = true);

    std::vector<double> GetNodalSolution()
    {
        return mFeSolution;
    }

    vtkSmartPointer<vtkImageData> GetVtkRegularGridSolution(boost::shared_ptr<RegularGrid<DIM, DIM> > pGrid, bool useVtkSampling = true);

    void SetSamplingDomainBounds(std::vector<double> domainBounds, double domainOutOfBoundsValue);

    void SetSamplingNetworkBounds(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork, double networkBounds);

    void SetMesh(boost::shared_ptr<HybridMesh<DIM, DIM> > pMesh);

    void Solve();

    void SetGuess(std::vector<double> guess)
    {
        mGuess = guess;
    }

    void SetUseSimpleNetonSolver(bool useNewton = true)
    {
        mUseNewton = useNewton;
    }
    void SetUseLinearSolveForGuess(bool useLinearSolve = true)
    {
        mUseLinearSolveForGuess = useLinearSolve;
    }

    void Update();

    void UpdateCellData();

    void Setup();
    
	void ReadSolution();

private:

    void Write();
};

#endif /* FINITEELEMENTSOLVER_HPP_ */
