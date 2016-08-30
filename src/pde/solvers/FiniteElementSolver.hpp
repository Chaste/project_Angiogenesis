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
#include "AbstractDiscreteContinuumSolver.hpp"
#include "DiscreteContinuumMesh.hpp"
#include "Part.hpp"
#include "DiscreteContinuumLinearEllipticPde.hpp"
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkUnstructuredGrid.h>

template<unsigned DIM>
class FiniteElementSolver : public AbstractDiscreteContinuumSolver<DIM>
{
    using AbstractDiscreteContinuumSolver<DIM>::Solve;

    // Keep the nodal solution
    std::vector<double> mFeSolution;

    // The solution in the form of a vtk grid
    vtkSmartPointer<vtkUnstructuredGrid> mFeVtkSolution;

    // The finite element mesh
    boost::shared_ptr<DiscreteContinuumMesh<DIM, DIM> > mpMesh;

    // Use the Chaste newton solver
    bool mUseNewton;

    // Use the linear solve as a guess
    bool mUseLinearSolveForGuess;

    // An initial solution guess
    std::vector<double> mGuess;

public:

    /*
     * Constructor
     */
    FiniteElementSolver();

    /*
     * Destructor
     */
    ~FiniteElementSolver();

    /*
     * Construct a new instance of the class and return a shared pointer to it.
     */
    static boost::shared_ptr<FiniteElementSolver<DIM> > Create();

    /*
     * Overridden solution at points
     */
    std::vector<double> GetSolutionAtPoints(std::vector<DimensionalChastePoint<DIM> > samplePoints);

    std::vector<double> GetSolutionAtGridPoints(boost::shared_ptr<RegularGrid<DIM, DIM> > pGrid);

    std::vector<double> GetNodalSolution();

    void ReadSolution();

    void Setup();

    void SetMesh(boost::shared_ptr<DiscreteContinuumMesh<DIM, DIM> > pMesh);

    void Solve();

    void SetGuess(std::vector<double> guess);

    void SetUseSimpleNetonSolver(bool useNewton = true);

    void SetUseLinearSolveForGuess(bool useLinearSolve = true);

    void Update();

    void UpdateCellData();

private:

    void Write();
};

#endif /* FINITEELEMENTSOLVER_HPP_ */
