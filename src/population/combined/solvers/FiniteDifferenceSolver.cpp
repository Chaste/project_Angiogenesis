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

#include "LinearSystem.hpp"
#include "ReplicatableVector.hpp"
#include "CaVesselSegment.hpp"
#include "CaVascularNetwork.hpp"
#include "HybridLinearEllipticPde.hpp"

#include "FiniteDifferenceSolver.hpp"
#include "Debug.hpp"

template<unsigned DIM>
FiniteDifferenceSolver<DIM>::FiniteDifferenceSolver()
    :   AbstractRegularGridHybridSolver<DIM>()
{

}

template<unsigned DIM>
boost::shared_ptr<FiniteDifferenceSolver<DIM> > FiniteDifferenceSolver<DIM>::Create()
{
    MAKE_PTR(FiniteDifferenceSolver<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
FiniteDifferenceSolver<DIM>::~FiniteDifferenceSolver()
{

}

template<unsigned DIM>
void FiniteDifferenceSolver<DIM>::Setup()
{
    // Set up the PDE
    this->mpPde->SetRegularGrid(this->mpRegularGrid);

    // Set up the boundary conditions
    mpBoundaryConditions = boost::shared_ptr<std::vector<std::pair<bool, double> > > (new std::vector<std::pair<bool, double> >(this->mpRegularGrid->GetNumberOfPoints()));
    for(unsigned idx=0; idx<this->mpRegularGrid->GetNumberOfPoints(); idx++)
    {
        (*mpBoundaryConditions)[idx] = std::pair<bool, double>(false, 0.0);
    }
    for(unsigned bound_index=0; bound_index<this->mDirichletBoundaryConditions.size(); bound_index++)
    {
        this->mDirichletBoundaryConditions[bound_index]->SetRegularGrid(this->mpRegularGrid);
    }

    // Set up the vtk solution grid
    AbstractRegularGridHybridSolver<DIM>::Setup();

    // Update the source strengths and boundary conditions;
    UpdateSourcesAndBoundaryConditions();
}

template<unsigned DIM>
void FiniteDifferenceSolver<DIM>::UpdateSourcesAndBoundaryConditions()
{
    // Update the PDE source strengths
    this->mpPde->UpdateDiscreteSourceStrengths();

    // Update the boundary conditions
    double grid_tolerance = this->mpRegularGrid->GetSpacing() * (std::sqrt(2.0) / 2.0);
    for(unsigned bound_index=0; bound_index<this->mDirichletBoundaryConditions.size(); bound_index++)
    {
        this->mDirichletBoundaryConditions[bound_index]->UpdateRegularGridBoundaryConditions(mpBoundaryConditions, grid_tolerance);
    }
}

template<unsigned DIM>
void FiniteDifferenceSolver<DIM>::Solve(bool writeSolution)
{
    // Set up the system
    unsigned number_of_points = this->mpRegularGrid->GetNumberOfPoints();
    unsigned extents_x = this->mpRegularGrid->GetExtents()[0];
    unsigned extents_y = this->mpRegularGrid->GetExtents()[1];
    unsigned extents_z = this->mpRegularGrid->GetExtents()[2];
    double spacing = this->mpRegularGrid->GetSpacing();
    double diffusion_term = this->mpPde->GetDiffusionConstant() / (spacing * spacing);

    LinearSystem linear_system(number_of_points, 7);
    for (unsigned i = 0; i < extents_z; i++) // Z
    {
        for (unsigned j = 0; j < extents_y; j++) // Y
        {
            for (unsigned k = 0; k < extents_x; k++) // X
            {
                unsigned grid_index = this->mpRegularGrid->Get1dGridIndex(k, j, i);

                c_vector<double, DIM> location = this->mpRegularGrid->GetLocation(k ,j, i);
                linear_system.AddToMatrixElement(grid_index, grid_index, this->mpPde->GetLinearInUTerm(grid_index) - 6.0 * diffusion_term);

                // Assume no flux on domain boundaries by default
                // No flux at x bottom
                if (k > 0)
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index - 1, diffusion_term);
                }
                else
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index, diffusion_term);
                }

                // No flux at x top
                if (k < extents_x - 1)
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index + 1, diffusion_term);
                }
                else
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index, diffusion_term);
                }

                // No flux at y bottom
                if (j > 0)
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index - extents_x, diffusion_term);
                }
                else
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index, diffusion_term);
                }

                // No flux at y top
                if (j < extents_y - 1)
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index + extents_x, diffusion_term);
                }
                else
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index, diffusion_term);
                }

                // No flux at z bottom
                if (i > 0)
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index - extents_x * extents_y, diffusion_term);
                }
                else
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index, diffusion_term);
                }

                // No flux at z top
                if (i < extents_z - 1)
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index + extents_x * extents_y, diffusion_term);
                }
                else
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index, diffusion_term);
                }
                linear_system.SetRhsVectorElement(grid_index, -this->mpPde->GetConstantInUTerm(grid_index));
            }
        }
    }

    // Apply the boundary conditions
    std::vector<unsigned> bc_indices;
    for(unsigned idx=0; idx<this->mpRegularGrid->GetNumberOfPoints(); idx++)
    {
        if((*mpBoundaryConditions)[idx].first)
        {
            bc_indices.push_back(idx);
            linear_system.SetRhsVectorElement(idx, (*mpBoundaryConditions)[idx].second);
        }
    }
    linear_system.ZeroMatrixRowsWithValueOnDiagonal(bc_indices, 1.0);

    // Solve the linear system
    linear_system.AssembleFinalLinearSystem();
    ReplicatableVector soln_repl(linear_system.Solve());

    // Populate the solution vector
    std::vector<double> solution(number_of_points, 0.0);
    for (unsigned row = 0; row < number_of_points; row++)
    {
        solution[row] = soln_repl[row];
    }

    std::map<std::string, std::vector<double> > data;
    data[this->mpPde->GetVariableName()] = solution;
    this->UpdateSolution(data);

    if (writeSolution)
    {
        this->Write();
    }
}

// Explicit instantiation
template class FiniteDifferenceSolver<2>;
template class FiniteDifferenceSolver<3>;
