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
#include "SimplePetscNonlinearSolver.hpp"

template<unsigned DIM>
PetscErrorCode HyrbidFiniteDifference_ComputeResidual(SNES snes, Vec solution_guess, Vec residual, void* pContext);

template<unsigned DIM>
PetscErrorCode HyrbidFiniteDifference_ComputeJacobian(SNES snes, Vec input, Mat* pJacobian, Mat* pPreconditioner, MatStructure* pMatStructure, void* pContext);

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
    // Set up the grid and PDE
    if(!this->mpRegularGrid)
    {
        EXCEPTION("This solver needs a regular grid to be set before calling Setup.");
    }

    if(!this->mpPde and !this->mpNonLinearPde)
    {
        EXCEPTION("This solver needs a PDE to be set before calling Setup.");
    }

    if(this->mCellPopulationIsSet)
    {
        this->mpRegularGrid->SetCellPopulation(*(this->mpCellPopulation));
    }

    if(this->mpNetwork)
    {
    	this->mpRegularGrid->SetVesselNetwork(this->mpNetwork);
    }

    if(this->mpPde)
    {
        this->mpPde->SetRegularGrid(this->mpRegularGrid);
    }
    else
    {
        this->mpNonLinearPde->SetRegularGrid(this->mpRegularGrid);
    }

    // Set up the boundary conditions
    mpBoundaryConditions = boost::shared_ptr<std::vector<std::pair<bool, double> > > (new std::vector<std::pair<bool, double> >(this->mpRegularGrid->GetNumberOfPoints()));
    for(unsigned idx=0; idx<this->mpRegularGrid->GetNumberOfPoints(); idx++)
    {
        (*mpBoundaryConditions)[idx] = std::pair<bool, double>(false, 0.0);
    }
    for(unsigned bound_index=0; bound_index<this->mBoundaryConditions.size(); bound_index++)
    {
        this->mBoundaryConditions[bound_index]->SetRegularGrid(this->mpRegularGrid);
    }

    // Set up the vtk solution grid
    AbstractRegularGridHybridSolver<DIM>::Setup();

    // Update the source strengths and boundary conditions;
    Update();

    this->IsSetupForSolve = true;
}

template<unsigned DIM>
void FiniteDifferenceSolver<DIM>::Update()
{
    // Update the PDE source strengths
    if(this->mpPde)
    {
        this->mpPde->UpdateDiscreteSourceStrengths();
    }
    else
    {
        this->mpNonLinearPde->UpdateDiscreteSourceStrengths();
    }

    // Update the boundary conditions
    for(unsigned bound_index=0; bound_index<this->mBoundaryConditions.size(); bound_index++)
    {
        this->mBoundaryConditions[bound_index]->UpdateRegularGridBoundaryConditions(mpBoundaryConditions);
    }
}

template<unsigned DIM>
void FiniteDifferenceSolver<DIM>::DoLinearSolve()
{
    // Set up the system
    unsigned number_of_points = this->mpRegularGrid->GetNumberOfPoints();
    unsigned extents_x = this->mpRegularGrid->GetExtents()[0];
    unsigned extents_y = this->mpRegularGrid->GetExtents()[1];
    unsigned extents_z = this->mpRegularGrid->GetExtents()[2];
    double spacing = this->mpRegularGrid->GetSpacing();
    double diffusion_term = 0.0;
    if(this->mpPde)
    {
        diffusion_term = this->mpPde->ComputeIsotropicDiffusionTerm() / (spacing * spacing);
    }
    else
    {
        diffusion_term = this->mpNonLinearPde->ComputeIsotropicDiffusionTerm() / (spacing * spacing);
    }

    LinearSystem linear_system(number_of_points, 7);
    for (unsigned i = 0; i < extents_z; i++) // Z
    {
        for (unsigned j = 0; j < extents_y; j++) // Y
        {
            for (unsigned k = 0; k < extents_x; k++) // X
            {
                unsigned grid_index = this->mpRegularGrid->Get1dGridIndex(k, j, i);

                c_vector<double, DIM> location = this->mpRegularGrid->GetLocation(k ,j, i);
                linear_system.AddToMatrixElement(grid_index, grid_index, this->mpPde->ComputeLinearInUCoeffInSourceTerm(grid_index) - 6.0 * diffusion_term);

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
                linear_system.SetRhsVectorElement(grid_index, -this->mpPde->ComputeConstantInUSourceTerm(grid_index));
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
    this->mPointSolution = std::vector<double>(number_of_points, 0.0);
    for (unsigned row = 0; row < number_of_points; row++)
    {
        this->mPointSolution[row] = soln_repl[row];
    }

    std::map<std::string, std::vector<double> > data;
    data[this->mLabel] = this->mPointSolution;
    this->UpdateSolution(data);

    if (this->mWriteSolution)
    {
        this->Write();
    }
}

template<unsigned DIM>
void FiniteDifferenceSolver<DIM>::Solve()
{
    if(!this->IsSetupForSolve)
    {
        Setup();
    }

    if(this->mpPde)
    {
        DoLinearSolve();
    }
    else
    {
        // Set up initial Guess
        unsigned number_of_points = this->mpRegularGrid->GetNumberOfPoints();
        //Vec initial_guess=PetscTools::CreateAndSetVec(number_of_points, 1.0);

        //SimplePetscNonlinearSolver solver_petsc;
        //int length = 7;
        //Vec answer_petsc = solver_petsc.Solve(&HyrbidFiniteDifference_ComputeResidual<DIM>,
         //                                     &HyrbidFiniteDifference_ComputeJacobian<DIM>, initial_guess, length, this);

        //ReplicatableVector soln_repl(answer_petsc);

        // Populate the solution vector
        this->mPointSolution = std::vector<double>(number_of_points, 0.0);
        for (unsigned row = 0; row < number_of_points; row++)
        {
           // this->mPointSolution[row] = soln_repl[row];
        }

        std::map<std::string, std::vector<double> > data;
        data[this->mLabel] = this->mPointSolution;
        this->UpdateSolution(data);

        if (this->mWriteSolution)
        {
            this->Write();
        }
    }
}

template<unsigned DIM>
PetscErrorCode HyrbidFiniteDifference_ComputeResidual(SNES snes, Vec solution_guess, Vec residual, void* pContext)
{
    FiniteDifferenceSolver<DIM>* solver = (FiniteDifferenceSolver<DIM>*) pContext;

    unsigned extents_x = solver->GetGrid()->GetExtents()[0];
    unsigned extents_y = solver->GetGrid()->GetExtents()[1];
    unsigned extents_z = solver->GetGrid()->GetExtents()[2];
    double spacing = solver->GetGrid()->GetSpacing();
    double diffusion_term = solver->GetNonLinearPde()->ComputeIsotropicDiffusionTerm() / (spacing * spacing);

    // Get the residual vector
    PetscVecTools::Zero(residual);
    for (unsigned i = 0; i < extents_z; i++) // Z
    {
        for (unsigned j = 0; j < extents_y; j++) // Y
        {
            for (unsigned k = 0; k < extents_x; k++) // X
            {
                unsigned grid_index = solver->GetGrid()->Get1dGridIndex(k, j, i);
                double grid_guess = PetscVecTools::GetElement(solution_guess, grid_index);

                c_vector<double, DIM> location = solver->GetGrid()->GetLocation(k ,j, i);
                PetscVecTools::AddToElement(residual, grid_index, grid_guess * (- 6.0 * diffusion_term) +
                                                solver->GetNonLinearPde()->ComputeNonlinearSourceTerm(grid_index, grid_guess));

                // Assume no flux on domain boundaries by default
                // No flux at x bottom
                if (k > 0)
                {
                    double neighbour_guess = PetscVecTools::GetElement(solution_guess, grid_index - 1);
                    PetscVecTools::AddToElement(residual, grid_index, neighbour_guess * diffusion_term);
                }
                else
                {
                    PetscVecTools::AddToElement(residual, grid_index, diffusion_term * grid_guess);
                }

                // No flux at x top
                if (k < extents_x - 1)
                {
                    double neighbour_guess = PetscVecTools::GetElement(solution_guess, grid_index + 1);
                    PetscVecTools::AddToElement(residual, grid_index, diffusion_term * neighbour_guess);
                }
                else
                {
                    PetscVecTools::AddToElement(residual, grid_index, diffusion_term * grid_guess);
                }

                // No flux at y bottom
                if (j > 0)
                {
                    double neighbour_guess = PetscVecTools::GetElement(solution_guess, grid_index - extents_x);
                    PetscVecTools::AddToElement(residual, grid_index, diffusion_term* neighbour_guess);
                }
                else
                {
                    PetscVecTools::AddToElement(residual, grid_index, diffusion_term * grid_guess);
                }

                // No flux at y top
                if (j < extents_y - 1)
                {
                    double neighbour_guess = PetscVecTools::GetElement(solution_guess, grid_index + extents_x);
                    PetscVecTools::AddToElement(residual, grid_index, diffusion_term* neighbour_guess);
                }
                else
                {
                    PetscVecTools::AddToElement(residual, grid_index, diffusion_term * grid_guess);
                }

                // No flux at z bottom
                if (i > 0)
                {
                    double neighbour_guess = PetscVecTools::GetElement(solution_guess, grid_index - extents_x * extents_y);
                    PetscVecTools::AddToElement(residual, grid_index, diffusion_term*neighbour_guess);
                }
                else
                {
                    PetscVecTools::AddToElement(residual, grid_index, diffusion_term * grid_guess);
                }

                // No flux at z top
                if (i < extents_z - 1)
                {
                    double neighbour_guess = PetscVecTools::GetElement(solution_guess, grid_index + extents_x * extents_y);
                    PetscVecTools::AddToElement(residual, grid_index, diffusion_term * neighbour_guess);
                }
                else
                {
                    PetscVecTools::AddToElement(residual, grid_index, diffusion_term * grid_guess);
                }
            }
        }
    }


    PetscVecTools::Finalise(residual);

    // Dirichlet Boundary conditions
    std::vector<unsigned> bc_indices;
    for(unsigned idx=0; idx<solver->GetGrid()->GetNumberOfPoints(); idx++)
    {
        if((*(solver->GetRGBoundaryConditions()))[idx].first)
        {
            PetscVecTools::SetElement(residual, idx, 0.0);
        }
    }

    PetscVecTools::Finalise(residual);

    return 0;
}

template<unsigned DIM>
PetscErrorCode HyrbidFiniteDifference_ComputeJacobian(SNES snes, Vec input, Mat* pJacobian, Mat* pPreconditioner, MatStructure* pMatStructure, void* pContext)
{
    Mat jacobian = *pJacobian;

    FiniteDifferenceSolver<DIM>* solver = (FiniteDifferenceSolver<DIM>*) pContext;
    PetscMatTools::Zero(jacobian);
    PetscMatTools::SwitchWriteMode(jacobian);

    unsigned extents_x = solver->GetGrid()->GetExtents()[0];
    unsigned extents_y = solver->GetGrid()->GetExtents()[1];
    unsigned extents_z = solver->GetGrid()->GetExtents()[2];
    double spacing = solver->GetGrid()->GetSpacing();
    double diffusion_term = solver->GetNonLinearPde()->ComputeIsotropicDiffusionTerm() / (spacing * spacing);

    // Get the residual vector
    for (unsigned i = 0; i < extents_z; i++) // Z
    {
        for (unsigned j = 0; j < extents_y; j++) // Y
        {
            for (unsigned k = 0; k < extents_x; k++) // X
            {
                unsigned grid_index = solver->GetGrid()->Get1dGridIndex(k, j, i);
                double grid_guess = PetscVecTools::GetElement(input, grid_index);

                c_vector<double, DIM> location = solver->GetGrid()->GetLocation(k ,j, i);
                PetscMatTools::AddToElement(jacobian, grid_index, grid_index, - 6.0 * diffusion_term +
                                            solver->GetNonLinearPde()->ComputeNonlinearSourceTermPrime(grid_index, grid_guess));

                // Assume no flux on domain boundaries by default
                // No flux at x bottom
                if (k > 0)
                {
                    PetscMatTools::AddToElement(jacobian, grid_index, grid_index - 1, diffusion_term);
                }
                else
                {
                    PetscMatTools::AddToElement(jacobian, grid_index, grid_index, diffusion_term);
                }

                // No flux at x top
                if (k < extents_x - 1)
                {
                    PetscMatTools::AddToElement(jacobian, grid_index, grid_index + 1, diffusion_term);
                }
                else
                {
                    PetscMatTools::AddToElement(jacobian, grid_index, grid_index, diffusion_term);
                }

                // No flux at y bottom
                if (j > 0)
                {
                    PetscMatTools::AddToElement(jacobian, grid_index, grid_index - extents_x, diffusion_term);
                }
                else
                {
                    PetscMatTools::AddToElement(jacobian, grid_index, grid_index, diffusion_term);
                }

                // No flux at y top
                if (j < extents_y - 1)
                {
                    PetscMatTools::AddToElement(jacobian, grid_index, grid_index + extents_x, diffusion_term);
                }
                else
                {
                    PetscMatTools::AddToElement(jacobian, grid_index, grid_index, diffusion_term);
                }

                // No flux at z bottom
                if (i > 0)
                {
                    PetscMatTools::AddToElement(jacobian, grid_index, grid_index - extents_x * extents_y, diffusion_term);
                }
                else
                {
                    PetscMatTools::AddToElement(jacobian, grid_index,grid_index, diffusion_term);
                }

                // No flux at z top
                if (i < extents_z - 1)
                {
                    PetscMatTools::AddToElement(jacobian, grid_index, grid_index + extents_x * extents_y, diffusion_term);
                }
                else
                {
                    PetscMatTools::AddToElement(jacobian, grid_index, grid_index, diffusion_term);
                }
            }
        }
    }

    // Apply the boundary conditions
    std::vector<unsigned> bc_indices;
    for(unsigned idx=0; idx<solver->GetGrid()->GetNumberOfPoints(); idx++)
    {
        if((*(solver->GetRGBoundaryConditions()))[idx].first)
        {
            bc_indices.push_back(idx);
        }
    }
    PetscMatTools::ZeroRowsWithValueOnDiagonal(jacobian, bc_indices, 1.0);

    PetscMatTools::Finalise(jacobian);
    return 0;
}


// Explicit instantiation
template class FiniteDifferenceSolver<2>;
template class FiniteDifferenceSolver<3>;
