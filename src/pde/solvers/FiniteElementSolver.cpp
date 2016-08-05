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

#include <math.h>
#include "SimpleLinearEllipticSolver.hpp"
#include "SimpleNonlinearEllipticSolver.hpp"
#include "SimpleNewtonNonlinearSolver.hpp"
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkProbeFilter.h>
#include <vtkImageData.h>
#include "VtkMeshWriter.hpp"
#include "FiniteElementSolver.hpp"
#include "VtkMeshReader.hpp"

template<unsigned DIM>
FiniteElementSolver<DIM>::FiniteElementSolver()
    : AbstractHybridSolver<DIM>(),
      mFeSolution(),
      mFeVtkSolution(),
      mpMesh(),
      mUseNewton(false),
      mUseLinearSolveForGuess(false),
      mGuess()
{

}

template<unsigned DIM>
FiniteElementSolver<DIM>::~FiniteElementSolver()
{

}

template <unsigned DIM>
boost::shared_ptr<FiniteElementSolver<DIM> > FiniteElementSolver<DIM>::Create()
{
    MAKE_PTR(FiniteElementSolver<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::Setup()
{

}

template<unsigned DIM>
void FiniteElementSolver<DIM>::Update()
{

}

template<unsigned DIM>
void FiniteElementSolver<DIM>::UpdateCellData()
{

}

template<unsigned DIM>
void FiniteElementSolver<DIM>::ReadSolution()
{
    mFeVtkSolution = vtkSmartPointer<vtkUnstructuredGrid>::New();
    VtkMeshReader<DIM,DIM> mesh_reader(this->mpOutputFileHandler->GetOutputDirectoryFullPath() + this->mFilename + ".vtu");
    vtkUnstructuredGrid* p_grid = mesh_reader.OutputMeshAsVtkUnstructuredGrid();
    mFeVtkSolution->DeepCopy(p_grid);
}

template<unsigned DIM>
std::vector<double> FiniteElementSolver<DIM>::GetSolutionAtGridPoints(boost::shared_ptr<RegularGrid<DIM, DIM> > pGrid)
{
    return GetSolutionAtPoints(pGrid->GetLocations());
}

template<unsigned DIM>
std::vector<double> FiniteElementSolver<DIM>::GetSolutionAtPoints(std::vector<c_vector<double, DIM> > samplePoints)
{
    if(!mFeVtkSolution)
    {
        ReadSolution();
    }
    
    // Sample the field at these locations
    vtkSmartPointer<vtkPolyData> p_polydata = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> p_points = vtkSmartPointer<vtkPoints>::New();
    p_points->SetNumberOfPoints(samplePoints.size());
    for(unsigned idx=0; idx< samplePoints.size(); idx++)
    {
        if(DIM==3)
        {
            p_points->SetPoint(idx, samplePoints[idx][0], samplePoints[idx][1], samplePoints[idx][2]);
        }
        else
        {
            p_points->SetPoint(idx, samplePoints[idx][0], samplePoints[idx][1], 0.0);
        }
    }
    p_polydata->SetPoints(p_points);

    vtkSmartPointer<vtkProbeFilter> p_probe_filter = vtkSmartPointer<vtkProbeFilter>::New();
    p_probe_filter->SetInputData(p_polydata);
    p_probe_filter->SetSourceData(mFeVtkSolution);
    p_probe_filter->Update();

    vtkSmartPointer<vtkPointData> p_point_data = p_probe_filter->GetOutput()->GetPointData();
    unsigned num_points = p_point_data->GetArray(this->mLabel.c_str())->GetNumberOfTuples();
    std::vector<double> results(num_points);
    for(unsigned idx=0; idx<num_points; idx++)
    {
        results[idx] = p_point_data->GetArray(this->mLabel.c_str())->GetTuple1(idx);
    }
    return results;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::SetMesh(boost::shared_ptr<HybridMesh<DIM, DIM> > pMesh)
{
    mpMesh = pMesh;
}

template<unsigned DIM>
std::vector<double> FiniteElementSolver<DIM>::GetNodalSolution()
{
    return mFeSolution;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::SetGuess(std::vector<double> guess)
{
    mGuess = guess;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::SetUseSimpleNetonSolver(bool useNewton)
{
    mUseNewton = useNewton;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::SetUseLinearSolveForGuess(bool useLinearSolve)
{
    mUseLinearSolveForGuess = useLinearSolve;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::Solve()
{
    boost::shared_ptr<BoundaryConditionsContainer<DIM, DIM, 1> > p_bcc =
            boost::shared_ptr<BoundaryConditionsContainer<DIM, DIM, 1> >(new BoundaryConditionsContainer<DIM, DIM, 1> );

    for(unsigned idx=0; idx<this->mBoundaryConditions.size(); idx++)
    {
        this->mBoundaryConditions[idx]->SetMesh(mpMesh);
        this->mBoundaryConditions[idx]->UpdateBoundaryConditionContainer(p_bcc);
    }

    // Do the solve
    // Check the type of pde
    if(this->mpPde and !this->mpNonLinearPde)
    {
        this->mpPde->SetUseRegularGrid(false);
        this->mpPde->SetMesh(mpMesh);
        this->mpPde->UpdateDiscreteSourceStrengths();

        SimpleLinearEllipticSolver<DIM, DIM> static_solver(mpMesh.get(), this->mpPde.get(), p_bcc.get());
        ReplicatableVector static_solution_repl(static_solver.Solve());
        mFeSolution = std::vector<double>(static_solution_repl.GetSize());
        for(unsigned idx = 0; idx < static_solution_repl.GetSize(); idx++)
        {
            mFeSolution[idx]= static_solution_repl[idx];
        }
    }
    else if(this->mpNonLinearPde)
    {
        this->mpNonLinearPde->SetUseRegularGrid(false);
        this->mpNonLinearPde->SetMesh(mpMesh);
        this->mpNonLinearPde->UpdateDiscreteSourceStrengths();

        if (this->mpPde)
        {
            this->mpPde->SetUseRegularGrid(false);
            this->mpPde->SetMesh(mpMesh);
            this->mpPde->UpdateDiscreteSourceStrengths();

            SimpleLinearEllipticSolver<DIM, DIM> static_solver(mpMesh.get(), this->mpPde.get(), p_bcc.get());
            ReplicatableVector static_solution_repl(static_solver.Solve());

            std::cout << "Completed Linear Solve: Starting Non-Linear Solve" << std::endl;
            std::vector<double> solution = std::vector<double>(static_solution_repl.GetSize());
            for(unsigned idx = 0; idx < static_solution_repl.GetSize(); idx++)
            {
                solution[idx]= static_solution_repl[idx];
                if(solution[idx]<0.0)
                {
                    solution[idx] = 0.0;
                }
                if(solution[idx]<1.0)
                {
                    solution[idx] = solution[idx] / 1.0;
                }
            }

            Vec initial_guess = PetscTools::CreateVec(mpMesh->GetNumNodes());
            for(unsigned idx=0; idx<solution.size();idx++)
            {
                PetscVecTools::SetElement(initial_guess, idx, solution[idx]);
            }
            PetscVecTools::Finalise(initial_guess);

            SimpleNonlinearEllipticSolver<DIM, DIM> solver(mpMesh.get(), this->mpNonLinearPde.get(), p_bcc.get());
            SimpleNewtonNonlinearSolver newton_solver;
            if(mUseNewton)
            {
                solver.SetNonlinearSolver(&newton_solver);
                newton_solver.SetTolerance(1e-5);
                newton_solver.SetWriteStats();
            }

            ReplicatableVector solution_repl(solver.Solve(initial_guess));
            mFeSolution = std::vector<double>(solution_repl.GetSize());
            for(unsigned idx = 0; idx < solution_repl.GetSize(); idx++)
            {
                mFeSolution[idx]= solution_repl[idx];
            }
            PetscTools::Destroy(initial_guess);

        }
        else
        {
            Vec initial_guess = PetscTools::CreateAndSetVec(mpMesh->GetNumNodes(), this->mBoundaryConditions[0]->GetValue());
            SimpleNonlinearEllipticSolver<DIM, DIM> solver(mpMesh.get(), this->mpNonLinearPde.get(), p_bcc.get());
            SimpleNewtonNonlinearSolver newton_solver;
            if(mUseNewton)
            {
                solver.SetNonlinearSolver(&newton_solver);
                newton_solver.SetTolerance(1e-5);
                newton_solver.SetWriteStats();
            }

            ReplicatableVector static_solution_repl(solver.Solve(initial_guess));

            mFeSolution = std::vector<double>(static_solution_repl.GetSize());
            for(unsigned idx = 0; idx < static_solution_repl.GetSize(); idx++)
            {
                mFeSolution[idx]= static_solution_repl[idx];
            }
            PetscTools::Destroy(initial_guess);
        }
    }
    else
    {
        EXCEPTION("PDE Type could not be identified, did you set a PDE?");
    }

    if(this->mWriteSolution)
    {
        this->Write();
    }
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::Write()
{
    // Write the output
    if(this->mFilename.empty())
    {
        this->mFilename = "solution";
    }

    if(!this->mpOutputFileHandler)
    {
        EXCEPTION("Output file handler not set");
    }

    VtkMeshWriter <DIM, DIM> mesh_writer(this->mpOutputFileHandler->GetRelativePath(), this->mFilename, false);
    if(mFeSolution.size() > 0)
    {
        mesh_writer.AddPointData(this->mLabel, mFeSolution);
    }
    mesh_writer.WriteFilesUsingMesh(*mpMesh);
    ReadSolution();
}

// Explicit instantiation
template class FiniteElementSolver<2> ;
template class FiniteElementSolver<3> ;
