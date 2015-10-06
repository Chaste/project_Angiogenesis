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
//Requires  "sudo aptitude install libvtk5-dev" or similar
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkProbeFilter.h>
#include "VtkMeshReader.hpp"
#include "VtkMeshWriter.hpp"
#include "ConstBoundaryCondition.hpp"
#include "FiniteElementSolver.hpp"
#include "VesselSurfaceGenerator.hpp"
#include "PlcMesh.hpp"
#include "CaVesselSegment.hpp"
#include "CaVascularNetwork.hpp"
#include "Debug.hpp"

template<unsigned DIM>
FiniteElementSolver<DIM>::FiniteElementSolver()
    : AbstractHybridSolver<DIM>(),
      mpDomain(),
      mGridSize(100.0),
      mMeshWriterPath(),
      mFeSolution(),
      mVesselRepresentation(VesselRepresentation::LINE),
      mpMesh()
{

}

template<unsigned DIM>
FiniteElementSolver<DIM>::~FiniteElementSolver()
{

}

template<unsigned DIM>
void FiniteElementSolver<DIM>::SetMeshWriterPath(std::string path)
{
    mMeshWriterPath = path;
}

template <unsigned DIM>
boost::shared_ptr<FiniteElementSolver<DIM> > FiniteElementSolver<DIM>::Create()
{
    MAKE_PTR(FiniteElementSolver<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::ReadSolution()
{
    mFeSolution = vtkSmartPointer<vtkUnstructuredGrid>::New();
    VtkMeshReader<DIM,DIM> mesh_reader(this->mWorkingDirectory + this->mFilename + ".vtu");
    vtkUnstructuredGrid* p_grid = mesh_reader.OutputMeshAsVtkUnstructuredGrid();
    mFeSolution->DeepCopy(p_grid);
}

template<unsigned DIM>
std::vector<double> FiniteElementSolver<DIM>::GetSolutionAtPoints(std::vector<c_vector<double, DIM> > samplePoints,
                                                                              const std::string& rSpeciesLabel)
{
    // Read the vtk mesh back in and sample it at the requested points
    if(!mFeSolution)
    {
        ReadSolution();
    }
    std::vector<double> sampled_solution(samplePoints.size(), 0.0);

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
    p_probe_filter->SetInput(p_polydata);
    p_probe_filter->SetSource(mFeSolution);
    p_probe_filter->Update();
    vtkSmartPointer<vtkPointData> p_point_data = p_probe_filter->GetOutput()->GetPointData();
    unsigned num_points = p_point_data->GetArray(rSpeciesLabel.c_str())->GetNumberOfTuples();
    for(unsigned idx=0; idx<num_points; idx++)
    {
        sampled_solution[idx] = p_point_data->GetArray(rSpeciesLabel.c_str())->GetTuple1(idx);
    }
    return sampled_solution;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::SetVesselRepresentation(VesselRepresentation::Value vesselRepresentation)
{
    mVesselRepresentation = vesselRepresentation;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::SetDomain(boost::shared_ptr<Part<DIM> > pDomain)
{
    mpDomain = pDomain;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::SetMaxElementArea(double maxElementArea)
{
    mGridSize = maxElementArea;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::SetMesh(boost::shared_ptr<PlcMesh<DIM, DIM> > pMesh)
{
    mpMesh = pMesh;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::Solve(bool writeSolution)
{
    // If there is a vessel network add it to the domain
    if (this->mpNetwork)
    {
        mpDomain->AddVesselNetwork(this->mpNetwork, this->mVesselRepresentation == VesselRepresentation::SURFACE);
    }

    // Mesh the domain
    if(!mpMesh)
    {
        mpMesh = PlcMesh<DIM, DIM>::Create();
        mpMesh->GenerateFromPart(mpDomain, mGridSize);
    }

    // Apply the dirichlet boundary conditions
    double node_distance_tolerance = 1.e-3;
    BoundaryConditionsContainer<DIM, DIM, 1> bcc;

    for(unsigned idx=0; idx<this->mDirichletBoundaryConditions.size(); idx++)
    {
        bool use_boundry_nodes = false;
        if(this->mDirichletBoundaryConditions[idx]->GetType() == BoundaryConditionType::OUTER ||
                this->mDirichletBoundaryConditions[idx]->GetType() == BoundaryConditionType::FACET)
        {
            use_boundry_nodes = true;
        }
        else if(this->mVesselRepresentation == VesselRepresentation::SURFACE)
        {
            if(this->mDirichletBoundaryConditions[idx]->GetType() == BoundaryConditionType::VESSEL_LINE
                    || this->mDirichletBoundaryConditions[idx]->GetType() == BoundaryConditionType::VESSEL_VOLUME)
            {

                use_boundry_nodes = true;
            }
        }

        if(!use_boundry_nodes)
        {
            typename PlcMesh<DIM, DIM>::NodeIterator iter = mpMesh->GetNodeIteratorBegin();
            while (iter != mpMesh->GetNodeIteratorEnd())
            {
                std::pair<bool,double> result = this->mDirichletBoundaryConditions[idx]->GetValue((*iter).GetPoint().rGetLocation(), node_distance_tolerance);
                if(result.first)
                {
                    ConstBoundaryCondition<DIM>* p_fixed_boundary_condition = new ConstBoundaryCondition<DIM>(result.second);
                    bcc.AddDirichletBoundaryCondition(&(*iter), p_fixed_boundary_condition, 0, false);
                }
                ++iter;
            }
        }
        else
        {
            typename PlcMesh<DIM, DIM>::BoundaryNodeIterator iter = mpMesh->GetBoundaryNodeIteratorBegin();
            while (iter < mpMesh->GetBoundaryNodeIteratorEnd())
            {
                std::pair<bool,double> result = this->mDirichletBoundaryConditions[idx]->GetValue((*iter)->GetPoint().rGetLocation(), node_distance_tolerance);
                if(result.first)
                {
                    ConstBoundaryCondition<DIM>* p_fixed_boundary_condition = new ConstBoundaryCondition<DIM>(result.second);
                    bcc.AddDirichletBoundaryCondition(*iter, p_fixed_boundary_condition);
                }
                ++iter;
            }
        }
    }

    // Do the solve
    SimpleLinearEllipticSolver<DIM, DIM> static_solver(mpMesh.get(), this->mpPde.get(), &bcc);
    ReplicatableVector static_solution_repl(static_solver.Solve());
    std::vector<double> output;
    for (unsigned idx = 0; idx < static_solution_repl.GetSize(); idx++)
    {
        output.push_back(static_solution_repl[idx]);
    }

    if(writeSolution)
    {
        this->Write(output, mpMesh);
    }
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::Write(std::vector<double> output, boost::shared_ptr<PlcMesh<DIM, DIM> > p_mesh)
{
    // Write the output
    std::string fname;
    if(!this->mFilename.empty())
    {
        fname = this->mFilename;
    }
    else
    {
        fname = "solution";
    }

    VtkMeshWriter <DIM, DIM> mesh_writer(mMeshWriterPath, fname, false);
    if(output.size() > 0)
    {
        mesh_writer.AddPointData(this->mpPde->GetVariableName(), output);
    }
    mesh_writer.WriteFilesUsingMesh(*p_mesh);
    ReadSolution();
}

// Explicit instantiation
template class FiniteElementSolver<2> ;
template class FiniteElementSolver<3> ;
