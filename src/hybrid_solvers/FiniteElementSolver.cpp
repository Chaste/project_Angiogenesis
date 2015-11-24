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
#include "Debug.hpp"


template<unsigned DIM>
FiniteElementSolver<DIM>::FiniteElementSolver()
    : AbstractHybridSolver<DIM>(),
      mFeSolution(),
      mFeVtkSolution(),
      mpMesh(),
      mDomainBounds(),
      mDomainOutOfBoundsValue(0.0),
      mpNetwork(),
      mNetworkBounds(0.0)
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
std::vector<double> FiniteElementSolver<DIM>::GetSolutionOnRegularGrid(boost::shared_ptr<RegularGrid<DIM, DIM> > pGrid, bool useVtkSampling)
{
    return GetSolutionAtPoints(pGrid->GetLocations(), useVtkSampling);
}

template<unsigned DIM>
std::vector<double> FiniteElementSolver<DIM>::GetSolutionAtPointsUsingVtk(std::vector<c_vector<double, DIM> > samplePoints)
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
    p_probe_filter->SetInput(p_polydata);
    p_probe_filter->SetSource(mFeVtkSolution);
    p_probe_filter->Update();

    vtkSmartPointer<vtkPointData> p_point_data = p_probe_filter->GetOutput()->GetPointData();
    unsigned num_points = p_point_data->GetArray(this->mpPde->GetVariableName().c_str())->GetNumberOfTuples();
    std::vector<double> results(num_points);
    for(unsigned idx=0; idx<num_points; idx++)
    {
        results[idx] = p_point_data->GetArray(this->mpPde->GetVariableName().c_str())->GetTuple1(idx);
    }
    return results;
}

template<unsigned DIM>
vtkSmartPointer<vtkImageData> FiniteElementSolver<DIM>::GetVtkRegularGridSolution(boost::shared_ptr<RegularGrid<DIM, DIM> > pGrid, bool useVtkSampling)
{
    vtkSmartPointer<vtkImageData> p_vtk_grid = vtkSmartPointer<vtkImageData>::New();
    p_vtk_grid->SetSpacing(pGrid->GetSpacing(), pGrid->GetSpacing(), pGrid->GetSpacing());
    std::vector<unsigned> extents = pGrid->GetExtents();
    if(DIM==3)
    {
        p_vtk_grid->SetDimensions(extents[0], extents[1], extents[2]);
    }
    else
    {
        p_vtk_grid->SetDimensions(extents[0], extents[1], 1);
    }

    c_vector<double, DIM> orign = pGrid->GetOrigin();
    if(DIM==3)
    {
        p_vtk_grid->SetOrigin(orign[0], orign[1], orign[2]);
    }
    else
    {
        p_vtk_grid->SetOrigin(orign[0], orign[1], 0.0);
    }

    std::vector<double> result = GetSolutionOnRegularGrid(pGrid);

    vtkSmartPointer<vtkDoubleArray> pPointData = vtkSmartPointer<vtkDoubleArray>::New();
    pPointData->SetNumberOfComponents(1);
    pPointData->SetNumberOfTuples(result.size());
    pPointData->SetName(this->mpPde->GetVariableName().c_str());
    for (unsigned i = 0; i < result.size(); i++)
    {
        pPointData->SetValue(i, result[i]);
    }
    p_vtk_grid->GetPointData()->AddArray(pPointData);
    return p_vtk_grid;
}

template<unsigned DIM>
std::vector<double> FiniteElementSolver<DIM>::GetSolutionAtPoints(std::vector<c_vector<double, DIM> > samplePoints, bool useVtkSampling)
{
    if(useVtkSampling)
    {
        return GetSolutionAtPointsUsingVtk(samplePoints);
    }
    else
    {
        std::vector<double> sampled_solution(samplePoints.size(), 0.0);
        double tolerance = 1.e-3;

        for(unsigned idx=0; idx<samplePoints.size(); idx++)
        {
            double solution_at_point = 0.0;

            if(mDomainBounds.size() > 0)
            {
                // TODO fill in
            }
            if(mNetworkBounds)
            {
                std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = this->mpNetwork->GetVesselSegments();
                for (unsigned jdx = 0; jdx <  segments.size(); jdx++)
                {
                    if (segments[jdx]->GetDistance(samplePoints[idx]) <= segments[jdx]->GetRadius() + tolerance)
                    {
                        sampled_solution[idx] = mNetworkBounds;
                        continue;
                    }
                }
            }

            unsigned elem_index = mpMesh->GetContainingElementIndex(ChastePoint<DIM>(samplePoints[idx]));
            Element<DIM,DIM>* p_containing_element = mpMesh->GetElement(elem_index);
            c_vector<double,DIM+1> weights = p_containing_element->CalculateInterpolationWeights(samplePoints[idx]);
            for (unsigned i=0; i<DIM+1; i++)
            {
                double nodal_value = mFeSolution[p_containing_element->GetNodeGlobalIndex(i)];
                solution_at_point += nodal_value * weights(i);
            }
            sampled_solution[idx] = solution_at_point;
        }
        return sampled_solution;
    }
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::SetMesh(boost::shared_ptr<HybridMesh<DIM, DIM> > pMesh)
{
    mpMesh = pMesh;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::SetSamplingDomainBounds(std::vector<double> domainBounds, double domainOutOfBoundsValue)
{
    mDomainBounds = domainBounds;
    mDomainOutOfBoundsValue = domainOutOfBoundsValue;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::SetSamplingNetworkBounds(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork, double networkBounds)
{
    mpNetwork = pNetwork;
    mNetworkBounds = networkBounds;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::Solve()
{
    this->mpPde->SetUseRegularGrid(false);
    this->mpPde->SetMesh(mpMesh);
    this->mpPde->UpdateDiscreteSourceStrengths();

    boost::shared_ptr<BoundaryConditionsContainer<DIM, DIM, 1> > p_bcc =
            boost::shared_ptr<BoundaryConditionsContainer<DIM, DIM, 1> >(new BoundaryConditionsContainer<DIM, DIM, 1> );

    for(unsigned idx=0; idx<this->mBoundaryConditions.size(); idx++)
    {
        this->mBoundaryConditions[idx]->SetMesh(mpMesh);
        this->mBoundaryConditions[idx]->UpdateBoundaryConditionContainer(p_bcc);
    }

    // Do the solve
    SimpleLinearEllipticSolver<DIM, DIM> static_solver(mpMesh.get(), this->mpPde.get(), p_bcc.get());
    ReplicatableVector static_solution_repl(static_solver.Solve());
    mFeSolution = std::vector<double>(static_solution_repl.GetSize());
    for(unsigned idx = 0; idx < static_solution_repl.GetSize(); idx++)
    {
        mFeSolution[idx]= static_solution_repl[idx];
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
    std::string fname;
    if(!this->mFilename.empty())
    {
        fname = this->mFilename;
    }
    else
    {
        fname = "solution";
    }

    if(!this->mpOutputFileHandler)
    {
        EXCEPTION("Output file handler not set");
    }
    VtkMeshWriter <DIM, DIM> mesh_writer(this->mpOutputFileHandler->GetRelativePath(), fname, false);
    if(mFeSolution.size() > 0)
    {
        mesh_writer.AddPointData(this->mpPde->GetVariableName(), mFeSolution);
    }
    mesh_writer.WriteFilesUsingMesh(*mpMesh);
    ReadSolution();
}

// Explicit instantiation
template class FiniteElementSolver<2> ;
template class FiniteElementSolver<3> ;
