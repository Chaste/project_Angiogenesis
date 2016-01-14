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

#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning for now (gcc4.3)
#include <vtkXMLImageDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkProbeFilter.h>

#include "AbstractRegularGridHybridSolver.hpp"

template<unsigned DIM>
AbstractRegularGridHybridSolver<DIM>::AbstractRegularGridHybridSolver()
    :   AbstractHybridSolver<DIM>(),
        mpVtkSolution(),
        mpRegularGrid()
{

}

template<unsigned DIM>
AbstractRegularGridHybridSolver<DIM>::~AbstractRegularGridHybridSolver()
{

}

template<unsigned DIM>
boost::shared_ptr<RegularGrid<DIM> > AbstractRegularGridHybridSolver<DIM>::GetGrid()
{
    if(!mpRegularGrid)
    {
        EXCEPTION("A regular grid has not been set.");
    }

    return mpRegularGrid;
}

template<unsigned DIM>
std::vector<double> AbstractRegularGridHybridSolver<DIM>::GetSolutionAtPoints(std::vector<c_vector<double, DIM> > samplePoints, bool useVtkSampling)
{
    if(!mpVtkSolution)
    {
        EXCEPTION("A VTK solution has not be set up, did you forget to call Setup prior to Solve.");
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
    p_probe_filter->SetSource(mpVtkSolution);
    p_probe_filter->Update();
    vtkSmartPointer<vtkPointData> p_point_data = p_probe_filter->GetOutput()->GetPointData();
    unsigned num_points = p_point_data->GetArray(this->mpPde->GetVariableName().c_str())->GetNumberOfTuples();
    for(unsigned idx=0; idx<num_points; idx++)
    {
        sampled_solution[idx] = p_point_data->GetArray(this->mpPde->GetVariableName().c_str())->GetTuple1(idx);
    }
    return sampled_solution;
}

template<unsigned DIM>
vtkSmartPointer<vtkImageData> AbstractRegularGridHybridSolver<DIM>::GetVtkSolution()
{
    if(!mpVtkSolution)
    {
        EXCEPTION("A VTK solution has not be set up, did you forget to call Setup prior to Solve.");
    }

    return mpVtkSolution;
}

template<unsigned DIM>
void AbstractRegularGridHybridSolver<DIM>::SetGrid(boost::shared_ptr<RegularGrid<DIM> > pGrid)
{
    mpRegularGrid = pGrid;
}

template<unsigned DIM>
bool AbstractRegularGridHybridSolver<DIM>::HasRegularGrid()
{
    return true;
}

template<unsigned DIM>
void AbstractRegularGridHybridSolver<DIM>::Setup()
{
    if(!mpRegularGrid)
    {
        EXCEPTION("Regular grid hybrid solvers need a grid before Setup can be called.");
    }

    // Set up the VTK solution
    mpVtkSolution = vtkSmartPointer<vtkImageData>::New();

    if(DIM==3)
    {
        mpVtkSolution->SetDimensions(mpRegularGrid->GetExtents()[0], mpRegularGrid->GetExtents()[1], mpRegularGrid->GetExtents()[2]);
    }
    else
    {
        mpVtkSolution->SetDimensions(mpRegularGrid->GetExtents()[0], mpRegularGrid->GetExtents()[1], 1);
    }

    mpVtkSolution->SetSpacing(mpRegularGrid->GetSpacing(), mpRegularGrid->GetSpacing(), mpRegularGrid->GetSpacing());

    if(DIM==3)
    {
        mpVtkSolution->SetOrigin(mpRegularGrid->GetOrigin()[0], mpRegularGrid->GetOrigin()[1], mpRegularGrid->GetOrigin()[2]);
    }
    else
    {
        mpVtkSolution->SetOrigin(mpRegularGrid->GetOrigin()[0], mpRegularGrid->GetOrigin()[1], 0.0);
    }
}

template<unsigned DIM>
void AbstractRegularGridHybridSolver<DIM>::UpdateSolution(std::map<std::string, std::vector<double> >& data)
{
    if(!mpVtkSolution)
    {
        EXCEPTION("A VTK solution has not be set up, did you forget to call Setup prior to Solve.");
    }
    std::map<std::string, std::vector<double> >::iterator iter;
    for (iter = data.begin(); iter != data.end(); ++iter)
    {
        vtkSmartPointer<vtkDoubleArray> pPointData = vtkSmartPointer<vtkDoubleArray>::New();
        pPointData->SetNumberOfComponents(1);
        pPointData->SetNumberOfTuples(iter->second.size());
        pPointData->SetName(iter->first.c_str());
        this->mPointSolution = std::vector<double>(iter->second.size());
        for (unsigned i = 0; i < iter->second.size(); i++)
        {
            pPointData->SetValue(i, data[iter->first][i]);
            this->mPointSolution[i] = data[iter->first][i];
        }
        mpVtkSolution->GetPointData()->AddArray(pPointData);
    }
}

template<unsigned DIM>
void AbstractRegularGridHybridSolver<DIM>::UpdateCellData()
{
    mpRegularGrid->SetCellPopulation(*(this->mpCellPopulation));
    std::vector<std::vector<CellPtr> > point_cell_map = mpRegularGrid->GetPointCellMap();
    for(unsigned idx=0; idx<point_cell_map.size(); idx++)
    {
        for(unsigned jdx=0; jdx<point_cell_map[idx].size(); jdx++)
        {
            point_cell_map[idx][jdx]->GetCellData()->SetItem(this->mLabel, this->mPointSolution[idx]);
        }
    }
}

template<unsigned DIM>
void AbstractRegularGridHybridSolver<DIM>::Write()
{
    if(!mpVtkSolution)
    {
        EXCEPTION("A VTK solution has not be set up, did you forget to call Setup prior to Solve.");
    }

    if(!this->mpOutputFileHandler)
    {
        EXCEPTION("An output file handler has not been set for the hybrid solver.");
    }

    vtkSmartPointer<vtkXMLImageDataWriter> pImageDataWriter = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    if(!this->mFilename.empty())
    {
        pImageDataWriter->SetFileName((this->mpOutputFileHandler->GetOutputDirectoryFullPath() + "/" + this->mFilename).c_str());
    }
    else
    {
        pImageDataWriter->SetFileName((this->mpOutputFileHandler->GetOutputDirectoryFullPath() + "/solution.vti").c_str());
    }

    pImageDataWriter->SetInput(mpVtkSolution);
    pImageDataWriter->Update();
    pImageDataWriter->Write();
}

// Explicit instantiation
template class AbstractRegularGridHybridSolver<2> ;
template class AbstractRegularGridHybridSolver<3> ;
