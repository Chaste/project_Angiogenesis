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
#include "SmartPointers.hpp"
#include "UblasVectorInclude.hpp"

#include "AbstractRegularGridHybridSolver.hpp"

template<unsigned DIM>
AbstractRegularGridHybridSolver<DIM>::AbstractRegularGridHybridSolver()
    :   AbstractHybridSolver<DIM>(),
        mpSolution(),
        mGridSize(1.0),
        mOrigin(zero_vector<double>(DIM)),
        mExtents(std::vector<unsigned>(3, 1))
{

}

template<unsigned DIM>
AbstractRegularGridHybridSolver<DIM>::~AbstractRegularGridHybridSolver()
{

}

template<unsigned DIM>
std::vector<double> AbstractRegularGridHybridSolver<DIM>::GetSolutionAtPoints(std::vector<c_vector<double, DIM> > samplePoints)
{
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
    p_probe_filter->SetSource(GetSolution());
    p_probe_filter->Update();
    vtkSmartPointer<vtkPointData> p_point_data = p_probe_filter->GetOutput()->GetPointData();
    unsigned num_points = p_point_data->GetArray("Default")->GetNumberOfTuples();
    for(unsigned idx=0; idx<num_points; idx++)
    {
        sampled_solution[idx] = p_point_data->GetArray("Default")->GetTuple1(idx);
    }
    return sampled_solution;
}

template<unsigned DIM>
unsigned AbstractRegularGridHybridSolver<DIM>::GetGridIndex(unsigned x_index, unsigned y_index, unsigned z_index)
{
    return x_index + mExtents[0] * y_index + mExtents[0] * mExtents[1] * z_index;
}

template<unsigned DIM>
c_vector<double, DIM> AbstractRegularGridHybridSolver<DIM>::GetLocation(unsigned x_index, unsigned y_index, unsigned z_index)
{
    c_vector<double, DIM> location;
    location[0] = double(x_index) * mGridSize + mOrigin[0];
    location[1] = double(y_index) * mGridSize + mOrigin[1];
    if(DIM==3)
    {
        location[2] = double(z_index) * mGridSize + mOrigin[2];
    }
    return location;
}

template<unsigned DIM>
void AbstractRegularGridHybridSolver<DIM>::SetOrigin(c_vector<double, DIM> origin)
{
    mOrigin = origin;
}

template<unsigned DIM>
vtkSmartPointer<vtkImageData> AbstractRegularGridHybridSolver<DIM>::GetSolution()
{
    return mpSolution;
}

template<unsigned DIM>
bool AbstractRegularGridHybridSolver<DIM>::IsOnBoundary(unsigned x_index, unsigned y_index, unsigned z_index)
{
    if(x_index == 0 || x_index == mExtents[0] - 1)
    {
        return true;
    }
    if(y_index == 0 || y_index == mExtents[1] - 1)
    {
        return true;
    }
    if(DIM==3)
    {
        if(z_index == 0 || z_index == mExtents[2] - 1)
        {
            return true;
        }
    }
    return false;
}

template<unsigned DIM>
void AbstractRegularGridHybridSolver<DIM>::SetExtents(std::vector<unsigned> extents)
{
    mExtents = extents;
}

template<unsigned DIM>
void AbstractRegularGridHybridSolver<DIM>::SetExtents(boost::shared_ptr<Part<DIM> > pPart, double gridSize)
{
    mGridSize = gridSize;
    c_vector<double, 2*DIM> spatial_extents = pPart->GetBoundingBox();
    mExtents = std::vector<unsigned>(DIM);
    mExtents[0] = unsigned((spatial_extents[1] - spatial_extents[0]) / mGridSize) + 1u;
    mExtents[1] = unsigned((spatial_extents[3] - spatial_extents[2]) / mGridSize) + 1u;
    if(DIM==3)
    {
        mExtents[2] = unsigned((spatial_extents[5] - spatial_extents[4]) / mGridSize) + 1u;
    }
    mOrigin[0] = spatial_extents[0];
    mOrigin[1] = spatial_extents[2];
    if(DIM==3)
    {
        mOrigin[2] = spatial_extents[4];
    }
}

template<unsigned DIM>
void AbstractRegularGridHybridSolver<DIM>::SetGridSize(double gridSize)
{
    mGridSize = gridSize;
}


template<unsigned DIM>
void AbstractRegularGridHybridSolver<DIM>::UpdateSolution(std::map<std::string, std::vector<double> >& data)
{
    mpSolution = vtkSmartPointer<vtkImageData>::New();
    if(DIM==3)
    {
        mpSolution->SetDimensions(mExtents[0], mExtents[1], mExtents[2]);
    }
    else
    {
        mpSolution->SetDimensions(mExtents[0], mExtents[1], 1);
    }

    mpSolution->SetSpacing(mGridSize, mGridSize, mGridSize);

    if(DIM==3)
    {
        mpSolution->SetOrigin(mOrigin[0], mOrigin[1], mOrigin[2]);
    }
    else
    {
        mpSolution->SetOrigin(mOrigin[0], mOrigin[1], 0.0);
    }

    std::map<std::string, std::vector<double> >::iterator iter;
    for (iter = data.begin(); iter != data.end(); ++iter)
    {
        vtkSmartPointer<vtkDoubleArray> pPointData = vtkSmartPointer<vtkDoubleArray>::New();
        pPointData->SetNumberOfComponents(1);
        pPointData->SetNumberOfTuples(iter->second.size());
        pPointData->SetName(iter->first.c_str());

        for (unsigned i = 0; i < iter->second.size(); i++)
        {
            pPointData->SetValue(i, data[iter->first][i]);
        }
        mpSolution->GetPointData()->AddArray(pPointData);
    }
}

template<unsigned DIM>
void AbstractRegularGridHybridSolver<DIM>::Write()
{
    if (!this->mWorkingDirectory.empty())
    {
        if (this->mpNetwork)
        {
            // Write the vessel network data
            this->mpNetwork->Write((this->mWorkingDirectory + "/vessel_network.vtp").c_str());
        }
        if (this->mpCellPopulation)
        {
           this->mpCellPopulation->Write((this->mWorkingDirectory + "/cell_population.vtp").c_str());
        }

        // Write the PDE solution
        vtkSmartPointer<vtkXMLImageDataWriter> pImageDataWriter = vtkSmartPointer<vtkXMLImageDataWriter>::New();
        if(!this->mFilename.empty())
        {
            pImageDataWriter->SetFileName((this->mWorkingDirectory + "/" + this->mFilename).c_str());
        }
        else
        {
            pImageDataWriter->SetFileName((this->mWorkingDirectory + "/solution.vti").c_str());
        }
        pImageDataWriter->SetInput(mpSolution);
        pImageDataWriter->Update();
        pImageDataWriter->Write();
    }
}

// Explicit instantiation
template class AbstractRegularGridHybridSolver<2> ;
template class AbstractRegularGridHybridSolver<3> ;
