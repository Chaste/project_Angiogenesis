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
#include "SmartPointers.hpp"
#include "UblasVectorInclude.hpp"

#include "AbstractRegularGridHybridSolver.hpp"

AbstractRegularGridHybridSolver::AbstractRegularGridHybridSolver()
    :   AbstractHybridSolver(),
        mpSolution(),
        mGridSize(1.0),
        mOrigin(zero_vector<double>(3)),
        mExtents(std::vector<unsigned>(3, 10))
{

}

AbstractRegularGridHybridSolver::~AbstractRegularGridHybridSolver()
{

}

unsigned AbstractRegularGridHybridSolver::GetGridIndex(unsigned x_index, unsigned y_index, unsigned z_index)
{
    return x_index + mExtents[0] * y_index + mExtents[0] * mExtents[1] * z_index;
}

c_vector<double, 3> AbstractRegularGridHybridSolver::GetLocation(unsigned x_index, unsigned y_index, unsigned z_index)
{
    c_vector<double, 3> location;
    location[0] = double(x_index) * mGridSize + mOrigin[0];
    location[1] = double(y_index) * mGridSize + mOrigin[1];
    location[2] = double(z_index) * mGridSize + mOrigin[2];
    return location;
}

void AbstractRegularGridHybridSolver::SetOrigin(c_vector<double, 3> origin)
{
    mOrigin = origin;
}

vtkSmartPointer<vtkImageData> AbstractRegularGridHybridSolver::GetSolution()
{
    return mpSolution;
}

bool AbstractRegularGridHybridSolver::IsOnBoundary(unsigned x_index, unsigned y_index, unsigned z_index)
{
    if(x_index == 0 || x_index == mExtents[0] - 1)
    {
        return true;
    }
    if(y_index == 0 || y_index == mExtents[1] - 1)
    {
        return true;
    }
    if(z_index == 0 || z_index == mExtents[2] - 1)
    {
        return true;
    }

    return false;
}

void AbstractRegularGridHybridSolver::SetExtents(std::vector<unsigned> extents)
{
    mExtents = extents;
}

void AbstractRegularGridHybridSolver::SetExtents(boost::shared_ptr<Part> pPart, double gridSize)
{
    mGridSize = gridSize;
    c_vector<double, 6> spatial_extents = pPart->GetBoundingBox();
    mExtents = std::vector<unsigned>(3);
    mExtents[0] = unsigned((spatial_extents[1] - spatial_extents[0]) / mGridSize) + 1u;
    mExtents[1] = unsigned((spatial_extents[3] - spatial_extents[2]) / mGridSize) + 1u;
    mExtents[2] = unsigned((spatial_extents[5] - spatial_extents[4]) / mGridSize) + 1u;
}

void AbstractRegularGridHybridSolver::SetGridSize(double gridSize)
{
    mGridSize = gridSize;
}

void SetOrigin(c_vector<double, 3> origin)
{

}

void AbstractRegularGridHybridSolver::UpdateSolution(std::map<std::string, std::vector<double> >& data)
{
    mpSolution = vtkSmartPointer<vtkImageData>::New();
    mpSolution->SetDimensions(mExtents[0], mExtents[1], mExtents[2]);
    mpSolution->SetSpacing(mGridSize, mGridSize, mGridSize);
    mpSolution->SetOrigin(mOrigin[0], mOrigin[1], mOrigin[2]);

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

void AbstractRegularGridHybridSolver::Write()
{
    if (!mWorkingDirectory.empty())
    {
        if (mpNetwork)
        {
            // Write the vessel network data
            mpNetwork->Write((mWorkingDirectory + "/vessel_network.vtp").c_str());
        }
        if (mpCellPopulation)
        {
           mpCellPopulation->Write((mWorkingDirectory + "/cell_population.vtp").c_str());
        }

        // Write the PDE solution
        vtkSmartPointer<vtkXMLImageDataWriter> pImageDataWriter = vtkSmartPointer<vtkXMLImageDataWriter>::New();
        pImageDataWriter->SetFileName((mWorkingDirectory + "/pde_solution.vti").c_str());
        pImageDataWriter->SetInput(mpSolution);
        pImageDataWriter->Update();
        pImageDataWriter->Write();
    }
}
