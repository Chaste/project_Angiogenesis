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
#include <vtkLineSource.h>
#include "SmartPointers.hpp"
#include "UblasVectorInclude.hpp"
#include "GeometryTools.hpp"

#include "AbstractRegularGridHybridSolver.hpp"

template<unsigned DIM>
AbstractRegularGridHybridSolver<DIM>::AbstractRegularGridHybridSolver()
    :   AbstractHybridSolver<DIM>(),
        mpRegularGridVtkSolution(),
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
    return mpRegularGrid;
}

template<unsigned DIM>
std::vector<double> AbstractRegularGridHybridSolver<DIM>::GetSolutionAtPoints(std::vector<c_vector<double, DIM> > samplePoints,
                                                                              const std::string& rSpeciesLabel)
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
    p_probe_filter->SetSource(GetVtkSolution());
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
std::vector<double> AbstractRegularGridHybridSolver<DIM>::GetSolutionOnLine(double sampleSpacing,
                                                                            c_vector<double, DIM> start_point,
                                                                            c_vector<double, DIM> end_point,
                                                                            const std::string& rSpeciesLabel)
{
    vtkSmartPointer<vtkLineSource> p_line = vtkSmartPointer<vtkLineSource>::New();
    p_line->SetPoint1(&start_point[0]);
    p_line->SetPoint2(&end_point[0]);
    p_line->SetResolution(sampleSpacing);
    vtkSmartPointer<vtkProbeFilter> p_probe_filter = vtkSmartPointer<vtkProbeFilter>::New();
    p_probe_filter->SetInputConnection(p_line->GetOutputPort());
    p_probe_filter->SetSource(mpRegularGridVtkSolution);
    p_probe_filter->Update();

    vtkSmartPointer<vtkPointData> p_point_data = p_probe_filter->GetOutput()->GetPointData();
    std::vector<double>line_data(p_point_data->GetNumberOfTuples());
    for(unsigned idx=0;idx<p_point_data->GetNumberOfTuples();idx++)
    {
        line_data[idx] = p_point_data->GetArray(rSpeciesLabel.c_str())->GetTuple1(idx);
    }
    return line_data;
}

template<unsigned DIM>
vtkSmartPointer<vtkImageData> AbstractRegularGridHybridSolver<DIM>::GetSolutionOnVolume(double sampleSpacing,
                                                          std::vector<unsigned> dimensions,
                                                          const std::string& rSpeciesLabel,
                                                          c_vector<double, DIM> origin)
{
    vtkSmartPointer<vtkImageData> p_sampling_grid = vtkSmartPointer<vtkImageData>::New();
    p_sampling_grid->SetSpacing(sampleSpacing, sampleSpacing, sampleSpacing);
    p_sampling_grid->SetOrigin(origin[0], origin[1], origin[2]);
    p_sampling_grid->SetDimensions(dimensions[0], dimensions[1], dimensions[2]);

    vtkSmartPointer<vtkProbeFilter> p_probe_filter = vtkSmartPointer<vtkProbeFilter>::New();
    p_probe_filter->SetInput(p_sampling_grid);
    p_probe_filter->SetSource(mpRegularGridVtkSolution);
    p_probe_filter->Update();
    return p_probe_filter->GetImageDataOutput();
}

template<unsigned DIM>
double AbstractRegularGridHybridSolver<DIM>::GetVolumeAverageSolution(const std::string& arrayName)
{
    unsigned num_points = mpRegularGridVtkSolution->GetPointData()->GetArray(0)->GetNumberOfTuples();
    double sum = 0.0;
    for(unsigned idx=0; idx<num_points; idx++)
    {
        sum+= mpRegularGridVtkSolution->GetPointData()->GetArray(arrayName.c_str())->GetTuple1(idx);
    }
    return sum / double(num_points);
}

template<unsigned DIM>
vtkSmartPointer<vtkImageData> AbstractRegularGridHybridSolver<DIM>::GetVtkSolution()
{
    return mpRegularGridVtkSolution;
}

template<unsigned DIM>
void AbstractRegularGridHybridSolver<DIM>::SetGridFromPart(boost::shared_ptr<Part<DIM> > pPart, double gridSize)
{
    mpRegularGrid = RegularGrid<DIM>::Create();
    mpRegularGrid->SetSpacing(gridSize);

    std::vector<unsigned> extents = std::vector<unsigned>(3);

    c_vector<double, 2*DIM> spatial_extents = pPart->GetBoundingBox();
    extents[0] = unsigned((spatial_extents[1] - spatial_extents[0]) / gridSize) + 1u;
    extents[1] = unsigned((spatial_extents[3] - spatial_extents[2]) / gridSize) + 1u;
    if(DIM==3)
    {
        extents[2] = unsigned((spatial_extents[5] - spatial_extents[4]) / gridSize) + 1u;
    }
    else
    {
        extents[2] = 1;
    }
    mpRegularGrid->SetExtents(extents);

    c_vector<double, DIM> origin;
    origin[0] = spatial_extents[0];
    origin[1] = spatial_extents[2];
    if(DIM==3)
    {
        origin[2] = spatial_extents[4];
    }
    mpRegularGrid->SetOrigin(origin);
}

template<unsigned DIM>
void AbstractRegularGridHybridSolver<DIM>::SetGrid(boost::shared_ptr<RegularGrid<DIM> > pGrid)
{
    mpRegularGrid = pGrid;
}

template<unsigned DIM>
void AbstractRegularGridHybridSolver<DIM>::UpdateSolution(std::map<std::string, std::vector<double> >& data)
{
    mpRegularGridVtkSolution = vtkSmartPointer<vtkImageData>::New();

    if(DIM==3)
    {
        mpRegularGridVtkSolution->SetDimensions(mpRegularGrid->GetExtents()[0], mpRegularGrid->GetExtents()[1], mpRegularGrid->GetExtents()[2]);
    }
    else
    {
        mpRegularGridVtkSolution->SetDimensions(mpRegularGrid->GetExtents()[0], mpRegularGrid->GetExtents()[1], 1);
    }

    mpRegularGridVtkSolution->SetSpacing(mpRegularGrid->GetSpacing(), mpRegularGrid->GetSpacing(), mpRegularGrid->GetSpacing());

    if(DIM==3)
    {
        mpRegularGridVtkSolution->SetOrigin(mpRegularGrid->GetOrigin()[0], mpRegularGrid->GetOrigin()[1], mpRegularGrid->GetOrigin()[2]);
    }
    else
    {
        mpRegularGridVtkSolution->SetOrigin(mpRegularGrid->GetOrigin()[0], mpRegularGrid->GetOrigin()[1], 0.0);
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
        mpRegularGridVtkSolution->GetPointData()->AddArray(pPointData);
    }
}

template<unsigned DIM>
void AbstractRegularGridHybridSolver<DIM>::Write()
{
    if (!this->mWorkingDirectory.empty())
    {
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

        pImageDataWriter->SetInput(mpRegularGridVtkSolution);
        pImageDataWriter->Update();
        pImageDataWriter->Write();
    }
}

template<unsigned DIM>
void AbstractRegularGridHybridSolver<DIM>::WriteHistograms(const std::string& arrayName, const std::string& fileName, double binSize, unsigned numberOfBins)
{
    std::vector<unsigned> bins(numberOfBins, 0);

    // populate the bins
    unsigned num_points = mpRegularGridVtkSolution->GetPointData()->GetArray(arrayName.c_str())->GetNumberOfTuples();
    for(unsigned idx=0; idx<num_points; idx++)
    {
        unsigned bin_label = std::floor(mpRegularGridVtkSolution->GetPointData()->GetArray(arrayName.c_str())->GetTuple1(idx)/ binSize);
        if(bin_label > numberOfBins)
        {
            bin_label = numberOfBins;
        }
        bins[bin_label]++;
    }

    std::ofstream output_file(fileName.c_str());
    output_file << "Bin , Value \n";
    if (output_file.is_open())
    {
        for(unsigned idx=0; idx< numberOfBins;idx++)
        {
            output_file << double(idx) * binSize << " , "<< bins[idx] << "\n";
        }
        output_file.close();
    }
}

// Explicit instantiation
template class AbstractRegularGridHybridSolver<2> ;
template class AbstractRegularGridHybridSolver<3> ;
