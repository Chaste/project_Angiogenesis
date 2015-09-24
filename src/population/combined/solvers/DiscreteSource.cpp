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

#include "Facet.hpp"
#include "DiscreteSource.hpp"
#include "SimpleCell.hpp"
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkProbeFilter.h>
#include <vtkLineSource.h>

template<unsigned DIM>
DiscreteSource<DIM>::DiscreteSource()
    :   mpNetwork(),
        mpCellPopulation(),
        mpDomain(),
        mpSolution(),
        mPoints(),
        mType(SourceType::POINT),
        mSourceStrength(SourceStrength::PRESCRIBED),
        mLabel("Default"),
        mValue(0.0),
        mIsLinearInSolution(false)
{

}

template<unsigned DIM>
DiscreteSource<DIM>::~DiscreteSource()
{

}

template<unsigned DIM>
SourceType::Value DiscreteSource<DIM>::GetType()
{
    return mType;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetSolution(vtkSmartPointer<vtkImageData>  pSolution)
{
    mpSolution = pSolution;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetCellPopulation(boost::shared_ptr<SimpleCellPopulation<DIM> > pCellPopulation)
{
    mpCellPopulation = pCellPopulation;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetDomain(boost::shared_ptr<Part<DIM> > pDomain)
{
    mpDomain = pDomain;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetPoint(c_vector<double, DIM> point)
{
    mPoints.push_back(point);
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetPoints(std::vector<c_vector<double, DIM> > points)
{
    mPoints = points;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetType(SourceType::Value boundaryType)
{
    mType = boundaryType;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetSource(SourceStrength::Value boundarySource)
{
    mSourceStrength = boundarySource;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetLabelName(const std::string& label)
{
    mLabel = label;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetValue(double value)
{
    mValue = value;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetIsLinearInSolution(bool isLinear)
{
    mIsLinearInSolution = isLinear;
}

template<unsigned DIM>
std::pair<bool, double> DiscreteSource<DIM>::GetValue(c_vector<double, DIM> location, double tolerance)
{
    // Check the source type
    if(mType == SourceType::POINT)
    {
        if(mPoints.size()==0)
        {
            EXCEPTION("A point is required for this type of source");
        }
        else
        {
            if(norm_2(location-mPoints[0]) < tolerance)
            {
                return std::pair<bool, double>(true, mValue);
            }
            else
            {
                return std::pair<bool, double>(false, mValue);
            }
        }
    }
    else if(mType == SourceType::MULTI_POINT)
    {
        if(mPoints.size()==0)
        {
            EXCEPTION("A point is required for this type of source");
        }
        else
        {
            for(unsigned idx=0; idx<mPoints.size(); idx++)
            {
                if(norm_2(location-mPoints[0]) < tolerance)
                {
                    return std::pair<bool, double>(true, mValue);
                }
            }
            return std::pair<bool, double>(false, mValue);
        }
    }

    else if(mType == SourceType::VESSEL_LINE)
    {
        if(!mpNetwork)
        {
            EXCEPTION("A vessel network is required for this type of source");
        }
        else
        {
            std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = this->mpNetwork->GetVesselSegments();
            for (unsigned idx = 0; idx <  segments.size(); idx++)
            {
                if (segments[idx]->GetDistance(location) <= tolerance)
                {
                    if(SourceStrength::PRESCRIBED)
                    {
                        return std::pair<bool, double>(true, mValue);
                    }
                    else
                    {
                        return std::pair<bool, double>(true, segments[idx]->template GetData<double>(mLabel));
                    }
                }
            }
            return std::pair<bool, double>(false, mValue);
        }
    }

    else if(mType == SourceType::CELL_POINT)
    {
        if(!mpCellPopulation)
        {
            EXCEPTION("A simple cell population is required for this type of source");
        }

        std::vector<boost::shared_ptr<SimpleCell<DIM> > > cells = mpCellPopulation->GetCells();
        for(unsigned idx=0; idx<cells.size(); idx++)
        {
            if (norm_2(cells[idx]->rGetLocation()-location)<tolerance)
            {
                return std::pair<bool, double>(true, mValue);
            }
        }
        return std::pair<bool, double>(false, mValue);
    }

    else if(mType == SourceType::SOLUTION)
    {
        if(!mpSolution)
        {
            EXCEPTION("A previous solution is required for this type of source");
        }
        std::vector<c_vector<double, DIM> > samplePoints(1, location);
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
        p_probe_filter->SetSource(mpSolution);
        p_probe_filter->Update();
        vtkSmartPointer<vtkPointData> p_point_data = p_probe_filter->GetOutput()->GetPointData();
        unsigned num_points = p_point_data->GetArray(mLabel.c_str())->GetNumberOfTuples();
        for(unsigned idx=0; idx<num_points; idx++)
        {
            sampled_solution[idx] = p_point_data->GetArray(mLabel.c_str())->GetTuple1(idx);
        }
        return std::pair<bool, double>(true, sampled_solution[0]);
    }
    return std::pair<bool, double>(false, mValue);
}

template<unsigned DIM>
std::vector<std::pair<bool, double> > DiscreteSource<DIM>::GetValues(std::vector<c_vector<double, DIM> > locations, double tolerance)
{
    std::vector<std::pair<bool, double> > result;
    for(unsigned idx=0; idx<locations.size();idx++)
    {
        result.push_back(GetValue(locations[idx],tolerance));
    }
    return result;
}

// Explicit instantiation
template class DiscreteSource<2>;
template class DiscreteSource<3>;
