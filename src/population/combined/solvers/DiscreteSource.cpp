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
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkProbeFilter.h>
#include <vtkLineSource.h>
#include "Debug.hpp"

template<unsigned DIM>
DiscreteSource<DIM>::DiscreteSource()
    :   mpNetwork(),
        mpCellPopulation(),
        mpSolution(),
        mPoints(),
        mType(SourceType::POINT),
        mSourceStrength(SourceStrength::PRESCRIBED),
        mLabel("Default"),
        mValue(0.0),
        mIsLinearInSolution(false),
        mMutationSpecificConsumptionRateMap()
{

}

template<unsigned DIM>
DiscreteSource<DIM>::~DiscreteSource()
{

}

template<unsigned DIM>
boost::shared_ptr<DiscreteSource<DIM> > DiscreteSource<DIM>::Create()
{
    MAKE_PTR(DiscreteSource<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
SourceType::Value DiscreteSource<DIM>::GetType()
{
    return mType;
}

template<unsigned DIM>
std::vector<double> DiscreteSource<DIM>::GetValues(std::vector<c_vector<double, DIM> > locations, double tolerance)
{
    std::vector<double> values(locations.size(), 0.0);

    // Check the source type
    if(mType == SourceType::POINT)
    {
        if(mPoints.size()==0)
        {
            EXCEPTION("A point is required for this type of source");
        }

        for(unsigned idx=0; idx<locations.size(); idx++)
        {
            for(unsigned jdx=0; jdx<mPoints.size(); jdx++)
            {
                if(norm_2(locations[idx]-mPoints[jdx]) < tolerance)
                {
                    values[idx] += mValue;
                }
            }
        }
    }
    else if(mType == SourceType::VESSEL)
    {
        if(!mpNetwork)
        {
            EXCEPTION("A vessel network is required for this type of source");
        }

        for(unsigned idx=0; idx<locations.size(); idx++)
        {
            std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = this->mpNetwork->GetVesselSegments();
            for (unsigned jdx = 0; jdx < segments.size(); jdx++)
            {
                if (segments[jdx]->GetDistance(locations[idx]) <= tolerance)
                {
                    if(SourceStrength::PRESCRIBED)
                    {
                        values[idx] += mValue;
                    }
                    else
                    {
                        values[idx] += segments[idx]->template GetData<double>(mLabel);
                    }
                }
            }
        }
    }
    else if(mType == SourceType::CELL)
    {
        if(!mpCellPopulation)
        {
            EXCEPTION("A cell population is required for this type of source");
        }

        boost::shared_ptr<AbstractCellProperty> apoptotic_property(new AbstractCellProperty);

        // Loop through all points
        for(unsigned idx=0; idx<locations.size(); idx++)
        {
            // Loop through all cells
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mpCellPopulation->Begin();
                 cell_iter != mpCellPopulation->End(); ++cell_iter)
            {
                if (norm_2(mpCellPopulation->GetLocationOfCellCentre(*cell_iter)-locations[idx])<tolerance)
                {
                    // If a mutation specific consumption rate has been specified
                    if(mMutationSpecificConsumptionRateMap.size()>0)
                    {
                        std::vector<std::pair<AbstractCellProperty, double > >::iterator it;

                        // If the cell is apoptotic
                        if ((*cell_iter)->template HasCellProperty<ApoptoticCellProperty>())
                        {
                            for (it = mMutationSpecificConsumptionRateMap.begin(); it != mMutationSpecificConsumptionRateMap.end(); it++)
                            {
                                if (it->first.IsSame(apoptotic_property))
                                {
                                    values[idx] += it->second;
                                    break;
                                }
                            }
                        }
                        else
                        {
                            for (it = mMutationSpecificConsumptionRateMap.begin(); it != mMutationSpecificConsumptionRateMap.end(); it++)
                            {
                                if ((*cell_iter)->GetMutationState()->IsSame(&(it->first)))
                                {
                                    values[idx] += it->second;
                                    break;
                                }
                            }
                        }

                    }
                    else
                    {
                        values[idx] += mValue;
                    }
                }
            }
        }
    }
    else if(mType == SourceType::SOLUTION)
    {
        if(!mpSolution)
        {
            EXCEPTION("A previous solution is required for this type of source");
        }

        // Sample the field at these locations
        vtkSmartPointer<vtkPolyData> p_polydata = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPoints> p_points = vtkSmartPointer<vtkPoints>::New();
        p_points->SetNumberOfPoints(locations.size());
        for(unsigned idx=0; idx< locations.size(); idx++)
        {
            if(DIM==3)
            {
                p_points->SetPoint(idx, locations[idx][0], locations[idx][1], locations[idx][2]);
            }
            else
            {
                p_points->SetPoint(idx, locations[idx][0], locations[idx][1], 0.0);
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
            values[idx] = p_point_data->GetArray(mLabel.c_str())->GetTuple1(idx);
        }
    }
    return values;
}

template<unsigned DIM>
bool DiscreteSource<DIM>::IsLinearInSolution()
{
    return mIsLinearInSolution;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetCellPopulation(AbstractCellPopulation<DIM>& rCellPopulation)
{
    mpCellPopulation = &rCellPopulation;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetIsLinearInSolution(bool isLinear)
{
    mIsLinearInSolution = isLinear;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetLabelName(const std::string& label)
{
    mLabel = label;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetMutationSpecificConsumptionRateMap(std::vector<std::pair<AbstractCellProperty, double > > mutationSpecificConsumptionRateMap)
{
    mMutationSpecificConsumptionRateMap = mutationSpecificConsumptionRateMap;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetPoints(std::vector<c_vector<double, DIM> > points)
{
    mPoints = points;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetSolution(vtkSmartPointer<vtkImageData>  pSolution)
{
    mpSolution = pSolution;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetSource(SourceStrength::Value boundarySource)
{
    mSourceStrength = boundarySource;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetType(SourceType::Value boundaryType)
{
    mType = boundaryType;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetValue(double value)
{
    mValue = value;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

// Explicit instantiation
template class DiscreteSource<2>;
template class DiscreteSource<3>;
