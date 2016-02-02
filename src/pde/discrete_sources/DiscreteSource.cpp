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

#include "DiscreteSource.hpp"
#include "AbstractCellPopulation.hpp"
#include "CaVascularNetwork.hpp"
#include "GeometryTools.hpp"

template<unsigned DIM>
DiscreteSource<DIM>::DiscreteSource()
    :   mpRegularGrid(),
        mpMesh(),
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
std::vector<double> DiscreteSource<DIM>::GetMeshValues()
{
    if(!mpMesh)
    {
        EXCEPTION("A mesh is required for this type of source");
    }
    return std::vector<double>();

    // TODO Implementations of getpointvalues etc.
}

template<unsigned DIM>
std::vector<double> DiscreteSource<DIM>::GetPointRegularGridValues()
{
    std::vector<double> values(mpRegularGrid->GetNumberOfPoints(), 0.0);

    if(mPoints.size()==0)
    {
        EXCEPTION("A point is required for this type of source");
    }

    // Loop through all points
    std::vector<std::vector<unsigned> > point_point_map = mpRegularGrid->GetPointPointMap(mPoints);

    for(unsigned idx=0; idx<point_point_map.size(); idx++)
    {
        values[idx] += mValue * point_point_map[idx].size();
    }
    return values;
}

template<unsigned DIM>
std::vector<double> DiscreteSource<DIM>::GetVesselRegularGridValues()
{
    std::vector<double> values(mpRegularGrid->GetNumberOfPoints(), 0.0);
    std::vector<std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > > point_segment_map = mpRegularGrid->GetPointSegmentMap();
    for(unsigned idx=0; idx<point_segment_map.size(); idx++)
    {
        for (unsigned jdx = 0; jdx < point_segment_map[idx].size(); jdx++)
        {
            if(mSourceStrength == SourceStrength::PRESCRIBED)
            {
                double point_volume = mpRegularGrid->GetSpacing() * mpRegularGrid->GetSpacing();
                if(DIM == 3)
                {
                    point_volume *= mpRegularGrid->GetSpacing();
                }

                double length_in_box = LengthOfLineInBox<DIM>(point_segment_map[idx][jdx]->GetNode(0)->GetLocationVector(),
                                                                             point_segment_map[idx][jdx]->GetNode(1)->GetLocationVector(),
                                                                             mpRegularGrid->GetLocationOf1dIndex(idx), mpRegularGrid->GetSpacing());
                values[idx] += mValue * length_in_box/ point_volume;
            }
            else
            {
                values[idx] += point_segment_map[idx][jdx]->template GetData<double>(mLabel);
            }
        }
    }
    return values;
}

template<unsigned DIM>
std::vector<double> DiscreteSource<DIM>::GetCellRegularGridValues()
{
    std::vector<double> values(mpRegularGrid->GetNumberOfPoints(), 0.0);
    std::vector<std::vector<CellPtr> > point_cell_map = mpRegularGrid->GetPointCellMap();
    for(unsigned idx=0; idx<point_cell_map.size(); idx++)
    {
        values[idx] += mValue * point_cell_map[idx].size();
    }
    return values;
}

template<unsigned DIM>
std::vector<double> DiscreteSource<DIM>::GetSolutionDependentRegularGridValues()
{
    std::vector<double> values(mpRegularGrid->GetNumberOfPoints(), 0.0);
    if(mpSolution.size() != mpRegularGrid->GetNumberOfPoints())
    {
        EXCEPTION("A solution sampled on the grid is required for this type of source");
    }
    for(unsigned idx=0; idx<mpRegularGrid->GetNumberOfPoints(); idx++)
    {
        values[idx] = mpSolution[idx];
    }
    return values;
}

template<unsigned DIM>
std::vector<double> DiscreteSource<DIM>::GetRegularGridValues()
{
    if(!mpRegularGrid)
    {
        EXCEPTION("A regular grid is required for this type of source");
    }

    // Check the source type
    if(mType == SourceType::POINT)
    {
        return GetPointRegularGridValues();
    }
    else if(mType == SourceType::VESSEL)
    {
        return GetVesselRegularGridValues();
    }
    else if(mType == SourceType::CELL)
    {
        return GetCellRegularGridValues();
    }
    else if(mType == SourceType::SOLUTION)
    {
        return GetSolutionDependentRegularGridValues();
    }
    else
    {
        EXCEPTION("Unknown type requested for discrete source");
    }

}

template<unsigned DIM>
bool DiscreteSource<DIM>::IsLinearInSolution()
{
    return mIsLinearInSolution;
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
void DiscreteSource<DIM>::SetPoints(std::vector<c_vector<double, DIM> > points)
{
    mPoints = points;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetMesh(boost::shared_ptr<HybridMesh<DIM, DIM> > pMesh)
{
    mpMesh = pMesh;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetRegularGrid(boost::shared_ptr<RegularGrid<DIM, DIM> > pRegularGrid)
{
    mpRegularGrid = pRegularGrid;
}

template<unsigned DIM>
void DiscreteSource<DIM>::SetSolution(std::vector<double> solution)
{
    mpSolution = solution;
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

// Explicit instantiation
template class DiscreteSource<2>;
template class DiscreteSource<3>;
