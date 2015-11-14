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

#include "Exception.hpp"
#include "RegularGrid.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RegularGrid<ELEMENT_DIM, SPACE_DIM>::RegularGrid() :
    mSpacing(1.0),
    mExtents(std::vector<unsigned>(3,10)),
    mOrigin(zero_vector<double>(SPACE_DIM)),
    mpNetwork(),
    mpCellPopulation(),
    mPointCellMap(),
    mPointSegmentMap()
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
boost::shared_ptr<RegularGrid<ELEMENT_DIM, SPACE_DIM> > RegularGrid<ELEMENT_DIM, SPACE_DIM>::Create()
{
    typedef RegularGrid<ELEMENT_DIM, SPACE_DIM> Reg_Grid_Templated;
    MAKE_PTR(Reg_Grid_Templated, pSelf);
    return pSelf;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RegularGrid<ELEMENT_DIM, SPACE_DIM>::~RegularGrid()
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RegularGrid<ELEMENT_DIM, SPACE_DIM>::GenerateFromPart(boost::shared_ptr<Part<SPACE_DIM> > pPart, double gridSize)
{
    mSpacing = gridSize;
    c_vector<double, 2*SPACE_DIM> spatial_extents = pPart->GetBoundingBox();
    mExtents[0] = unsigned((spatial_extents[1] - spatial_extents[0]) / gridSize) + 1u;
    mExtents[1] = unsigned((spatial_extents[3] - spatial_extents[2]) / gridSize) + 1u;
    if(SPACE_DIM==3)
    {
        mExtents[2] = unsigned((spatial_extents[5] - spatial_extents[4]) / gridSize) + 1u;
    }
    else
    {
        mExtents[2] = 1;
    }

    c_vector<double, SPACE_DIM> origin;
    mOrigin[0] = spatial_extents[0];
    mOrigin[1] = spatial_extents[2];
    if(SPACE_DIM==3)
    {
        mOrigin[2] = spatial_extents[4];
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned RegularGrid<ELEMENT_DIM, SPACE_DIM>::Get1dGridIndex(unsigned x_index, unsigned y_index, unsigned z_index)
{
    return x_index + mExtents[0] * y_index + mExtents[0] * mExtents[1] * z_index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::vector<unsigned> > RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetPointPointMap(std::vector<c_vector<double, SPACE_DIM> > inputPoints)
{
    std::vector<std::vector<unsigned> > point_point_map(GetNumberOfPoints());
    for(unsigned idx=0; idx<GetNumberOfPoints(); idx++)
    {
        std::vector<unsigned> point_index_vector;

        // Loop through all input points
        for(unsigned jdx=0; jdx<inputPoints.size(); jdx++)
        {
            if (IsLocationInPointVolume(inputPoints[jdx], idx))
            {
                point_index_vector.push_back(jdx);
            }
        }
    }
    return point_point_map;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::vector<CellPtr> > RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetPointCellMap(bool update)
{
    if(!update)
    {
        return mPointCellMap;
    }

    if(!mpCellPopulation)
    {
        EXCEPTION("A cell population has not been set. Can not create a cell point map.");
    }

    // Loop over all points and cells and associate cells with the points
    mPointCellMap = std::vector<std::vector<CellPtr> >(GetNumberOfPoints());
    for(unsigned idx=0; idx<GetNumberOfPoints(); idx++)
    {
        std::vector<CellPtr> cell_vector;

        // Loop through all cells
        for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = mpCellPopulation->Begin();
             cell_iter != mpCellPopulation->End(); ++cell_iter)
        {
            if (IsLocationInPointVolume(mpCellPopulation->GetLocationOfCellCentre(*cell_iter), idx))
            {
                cell_vector.push_back(*cell_iter);
            }
        }
    }
    return mPointCellMap;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::vector<boost::shared_ptr<CaVesselSegment<SPACE_DIM> > > > RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetPointSegmentMap(bool update, bool useVesselSurface)
{
    if(!update)
    {
        return mPointSegmentMap;
    }

    if(!mpNetwork)
    {
        EXCEPTION("A vessel network has not been set. Can not create a vessel point map.");
    }

    // Loop over all points and segments and associate segments with the points
    mPointSegmentMap = std::vector<std::vector<boost::shared_ptr<CaVesselSegment<SPACE_DIM> > > >(GetNumberOfPoints());
    for(unsigned idx=0; idx<GetNumberOfPoints(); idx++)
    {
        std::vector<boost::shared_ptr<CaVesselSegment<SPACE_DIM> > > segment_vector;
        std::vector<boost::shared_ptr<CaVesselSegment<SPACE_DIM> > > segments = mpNetwork->GetVesselSegments();
        for (unsigned jdx = 0; jdx < segments.size(); jdx++)
        {
            if(!useVesselSurface)
            {
                if (segments[jdx]->GetDistance(GetLocationOf1dIndex(idx)) <= mSpacing/2.0)
                {
                    segment_vector.push_back(segments[jdx]);
                }
            }
            else
            {
                if(segments[jdx]->GetDistance(GetLocationOf1dIndex(idx)) <= segments[jdx]->GetRadius() +  mSpacing/2.0)
                {
                    segment_vector.push_back(segments[jdx]);
                }
            }
        }
    }
    return mPointSegmentMap;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned> RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetExtents()
{
    return mExtents;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetLocation(unsigned x_index, unsigned y_index, unsigned z_index)
{
    c_vector<double, SPACE_DIM> location;
    location[0] = double(x_index) * mSpacing + mOrigin[0];
    location[1] = double(y_index) * mSpacing + mOrigin[1];
    if(SPACE_DIM == 3)
    {
        location[2] = double(z_index) * mSpacing + mOrigin[2];
    }
    return location;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetLocationOf1dIndex(unsigned grid_index)
{
    unsigned mod_z = grid_index % (mExtents[0] * mExtents[1]);
    unsigned z_index = (grid_index - mod_z) / (mExtents[0] * mExtents[1]);
    unsigned mod_y = mod_z % mExtents[0];
    unsigned y_index = (mod_z - mod_y) / mExtents[0];
    unsigned x_index = mod_y;
    return GetLocation(x_index, y_index, z_index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<c_vector<double, SPACE_DIM> > RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetLocations()
{
    std::vector<c_vector<double, SPACE_DIM> > locations(GetNumberOfPoints());
    for(unsigned idx=0; idx<GetNumberOfPoints(); idx++)
    {
        locations[idx] = GetLocationOf1dIndex(idx);
    }
    return locations;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetNumberOfPoints()
{
    return mExtents[0] * mExtents[1] * mExtents[2];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetOrigin()
{
    return mOrigin;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetSpacing()
{
    return mSpacing;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool RegularGrid<ELEMENT_DIM, SPACE_DIM>::IsOnBoundary(unsigned grid_index)
{
    unsigned mod_z = grid_index % (mExtents[0] * mExtents[1]);
    unsigned z_index = (grid_index - mod_z) / (mExtents[0] * mExtents[1]);
    unsigned mod_y = mod_z % mExtents[0];
    unsigned y_index = (mod_z - mod_y) / mExtents[0];
    unsigned x_index = mod_y;
    return IsOnBoundary(x_index, y_index, z_index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool RegularGrid<ELEMENT_DIM, SPACE_DIM>::IsOnBoundary(unsigned x_index, unsigned y_index, unsigned z_index)
{
    if(x_index == 0 || x_index == mExtents[0] - 1)
    {
        return true;
    }
    if(y_index == 0 || y_index == mExtents[1] - 1)
    {
        return true;
    }
    if(ELEMENT_DIM == 3)
    {
        if(z_index == 0 || z_index == mExtents[2] - 1)
        {
            return true;
        }
    }
    return false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool RegularGrid<ELEMENT_DIM, SPACE_DIM>::IsLocationInPointVolume(c_vector<double, SPACE_DIM> point, unsigned gridIndex)
{
    bool point_in_box = false;
    c_vector<double, SPACE_DIM> location = GetLocationOf1dIndex(gridIndex);
    if(point[0] >= location[0] -mSpacing/2.0 && point[0] <= location [0] + mSpacing/2.0)
    {
        if(point[1] >= location[1] -mSpacing/2.0 && point[1] <= location [1] + mSpacing/2.0)
        {
            if(SPACE_DIM==3)
            {
                if(point[2] >= location[2] -mSpacing/2.0 && point[2] <= location [2] + mSpacing/2.0)
                {
                    return true;
                }
            }
            else
            {
                return true;
            }
        }
    }
    return point_in_box;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RegularGrid<ELEMENT_DIM, SPACE_DIM>::SetCellPopulation(AbstractCellPopulation<SPACE_DIM>& rCellPopulation)
{
    mpCellPopulation = &rCellPopulation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RegularGrid<ELEMENT_DIM, SPACE_DIM>::SetExtents(std::vector<unsigned> extents)
{
    if(extents.size()<3)
    {
        EXCEPTION("The extents should be of dimension 3, regardless of element or space dimension");
    }
    mExtents = extents;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RegularGrid<ELEMENT_DIM, SPACE_DIM>::SetOrigin(c_vector<double, SPACE_DIM> origin)
{
    mOrigin = origin;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RegularGrid<ELEMENT_DIM, SPACE_DIM>::SetSpacing(double spacing)
{
    mSpacing = spacing;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RegularGrid<ELEMENT_DIM, SPACE_DIM>::SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<SPACE_DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

// Explicit instantiation
template class RegularGrid<2> ;
template class RegularGrid<3> ;
