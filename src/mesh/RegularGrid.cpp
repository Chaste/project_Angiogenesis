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
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning for now (gcc4.3)
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkProbeFilter.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RegularGrid<ELEMENT_DIM, SPACE_DIM>::RegularGrid() :
        mSpacing(1.0),
        mExtents(std::vector<unsigned>(3, 10)),
        mOrigin(zero_vector<double>(SPACE_DIM)),
        mpNetwork(),
        mpCellPopulation(),
        mPointCellMap(),
        mPointNodeMap(),
        mPointSegmentMap(),
        mpVtkGrid(),
        mVtkGridIsSetUp(false),
        mNeighbourData()
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
void RegularGrid<ELEMENT_DIM, SPACE_DIM>::CalculateNeighbourData()
{
    mNeighbourData = std::vector<std::vector<unsigned> >(GetNumberOfPoints());
    for (unsigned kdx = 0; kdx < mExtents[2]; kdx++)
    {
        for (unsigned jdx = 0; jdx < mExtents[1]; jdx++)
        {
            for (unsigned idx = 0; idx < mExtents[0]; idx++)
            {
                unsigned index = Get1dGridIndex(idx, jdx, kdx);
                if (idx > 0)
                {
                    mNeighbourData[index].push_back(Get1dGridIndex(idx - 1, jdx, kdx));
                }
                if (idx < mExtents[0] - 1)
                {
                    mNeighbourData[index].push_back(Get1dGridIndex(idx + 1, jdx, kdx));
                }
                if (jdx > 0)
                {
                    mNeighbourData[index].push_back(Get1dGridIndex(idx, jdx - 1, kdx));
                }
                if (jdx < mExtents[1] - 1)
                {
                    mNeighbourData[index].push_back(Get1dGridIndex(idx, jdx + 1, kdx));
                }
                if (kdx > 0)
                {
                    mNeighbourData[index].push_back(Get1dGridIndex(idx, jdx, kdx - 1));
                }
                if (kdx < mExtents[2] - 1)
                {
                    mNeighbourData[index].push_back(Get1dGridIndex(idx, jdx, kdx + 1));
                }
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::vector<std::vector<unsigned> >& RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetNeighbourData()
{
    if (mNeighbourData.size() == 0 or mNeighbourData.size() != GetNumberOfPoints())
    {
        CalculateNeighbourData();
    }
    return mNeighbourData;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetNearestGridIndex(c_vector<double, SPACE_DIM> location)
{
    unsigned x_index = round((location[0] - mOrigin[0]) / mSpacing);
    unsigned y_index = round((location[1] - mOrigin[1]) / mSpacing);
    unsigned z_index = 0;
    if (SPACE_DIM == 3)
    {
        z_index = round((location[2] - mOrigin[2]) / mSpacing);
    }
    return Get1dGridIndex(x_index, y_index, z_index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RegularGrid<ELEMENT_DIM, SPACE_DIM>::GenerateFromPart(boost::shared_ptr<Part<SPACE_DIM> > pPart, double gridSize)
{
    mSpacing = gridSize;
    c_vector<double, 2 * SPACE_DIM> spatial_extents = pPart->GetBoundingBox();
    mExtents[0] = unsigned((spatial_extents[1] - spatial_extents[0]) / gridSize) + 1u;
    mExtents[1] = unsigned((spatial_extents[3] - spatial_extents[2]) / gridSize) + 1u;
    if (SPACE_DIM == 3)
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
    if (SPACE_DIM == 3)
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
void RegularGrid<ELEMENT_DIM, SPACE_DIM>::SetPointValues(std::vector<double> pointSolution)
{
    mPointSolution = pointSolution;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::vector<unsigned> > RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetPointPointMap(
        std::vector<c_vector<double, SPACE_DIM> > inputPoints)
{
    double origin_x = mOrigin[0];
    double origin_y = mOrigin[1];
    double origin_z = 0.0;
    if (SPACE_DIM == 3)
    {
        origin_z = mOrigin[2];
    }

    std::vector<std::vector<unsigned> > point_point_map(GetNumberOfPoints());
    for (unsigned idx = 0; idx < inputPoints.size(); idx++)
    {
        unsigned x_index = round((inputPoints[idx][0] - origin_x) / mSpacing);
        unsigned y_index = round((inputPoints[idx][1] - origin_y) / mSpacing);
        unsigned z_index = 0;
        if (SPACE_DIM == 3)
        {
            z_index = round((inputPoints[idx][2] - origin_z) / mSpacing);
        }

        if (x_index <= mExtents[0] && y_index <= mExtents[1] && z_index <= mExtents[2])
        {
            unsigned grid_index = x_index + y_index * mExtents[0] + z_index * mExtents[0] * mExtents[1];
            point_point_map[grid_index].push_back(idx);
        }
    }
    return point_point_map;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> RegularGrid<ELEMENT_DIM, SPACE_DIM>::InterpolateGridValues(
        std::vector<c_vector<double, SPACE_DIM> > locations, std::vector<double> values, bool useVtk)
{
    std::vector<double> sampled_values(locations.size(), 0.0);

    if (!useVtk)
    {
        double origin_x = mOrigin[0];
        double origin_y = mOrigin[1];
        double origin_z = 0.0;
        if (SPACE_DIM == 3)
        {
            origin_z = mOrigin[2];
        }

        // Use a simple version
        for (unsigned idx = 0; idx < locations.size(); idx++)
        {
            // Get the nearest front bottom left point and distances
            unsigned x_index = unsigned((locations[idx][0] - origin_x) / mSpacing);
            double del_x = locations[idx][0] - (double(x_index) * mSpacing + origin_x);
            unsigned y_index = unsigned((locations[idx][1] - origin_y) / mSpacing);
            double del_y = locations[idx][1] - (double(y_index) * mSpacing + origin_y);
            unsigned z_index = 0;
            double del_z = 0.0;
            if (SPACE_DIM == 3)
            {
                z_index = unsigned((locations[idx][2] - origin_z) / mSpacing);
                del_z = locations[idx][2] - (double(z_index) * mSpacing + origin_z);
            }

            if (x_index < mExtents[0] && y_index < mExtents[1] && z_index < mExtents[2])
            {
                double p000 = values[x_index + y_index * mExtents[0] + z_index * mExtents[0] * mExtents[1]];
                double p100 =
                        x_index < (mExtents[0] - 1) ?
                                values[(x_index + 1) + y_index * mExtents[0] + z_index * mExtents[0] * mExtents[1]] :
                                p000;
                double p010 =
                        y_index < (mExtents[1] - 1) ?
                                values[x_index + (y_index + 1) * mExtents[0] + z_index * mExtents[0] * mExtents[1]] :
                                p000;
                double p110 = 0.0;
                if (x_index == mExtents[0] - 1 && y_index == mExtents[1] - 1)
                {
                    p110 = p000;
                }
                else if (x_index == mExtents[0] - 1)
                {
                    p110 = p010;
                }
                else if (y_index == mExtents[1] - 1)
                {
                    p110 = p100;
                }
                else
                {
                    p110 = values[(x_index + 1) + (y_index + 1) * mExtents[0] + z_index * mExtents[0] * mExtents[1]];
                }
                double p001 = 0.0;
                double p101 = 0.0;
                double p111 = 0.0;
                double p011 = 0.0;
                if (z_index < mExtents[2] - 1)
                {
                    p001 = values[x_index + y_index * mExtents[0] + (z_index + 1) * mExtents[0] * mExtents[1]];
                    p101 = x_index < (mExtents[0] - 1) ?
                            values[(x_index + 1) + y_index * mExtents[0] + (z_index + 1) * mExtents[0] * mExtents[1]] :
                            p001;
                    p011 = y_index < (mExtents[1] - 1) ?
                            values[x_index + (y_index + 1) * mExtents[0] + (z_index + 1) * mExtents[0] * mExtents[1]] :
                            p001;
                    if (x_index == mExtents[0] - 1 && y_index == mExtents[1] - 1)
                    {
                        p111 = p001;
                    }
                    else if (x_index == mExtents[0] - 1)
                    {
                        p111 = p011;
                    }
                    else if (y_index == mExtents[1] - 1)
                    {
                        p111 = p101;
                    }
                    else
                    {
                        p111 = values[(x_index + 1) + (y_index + 1) * mExtents[0]
                                + (z_index + 1) * mExtents[0] * mExtents[1]];
                    }
                }
                else
                {
                    p001 = p000;
                    p101 = p100;
                    p111 = p110;
                    p011 = p010;
                }
                double c0 = p000;
                double c1 = p100 - p000;
                double c2 = p010 - p000;
                double c3 = p001 - p000;
                double c4 = p110 - (p010 + p100) + p000;
                double c5 = p011 - p001 - p010 + p000;
                double c6 = p101 - p001 - p100 + p000;
                double c7 = p111 - p011 - p101 - p110 + p100 + p001 + p010 - p000;
                sampled_values[idx] = c0 + c1 * del_x + c2 * del_y + c3 * del_z + c4 * del_x * del_y
                        + c5 * del_y * del_z + c6 * del_z * del_x + c7 * del_x * del_y * del_z;
            }
            else
            {
                EXCEPTION("Sample point is outside grid.");
            }
        }

    }
    else
    {
        if (!mVtkGridIsSetUp)
        {
            SetUpVtkGrid();
        }

        vtkSmartPointer<vtkDoubleArray> pPointData = vtkSmartPointer<vtkDoubleArray>::New();
        pPointData->SetNumberOfComponents(1);
        pPointData->SetNumberOfTuples(GetNumberOfPoints());
        const std::string ny_name = "test";
        pPointData->SetName(ny_name.c_str());
        for (unsigned idx = 0; idx < GetNumberOfPoints(); idx++)
        {
            pPointData->SetValue(idx, values[idx]);
        }
        mpVtkGrid->GetPointData()->AddArray(pPointData);
        mpVtkGrid->Update();

        // Sample the field at these locations
        vtkSmartPointer<vtkPolyData> p_polydata = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPoints> p_points = vtkSmartPointer<vtkPoints>::New();
        p_points->SetNumberOfPoints(locations.size());
        for (unsigned idx = 0; idx < locations.size(); idx++)
        {
            if (SPACE_DIM == 3)
            {
                p_points->SetPoint(idx, locations[idx][0], locations[idx][1], locations[idx][2]);
            }
            else
            {
                p_points->SetPoint(idx, locations[idx][0], locations[idx][1], 0.0);
            }
        }
        p_polydata->SetPoints(p_points);
        p_polydata->Update();

        vtkSmartPointer<vtkProbeFilter> p_probe_filter = vtkSmartPointer<vtkProbeFilter>::New();
        p_probe_filter->SetInput(p_polydata);
        p_probe_filter->SetSource(mpVtkGrid);
        p_probe_filter->Update();
        vtkSmartPointer<vtkPointData> p_point_data = p_probe_filter->GetPolyDataOutput()->GetPointData();

        unsigned num_points = p_point_data->GetArray("test")->GetNumberOfTuples();
        for (unsigned idx = 0; idx < num_points; idx++)
        {
            sampled_values[idx] = p_point_data->GetArray("test")->GetTuple1(idx);
        }
    }
    return sampled_values;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::vector<std::vector<boost::shared_ptr<VascularNode<SPACE_DIM> > > >& RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetPointNodeMap(
        bool update)
{
    if (!update)
    {
        return mPointNodeMap;
    }

    if (!mpNetwork)
    {
        EXCEPTION("A vessel network has not been set. Can not create a point node map.");
    }

    // Loop over all nodes and associate cells with the points
    double origin_x = mOrigin[0];
    double origin_y = mOrigin[1];
    double origin_z = 0.0;
    if (SPACE_DIM == 3)
    {
        origin_z = mOrigin[2];
    }

    mPointNodeMap = std::vector<std::vector<boost::shared_ptr<VascularNode<SPACE_DIM> > > >(GetNumberOfPoints());
    std::vector<boost::shared_ptr<VascularNode<SPACE_DIM> > > nodes = mpNetwork->GetNodes();

    for (unsigned idx = 0; idx < nodes.size(); idx++)
    {
        c_vector<double, SPACE_DIM> location = nodes[idx]->GetLocationVector();
        unsigned x_index = round((location[0] - origin_x) / mSpacing);
        unsigned y_index = round((location[1] - origin_y) / mSpacing);
        unsigned z_index = 0;
        if (SPACE_DIM == 3)
        {
            z_index = round((location[2] - origin_z) / mSpacing);
        }

        if (x_index <= mExtents[0] && y_index <= mExtents[1] && z_index <= mExtents[2])
        {
            unsigned grid_index = x_index + y_index * mExtents[0] + z_index * mExtents[0] * mExtents[1];
            mPointNodeMap[grid_index].push_back(nodes[idx]);
        }
    }
    return mPointNodeMap;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::vector<std::vector<CellPtr> >& RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetPointCellMap(bool update)
{
    if (!update)
    {
        return mPointCellMap;
    }

    if (!mpCellPopulation)
    {
        EXCEPTION("A cell population has not been set. Can not create a cell point map.");
    }

    // Loop over all cells and associate cells with the points
    double origin_x = mOrigin[0];
    double origin_y = mOrigin[1];
    double origin_z = 0.0;
    if (SPACE_DIM == 3)
    {
        origin_z = mOrigin[2];
    }

    mPointCellMap = std::vector<std::vector<CellPtr> >(GetNumberOfPoints());
    for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = mpCellPopulation->Begin();
            cell_iter != mpCellPopulation->End(); ++cell_iter)
    {
        c_vector<double, SPACE_DIM> location = mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
        unsigned x_index = round((location[0] - origin_x) / mSpacing);
        unsigned y_index = round((location[1] - origin_y) / mSpacing);
        unsigned z_index = 0;
        if (SPACE_DIM == 3)
        {
            z_index = round((location[2] - origin_z) / mSpacing);
        }

        if (x_index <= mExtents[0] && y_index <= mExtents[1] && z_index <= mExtents[2])
        {
            unsigned grid_index = x_index + y_index * mExtents[0] + z_index * mExtents[0] * mExtents[1];
            mPointCellMap[grid_index].push_back(*cell_iter);
        }

    }
    return mPointCellMap;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::vector<boost::shared_ptr<CaVesselSegment<SPACE_DIM> > > > RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetPointSegmentMap(
        bool update, bool useVesselSurface)
{
    if (!update)
    {
        return mPointSegmentMap;
    }

    if (!mpNetwork)
    {
        EXCEPTION("A vessel network has not been set. Can not create a vessel point map.");
    }

    // Loop over all points and segments and associate segments with the points
    mPointSegmentMap = std::vector<std::vector<boost::shared_ptr<CaVesselSegment<SPACE_DIM> > > >(GetNumberOfPoints());
    std::vector<boost::shared_ptr<CaVesselSegment<SPACE_DIM> > > segments = mpNetwork->GetVesselSegments();
    for (unsigned jdx = 0; jdx < segments.size(); jdx++)
    {
        for (unsigned idx = 0; idx < GetNumberOfPoints(); idx++)
        {
            if (!useVesselSurface)
            {
                if (segments[jdx]->GetDistance(GetLocationOf1dIndex(idx)) < sqrt(1.0 / 2.0) * mSpacing)
                {
                    mPointSegmentMap[idx].push_back(segments[jdx]);
                }
            }
            else
            {
                if (segments[jdx]->GetDistance(GetLocationOf1dIndex(idx))
                        < segments[jdx]->GetRadius() + sqrt(1.0 / 2.0) * mSpacing)
                {
                    mPointSegmentMap[idx].push_back(segments[jdx]);
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
c_vector<double, SPACE_DIM> RegularGrid<ELEMENT_DIM, SPACE_DIM>::GetLocation(unsigned x_index, unsigned y_index,
                                                                             unsigned z_index)
{
    c_vector<double, SPACE_DIM> location;
    location[0] = double(x_index) * mSpacing + mOrigin[0];
    location[1] = double(y_index) * mSpacing + mOrigin[1];
    if (SPACE_DIM == 3)
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
    for (unsigned idx = 0; idx < GetNumberOfPoints(); idx++)
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
    if (x_index == 0 || x_index == mExtents[0] - 1)
    {
        return true;
    }
    if (y_index == 0 || y_index == mExtents[1] - 1)
    {
        return true;
    }
    if (ELEMENT_DIM == 3)
    {
        if (z_index == 0 || z_index == mExtents[2] - 1)
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
    unsigned mod_z = gridIndex % (mExtents[0] * mExtents[1]);
    unsigned z_index = (gridIndex - mod_z) / (mExtents[0] * mExtents[1]);
    unsigned mod_y = mod_z % mExtents[0];
    unsigned y_index = (mod_z - mod_y) / mExtents[0];
    unsigned x_index = mod_y;

    double loc_x = double(x_index) * mSpacing + mOrigin[0];
    double loc_y = double(y_index) * mSpacing + mOrigin[1];
    double loc_z = 0.0;
    if (SPACE_DIM == 3)
    {
        loc_z = double(z_index) * mSpacing + mOrigin[2];
    }

    c_vector<double, SPACE_DIM> location = GetLocation(x_index, y_index, z_index);
    if (point[0] >= loc_x - mSpacing / 2.0 && point[0] <= loc_x + mSpacing / 2.0)
    {
        if (point[1] >= loc_y - mSpacing / 2.0 && point[1] <= loc_y + mSpacing / 2.0)
        {
            if (SPACE_DIM == 3)
            {
                if (point[2] >= loc_z - mSpacing / 2.0 && point[2] <= loc_z + mSpacing / 2.0)
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
    if (extents.size() < 3)
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
void RegularGrid<ELEMENT_DIM, SPACE_DIM>::SetUpVtkGrid()
{
    // Set up a VTK grid
    mpVtkGrid = vtkSmartPointer<vtkImageData>::New();
    if (SPACE_DIM == 3)
    {
        mpVtkGrid->SetDimensions(mExtents[0], mExtents[1], mExtents[2]);
    }
    else
    {
        mpVtkGrid->SetDimensions(mExtents[0], mExtents[1], 1);
    }
    mpVtkGrid->SetSpacing(mSpacing, mSpacing, mSpacing);

    if (SPACE_DIM == 3)
    {
        mpVtkGrid->SetOrigin(mOrigin[0], mOrigin[1], mOrigin[2]);
    }
    else
    {
        mpVtkGrid->SetOrigin(mOrigin[0], mOrigin[1], 0.0);
    }
    mVtkGridIsSetUp = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RegularGrid<ELEMENT_DIM, SPACE_DIM>::SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<SPACE_DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RegularGrid<ELEMENT_DIM, SPACE_DIM>::Write(boost::shared_ptr<OutputFileHandler> pFileHandler)
{
    if (!mVtkGridIsSetUp)
    {
        SetUpVtkGrid();
    }

    if (mPointSolution.size() == GetNumberOfPoints())
    {
        vtkSmartPointer<vtkDoubleArray> pPointData = vtkSmartPointer<vtkDoubleArray>::New();
        pPointData->SetNumberOfComponents(1);
        pPointData->SetNumberOfTuples(mPointSolution.size());
        pPointData->SetName("Point Values");
        for (unsigned i = 0; i < mPointSolution.size(); i++)
        {
            pPointData->SetValue(i, mPointSolution[i]);
        }
        mpVtkGrid->GetPointData()->AddArray(pPointData);
    }

    vtkSmartPointer<vtkXMLImageDataWriter> pImageDataWriter = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    pImageDataWriter->SetFileName((pFileHandler->GetOutputDirectoryFullPath() + "/grid.vti").c_str());
    pImageDataWriter->SetInput(mpVtkGrid);
    pImageDataWriter->Update();
    pImageDataWriter->Write();
}

// Explicit instantiation
template class RegularGrid<2> ;
template class RegularGrid<3> ;
