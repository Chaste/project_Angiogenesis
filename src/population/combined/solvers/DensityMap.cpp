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
#include <vtkBox.h>
#include "CaVesselSegment.hpp"
#include "ChastePoint.hpp"
#include "SimpleCell.hpp"

#include "DensityMap.hpp"

DensityMap::DensityMap()

    :   AbstractRegularGridHybridSolver()
{

}

boost::shared_ptr<DensityMap> DensityMap::Create()
{
    MAKE_PTR(DensityMap, pSelf);
    return pSelf;
}

DensityMap::~DensityMap()
{

}

bool DensityMap::IsPointInBox(c_vector<double, 3> point, c_vector<double, 3> location, double spacing)
{
    bool point_in_box = false;
    if(point[0] >= location[0] -spacing/2.0 && point[0] <= location [0] + spacing/2.0)
    {
        if(point[1] >= location[1] -spacing/2.0 && point[1] <= location [1] + spacing/2.0)
        {
            if(point[2] >= location[2] -spacing/2.0 && point[2] <= location [2] + spacing/2.0)
            {
                return true;
            }
        }
    }
    return point_in_box;
}

double DensityMap::LengthOfLineInBox(c_vector<double, 3> start_point, c_vector<double, 3> end_point, c_vector<double, 3> location, double spacing)
{
    // If the line is fully in the box return its length
    bool point1_in_box = IsPointInBox(start_point, location, spacing);
    bool point2_in_box = IsPointInBox(end_point, location, spacing);
    if(point1_in_box && point2_in_box)
    {
        return norm_2(end_point - start_point);
    }
    else
    {
        c_vector<double,6> bounds;
        bounds[0] = location[0] - spacing/2.0;
        bounds[1] = location[0] + spacing/2.0;
        bounds[2] = location[1] - spacing/2.0;
        bounds[3] = location[1] + spacing/2.0;
        bounds[4] = location[2] - spacing/2.0;
        bounds[5] = location[2] + spacing/2.0;

        double t1;
        double t2;
        int plane1;
        int plane2;
        c_vector<double,3> intercept_1;
        c_vector<double,3> intercept_2;

        int in_box = vtkBox::IntersectWithLine(&bounds[0], &start_point[0], &end_point[0], t1, t2, &intercept_1[0], &intercept_2[0], plane1, plane2);

        if(point1_in_box)
        {
            return norm_2(intercept_2 - start_point);
        }

        if(point2_in_box)
        {
            return norm_2(intercept_1 - end_point);
        }

        if(in_box)
        {
            return norm_2(intercept_2 - intercept_1);
        }
        else
        {
            return 0.0;
        }
    }
}

void DensityMap::Solve(bool writeSolution)
{
    unsigned number_of_points = mExtents[0] * mExtents[1] * mExtents[2];
    std::vector<double> vessel_solution(number_of_points, 0.0);
    std::vector<double> cell_solution(number_of_points, 0.0);

    std::vector<c_vector<double, 3> > cell_locations;
    if(mpCellPopulation)
    {
        std::vector<boost::shared_ptr<SimpleCell> > cells = mpCellPopulation->GetCells();
        for(unsigned idx = 0; idx < cells.size(); idx++)
        {
            cell_locations.push_back(cells[idx]->rGetLocation());
        }
    }

    std::vector<boost::shared_ptr<CaVesselSegment<3> > > segments;

    unsigned grid_index;
    if (mpNetwork || mpCellPopulation)
    {
        if(mpNetwork)
        {
            segments = mpNetwork->GetVesselSegments();
        }
        for (unsigned i = 0; i < mExtents[2]; i++) // Z
        {
            for (unsigned j = 0; j < mExtents[1]; j++) // Y
            {
                for (unsigned k = 0; k < mExtents[0]; k++) // X
                {
                    grid_index = k + mExtents[0] * j + mExtents[0] * mExtents[1] * i;

                    // If the vessel is in the box add it's length
                    c_vector<double, 3> location;
                    location[0] = double(k) * mGridSize + mOrigin[0];
                    location[1] = double(j) * mGridSize + mOrigin[1];
                    location[2] = double(i) * mGridSize + mOrigin[2];
                    if(mpNetwork)
                    {
                        for (unsigned idx = 0; idx <  segments.size(); idx++)
                        {
                            vessel_solution[grid_index] += LengthOfLineInBox(segments[idx]->GetNode(0)->GetLocationVector(),
                                                                             segments[idx]->GetNode(1)->GetLocationVector(),
                                                                             location,
                                                                             mGridSize);
                        }
                    }
                    if(mpCellPopulation)
                    {
                        for (unsigned idx=0; idx<cell_locations.size(); idx++)
                        {
                            if(IsPointInBox(cell_locations[idx], location, mGridSize))
                            {
                                cell_solution[grid_index] += 1.0;
                            }
                        }
                    }
                }
            }
        }
    }

    std::map<std::string, std::vector<double> > data;
    data["CellDensity"] = cell_solution;
    data["VesselDensity"] = vessel_solution;
    UpdateSolution(data);

    if (writeSolution)
    {
        Write();
    }
}
