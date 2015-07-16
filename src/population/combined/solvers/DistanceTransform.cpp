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

#include "CaVesselSegment.hpp"
#include "ChastePoint.hpp"
#include "SimpleCell.hpp"

#include "DistanceTransform.hpp"

DistanceTransform::DistanceTransform()

    :   AbstractRegularGridHybridSolver()
{

}

boost::shared_ptr<DistanceTransform> DistanceTransform::Create()
{
    MAKE_PTR(DistanceTransform, pSelf);
    return pSelf;
}

DistanceTransform::~DistanceTransform()
{

}

void DistanceTransform::Solve(bool writeSolution)
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

    unsigned grid_index;
    if (mpNetwork || mpCellPopulation)
    {
        std::vector<boost::shared_ptr<CaVesselSegment<3> > > segments;
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

                    // If the grid point is crossed by a vessel segment add the segment's contribution
                    c_vector<double, 3> location;
                    location[0] = double(k) * mGridSize + mOrigin[0];
                    location[1] = double(j) * mGridSize + mOrigin[1];
                    location[2] = double(i) * mGridSize + mOrigin[2];
                    if(mpNetwork)
                    {
                        double min_distance = 1.e6;
                        for (unsigned idx = 0; idx <  segments.size(); idx++)
                        {
                            double seg_dist = segments[idx]->GetDistance(location);
                            if(seg_dist < min_distance)
                            {
                                min_distance = seg_dist;
                            }
                        }
                        vessel_solution[grid_index] = min_distance;
                    }
                    if(mpCellPopulation)
                    {
                        double min_dist = 1.e6;
                        for (unsigned index=0; index<cell_locations.size(); index++)
                        {
                            double dist = norm_2(location - cell_locations[index]);
                            if(dist < min_dist)
                            {
                                min_dist = dist;
                            }
                        }
                        cell_solution[grid_index] = min_dist;
                    }
                }
            }
        }
    }

    std::map<std::string, std::vector<double> > data;
    data["CellDistance"] = cell_solution;
    data["VesselDistance"] = vessel_solution;
    UpdateSolution(data);

    if (writeSolution)
    {
        Write();
    }
}
