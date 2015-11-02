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
//#include "SimpleCell.hpp"

#include "DensityMap.hpp"

template<unsigned DIM>
DensityMap<DIM>::DensityMap()
    :   AbstractRegularGridHybridSolver<DIM>()
{

}

template<unsigned DIM>
boost::shared_ptr<DensityMap<DIM> > DensityMap<DIM>::Create()
{
    MAKE_PTR(DensityMap<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
DensityMap<DIM>::~DensityMap()
{

}

template<unsigned DIM>
bool DensityMap<DIM>::IsPointInBox(c_vector<double, DIM> point, c_vector<double, DIM> location, double spacing)
{
    bool point_in_box = false;
    if(point[0] >= location[0] -spacing/2.0 && point[0] <= location [0] + spacing/2.0)
    {
        if(point[1] >= location[1] -spacing/2.0 && point[1] <= location [1] + spacing/2.0)
        {
            if(DIM == 3)
            {
                if(point[2] >= location[2] -spacing/2.0 && point[2] <= location [2] + spacing/2.0)
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

template<unsigned DIM>
double DensityMap<DIM>::LengthOfLineInBox(c_vector<double, DIM> start_point, c_vector<double, DIM> end_point, c_vector<double, DIM> location, double spacing)
{
    if(DIM==2)
    {
        EXCEPTION("Line in box method is currently 3D only");
    }

    // If the line is fully in the box return its length
    bool point1_in_box = IsPointInBox(start_point, location, spacing);
    bool point2_in_box = IsPointInBox(end_point, location, spacing);
    if(point1_in_box && point2_in_box)
    {
        return norm_2(end_point - start_point);
    }
    else
    {
        c_vector<double,2*DIM> bounds;
        bounds[0] = location[0] - spacing/2.0;
        bounds[1] = location[0] + spacing/2.0;
        bounds[2] = location[1] - spacing/2.0;
        bounds[3] = location[1] + spacing/2.0;
        if(DIM==3)
        {
            bounds[4] = location[2] - spacing/2.0;
            bounds[5] = location[2] + spacing/2.0;
        }

        double t1;
        double t2;
        int plane1;
        int plane2;
        c_vector<double,DIM> intercept_1;
        c_vector<double,DIM> intercept_2;

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

template<unsigned DIM>
void DensityMap<DIM>::Solve(bool writeSolution)
{
    unsigned number_of_points = this->mExtents[0] * this->mExtents[1] * this->mExtents[2];
    std::vector<double> vessel_solution(number_of_points, 0.0);
    std::vector<double> cell_solution(number_of_points, 0.0);

    std::vector<c_vector<double, DIM> > cell_locations;
//    if(this->mpCellPopulation)
//    {
//        std::vector<boost::shared_ptr<SimpleCell<DIM> > > cells = this->mpCellPopulation->GetCells();
//        for(unsigned idx = 0; idx < cells.size(); idx++)
//        {
//            cell_locations.push_back(cells[idx]->rGetLocation());
//        }
//    }
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments;

    unsigned grid_index;
    if (this->mpNetwork)// || this->mpCellPopulation)
    {
        if(this->mpNetwork)
        {
            segments = this->mpNetwork->GetVesselSegments();
        }
        for (unsigned i = 0; i < this->mExtents[2]; i++) // Z
        {
            for (unsigned j = 0; j < this->mExtents[1]; j++) // Y
            {
                for (unsigned k = 0; k < this->mExtents[0]; k++) // X
                {
                    grid_index = k + this->mExtents[0] * j + this->mExtents[0] * this->mExtents[1] * i;
                    // If the vessel is in the box add it's length
                    c_vector<double, DIM> location;
                    location[0] = double(k) * this->mGridSize + this->mOrigin[0];
                    location[1] = double(j) * this->mGridSize + this->mOrigin[1];
                    if(DIM==3)
                    {
                        location[2] = double(i) * this->mGridSize + this->mOrigin[2];
                    }

                    if(this->mpNetwork)
                    {
                        for (unsigned idx = 0; idx <  segments.size(); idx++)
                        {
                            vessel_solution[grid_index] += LengthOfLineInBox(segments[idx]->GetNode(0)->GetLocationVector(),
                                                                             segments[idx]->GetNode(1)->GetLocationVector(),
                                                                             location, this->mGridSize);
                        }
                        vessel_solution[grid_index] /= (std::pow(this->mGridSize,3));
                    }
//                    if(this->mpCellPopulation)
//                    {
//                        for (unsigned idx=0; idx<cell_locations.size(); idx++)
//                        {
//                            if(IsPointInBox(cell_locations[idx], location, this->mGridSize))
//                            {
//                                cell_solution[grid_index] += 1.0;
//                            }
//                        }
//                        cell_solution[grid_index] /= (std::pow(this->mGridSize,3));
//                    }
                }
            }
        }
    }
    std::map<std::string, std::vector<double> > data;
    data["CellDensity"] = cell_solution;
    data["VesselDensity"] = vessel_solution;
    this->UpdateSolution(data);

    if (writeSolution)
    {
        this->Write();
    }
}

// Explicit instantiation
template class DensityMap<2> ;
template class DensityMap<3> ;
