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
//#include "SimpleCell.hpp"

#include "DistanceTransform.hpp"

template<unsigned DIM>
DistanceTransform<DIM>::DistanceTransform()

    :   AbstractRegularGridHybridSolver<DIM>(),
        mpNetwork()
{

}

template<unsigned DIM>
boost::shared_ptr<DistanceTransform<DIM> > DistanceTransform<DIM>::Create()
{
    MAKE_PTR(DistanceTransform, pSelf);
    return pSelf;
}

template<unsigned DIM>
DistanceTransform<DIM>::~DistanceTransform()
{

}

template<unsigned DIM>
void DistanceTransform<DIM>::SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

template<unsigned DIM>
void DistanceTransform<DIM>::Solve(bool writeSolution)
{
    unsigned number_of_points = this->mpRegularGrid->GetNumberOfPoints();
    unsigned extents_x = this->mpRegularGrid->GetExtents()[0];
    unsigned extents_y = this->mpRegularGrid->GetExtents()[1];
    unsigned extents_z = this->mpRegularGrid->GetExtents()[2];

    std::vector<double> vessel_solution(number_of_points, 0.0);

    if (this->mpNetwork)
    {
        std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments;
        if(this->mpNetwork)
        {
            segments = this->mpNetwork->GetVesselSegments();
        }
        for (unsigned i = 0; i < extents_z; i++) // Z
        {
            for (unsigned j = 0; j < extents_y; j++) // Y
            {
                for (unsigned k = 0; k < extents_x; k++) // X
                {
                    unsigned grid_index = this->mpRegularGrid->Get1dGridIndex(k, j, i);
                    c_vector<double, DIM> location = this->mpRegularGrid->GetLocation(k ,j, i);
                    if(this->mpNetwork)
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
                }
            }
        }
    }

    std::map<std::string, std::vector<double> > data;
    data["VesselDistance"] = vessel_solution;
    this->UpdateSolution(data);

    if (writeSolution)
    {
        this->Write();
    }
}

// Explicit instantiation
template class DistanceTransform<2> ;
template class DistanceTransform<3> ;
