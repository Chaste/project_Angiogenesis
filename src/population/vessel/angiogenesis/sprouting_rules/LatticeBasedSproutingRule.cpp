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

#include <stdio.h>
#include "RandomNumberGenerator.hpp"
#include "CaVesselSegment.hpp"
#include "CaVessel.hpp"
#include "LatticeBasedSproutingRule.hpp"

template<unsigned DIM>
LatticeBasedSproutingRule<DIM>::LatticeBasedSproutingRule()
    : AbstractSproutingRule<DIM>(),
      mpGrid(),
      mpSolver(),
      mCellMotility(60.0), // um^2/hour
      mCellChemotacticParameter(2.0e4 * 60.0), // um^2/hour/nM
      mMaxSproutingProbability(0.00025 * 60.0), //hour^-1
      mHalfMaxSproutingProbability(0.5) // nM
{

}

template<unsigned DIM>
LatticeBasedSproutingRule<DIM>::~LatticeBasedSproutingRule()
{

}

template<unsigned DIM>
void LatticeBasedSproutingRule<DIM>::SetGrid(boost::shared_ptr<RegularGrid<DIM> > pGrid)
{
    mpGrid = pGrid;
}

template<unsigned DIM>
void LatticeBasedSproutingRule<DIM>::SetHybridSolver(boost::shared_ptr<AbstractRegularGridHybridSolver<DIM> > pSolver)
{
    mpSolver = pSolver;
}

template<unsigned DIM>
void LatticeBasedSproutingRule<DIM>::SetCellMotilityParameter(double cellMotility)
{
    mCellMotility = cellMotility;
}

template<unsigned DIM>
void LatticeBasedSproutingRule<DIM>::SetCellChemotacticParameter(double cellChemotacticParameter)
{
    mCellChemotacticParameter = cellChemotacticParameter;
}

template<unsigned DIM>
void LatticeBasedSproutingRule<DIM>::SetMaxSproutingProbability(double maxSproutingProbability)
{
    mMaxSproutingProbability = maxSproutingProbability;
}

template<unsigned DIM>
void LatticeBasedSproutingRule<DIM>::SetHalfMaxSproutingProbability(double halfMaxSproutingProbability)
{
    mHalfMaxSproutingProbability = halfMaxSproutingProbability;
}

template <unsigned DIM>
boost::shared_ptr<LatticeBasedSproutingRule<DIM> > LatticeBasedSproutingRule<DIM>::Create()
{
    MAKE_PTR(LatticeBasedSproutingRule<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
std::vector<c_vector<double, DIM> > LatticeBasedSproutingRule<DIM>::GetSproutDirections(const std::vector<boost::shared_ptr<VascularNode<DIM> > >& rNodes)
{
    if(!this->mpGrid)
    {
        EXCEPTION("A regular grid is required for this type of sprouting rule.");
    }

    if(!mpSolver)
    {
        EXCEPTION("A hybrid solver is required for this type of sprouting rule.");
    }

    if(mpSolver->GetGrid() != this->mpGrid)
    {
        EXCEPTION("The vessel grid and hybrid solver grid must be the same for this type of sprouting rule.");
    }

    // Set up the output directions vector
    std::vector<c_vector<double, DIM> > directions(rNodes.size(), zero_vector<double>(DIM));

    // Get the point-node map from the regular grid
    std::vector<std::vector<boost::shared_ptr<VascularNode<DIM> > > > point_node_map = this->mpGrid->GetPointNodeMap();

    // Get the neighbour data from the regular grid
    std::vector<std::vector<unsigned> > neighbour_indices = this->mpGrid->GetNeighbourData();

    // Get the vegf values from the hybrid solver
    std::vector<double> vegf_values = mpSolver->GetPointSolution();

    // Loop over all nodes, if they sprout set sprout the direction
    for(unsigned idx = 0; idx < rNodes.size(); idx++)
    {
        // Only non tip and branch nodes can sprout
        if(rNodes[idx]->GetNumberOfSegments()!=2)
        {
            continue;
        }

        // If the node is too close to an existing branch/sprout don't sprout
        if(this->mVesselEndCutoff > 0.0)
        {
            if(rNodes[idx]->GetVesselSegment(0)->GetVessel()->GetClosestEndNodeDistance(rNodes[idx]->GetLocationVector())< this->mVesselEndCutoff)
            {
                continue;
            }
            if(rNodes[idx]->GetVesselSegment(1)->GetVessel()->GetClosestEndNodeDistance(rNodes[idx]->GetLocationVector())< this->mVesselEndCutoff)
            {
                continue;
            }
        }

        // Get the grid index of the node
        unsigned x_index = round((rNodes[idx]->GetLocationVector()[0] - this->mpGrid->GetOrigin()[0]) / this->mpGrid->GetSpacing());
        unsigned y_index = round((rNodes[idx]->GetLocationVector()[1] - this->mpGrid->GetOrigin()[1]) / this->mpGrid->GetSpacing());
        unsigned z_index = 0;
        if (DIM == 3)
        {
            z_index = round((rNodes[idx]->GetLocationVector()[2] - this->mpGrid->GetOrigin()[2]) / this->mpGrid->GetSpacing());
        }
        unsigned grid_index = this->mpGrid->Get1dGridIndex(x_index, y_index, z_index);

        // Get the probability of moving into each of the neighbour sites
        std::vector<double> probability_of_moving(neighbour_indices[grid_index].size(), 0.0);
        for(unsigned jdx=0; jdx<neighbour_indices[grid_index].size(); jdx++)
        {
            // make sure that tip cell does not try to move into a location already occupied by the vessel that it comes from
            // i.e. that it doesn't loop back around
            c_vector<double, DIM> neighbour_location = this->mpGrid->GetLocationOf1dIndex(neighbour_indices[grid_index][jdx]);
            bool sprout_already_attached_to_vessel_at_location = false;

            for (unsigned seg_index = 0; seg_index < rNodes[idx]->GetNumberOfSegments(); seg_index++)
            {
                if(rNodes[idx]->GetVesselSegment(seg_index)->GetOppositeNode(rNodes[idx])->IsCoincident(ChastePoint<DIM>(neighbour_location)))
                {
                     sprout_already_attached_to_vessel_at_location = true;
                     break;
                }
            }

            //ensure that the new sprout would not try to cross a vessel which is oriented diagonally
            bool vessel_crosses_line_segment = this->mpVesselNetwork->VesselCrossesLineSegment(neighbour_location, rNodes[idx]->GetLocationVector());
            if (!vessel_crosses_line_segment && !sprout_already_attached_to_vessel_at_location)
            {
                // Get the vegf difference
                double vegf = vegf_values[grid_index];
                double VEGF_diff = (vegf_values[neighbour_indices[grid_index][jdx]] - vegf);
                double dt = SimulationTime::Instance()->GetTimeStep();
                double dij = norm_2(rNodes[idx]->GetLocationVector() - neighbour_location);
                probability_of_moving[idx] = ((mCellMotility * dt)/(2.0*dij*dij))*(1.0 + mCellChemotacticParameter*VEGF_diff/(2.0*mCellMotility));
                if (probability_of_moving[jdx] < 0.0)
                {
                    probability_of_moving[jdx] = 0.0;
                }
            }
        }

        if(std::fabs(std::accumulate(probability_of_moving.begin(),probability_of_moving.end(),0.0)) > 1e-16)
        {
            // use roulette-wheel style selection to select which location the tip will move into
            unsigned location_index = neighbour_indices[grid_index].size()*2;
            std::vector<double> cumulativeProbabilityVector(probability_of_moving.size());
            std::partial_sum(probability_of_moving.begin(), probability_of_moving.end(), cumulativeProbabilityVector.begin());
            if (cumulativeProbabilityVector.back() > 1.0)
            {
                std::string message;
                message = "Cumulative probability of tip cell moving is greater than one (" + boost::lexical_cast<std::string>(cumulativeProbabilityVector.back()) + "). Reduce time-step accordingly.";
                EXCEPTION(message);
            }

            double random_number = RandomNumberGenerator::Instance()->ranf();
            for (unsigned ind = 0; ind < cumulativeProbabilityVector.size(); ind++)
            {
                if (random_number <= cumulativeProbabilityVector[ind])
                {
                    location_index = ind;
                    break;
                }
            }
            c_vector<double, DIM> sprout_location = this->mpGrid->GetLocationOf1dIndex(neighbour_indices[grid_index][location_index]);
            directions[idx] = (sprout_location - rNodes[idx]->GetLocationVector()) / norm_2(sprout_location - rNodes[idx]->GetLocationVector());
        }
    }
    return directions;
}

// Explicit instantiation
template class LatticeBasedSproutingRule<2> ;
template class LatticeBasedSproutingRule<3> ;
