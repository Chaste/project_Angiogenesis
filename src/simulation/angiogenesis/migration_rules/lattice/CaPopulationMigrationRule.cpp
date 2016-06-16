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

#include "GeometryTools.hpp"
#include "CaPopulationMigrationRule.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned DIM>
CaPopulationMigrationRule<DIM>::CaPopulationMigrationRule()
    : AbstractMigrationRule<DIM>(),
      mpGrid(),
      mMovementProbability(0.01),
      mVolumeFractionMap(),
      mpCellPopulation()
{

}

template <unsigned DIM>
boost::shared_ptr<CaPopulationMigrationRule<DIM> > CaPopulationMigrationRule<DIM>::Create()
{
    MAKE_PTR(CaPopulationMigrationRule<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
CaPopulationMigrationRule<DIM>::~CaPopulationMigrationRule()
{

}

template<unsigned DIM>
void CaPopulationMigrationRule<DIM>::SetCellPopulation(boost::shared_ptr<CaBasedCellPopulation<DIM> > cell_population)
{
    mpCellPopulation = cell_population;
}

template<unsigned DIM>
void CaPopulationMigrationRule<DIM>::SetVolumeFraction(boost::shared_ptr<AbstractCellMutationState> mutation_state, double volume_fraction)
{
    if(volume_fraction >1.0)
    {
        EXCEPTION("Specified volume fractions should not be greater than 1.");
    }

    typedef std::map<boost::shared_ptr<AbstractCellMutationState> , double>::iterator it_type; it_type iterator;
    for(iterator = mVolumeFractionMap.begin(); iterator != mVolumeFractionMap.end(); iterator++)
    {
        if (iterator->first->IsSame(mutation_state))
        {
            iterator->second = volume_fraction;
            break;
        }
    }

    // if mutation state does not exist in map yet then add it to the map
    if (iterator == mVolumeFractionMap.end())
    {
        mVolumeFractionMap[mutation_state] = volume_fraction;
    }
}

template<unsigned DIM>
double CaPopulationMigrationRule<DIM>::GetOccupyingVolumeFraction(boost::shared_ptr<AbstractCellMutationState> mutation_state)
{
    typedef std::map<boost::shared_ptr<AbstractCellMutationState> , double>::iterator it_type; it_type iterator;
    for(iterator = mVolumeFractionMap.begin(); iterator != mVolumeFractionMap.end(); iterator++)
    {
        if (iterator->first->IsSame(mutation_state))
        {
            return iterator->second;
        }
    }

    // if a map is not provided or if the prescribed mutation state is not in the map then the
    // occupying volume fraction is 1.
    return 1;
}

template<unsigned DIM>
unsigned CaPopulationMigrationRule<DIM>::GetMaximumCarryingCapacity(boost::shared_ptr<AbstractCellMutationState> mutation_state)
{
    return unsigned(1.0/GetOccupyingVolumeFraction(mutation_state));
}

template<unsigned DIM>
double CaPopulationMigrationRule<DIM>::GetOccupiedVolumeFraction(unsigned index)
{
    std::set<CellPtr> cells = mpCellPopulation->GetCellsUsingLocationIndex(index);
    std::set<CellPtr>::iterator it;

    double occupied_volume_fraction = 0;

    for (it = cells.begin(); it != cells.end(); ++it)
    {
        occupied_volume_fraction += GetOccupyingVolumeFraction((*it)->GetMutationState());
    }

    assert(occupied_volume_fraction <= 1.0);
    return occupied_volume_fraction;
}

template<unsigned DIM>
std::vector<double> CaPopulationMigrationRule<DIM>::GetNeighbourMovementProbabilities(boost::shared_ptr<VascularNode<DIM> > pNode,
                                                       std::vector<unsigned> neighbourIndices, unsigned gridIndex)
{
    std::vector<double> probability_of_moving(neighbourIndices.size(), 0.0);
    for(unsigned idx=0; idx<neighbourIndices.size(); idx++)
    {
        // Make sure that tip cell does not try to move into a location already occupied by the vessel that it comes from
        c_vector<double, DIM> neighbour_location = this->mpGrid->GetLocationOf1dIndex(neighbourIndices[idx]);

        bool already_attached = false;
        for (unsigned seg_index = 0; seg_index < pNode->GetNumberOfSegments(); seg_index++)
        {
            if(pNode->GetVesselSegment(seg_index)->GetOppositeNode(pNode)->IsCoincident(ChastePoint<DIM>(neighbour_location)))
            {
                already_attached = true;
                break;
            }
        }

        // Also ensure that the new location would not try to cross a vessel which is oriented diagonally
        // TODO: Very slow, bottleneck
//        if(already_attached or this->mpVesselNetwork->VesselCrossesLineSegment(neighbour_location, pNode->GetLocationVector()))
//        {
//            continue;
//        }

        if(already_attached)
        {
            continue;
        }

        // Simple rule, equal probability for all directions
        probability_of_moving[idx] = mMovementProbability * SimulationTime::Instance()->GetTimeStep();
    }
    return probability_of_moving;
}

template<unsigned DIM>
int CaPopulationMigrationRule<DIM>::GetNeighbourMovementIndex(std::vector<double> movementProbabilities,
                                                                   std::vector<unsigned> neighbourIndices)
{
    int location_index = -1;

    // Check that the cumulative movement probability is less than one, otherwise our time-step is too large
    std::vector<double> cumulativeProbabilityVector(movementProbabilities.size());
    std::partial_sum(movementProbabilities.begin(), movementProbabilities.end(), cumulativeProbabilityVector.begin());

    if (cumulativeProbabilityVector.back() > 1.0)
    {
        EXCEPTION("Cumulative probability of tip cell moving is greater than one");
    }

    // Use roulette-wheel style selection to select which location the tip will move into
    double cumulativeProbability = cumulativeProbabilityVector.back();
    double random_number = RandomNumberGenerator::Instance()->ranf();

    // If we move, choose a node to go to
    if(random_number < cumulativeProbability)
    {
        for (unsigned ind = 0; ind < cumulativeProbabilityVector.size(); ind++)
        {
            if (random_number <= cumulativeProbabilityVector[ind])
            {
                location_index = neighbourIndices[ind];
                break;
            }
        }
    }
    return location_index;
}

template<unsigned DIM>
void CaPopulationMigrationRule<DIM>::SetGrid(boost::shared_ptr<RegularGrid<DIM> > pGrid)
{
    mpGrid = pGrid;
}

template<unsigned DIM>
void CaPopulationMigrationRule<DIM>::SetMovementProbability(double movementProbability)
{
    mMovementProbability = movementProbability;
}

template<unsigned DIM>
std::vector<int> CaPopulationMigrationRule<DIM>::GetIndices(const std::vector<boost::shared_ptr<VascularNode<DIM> > >& rNodes)
{
    if(!this->mpGrid)
    {
        EXCEPTION("A regular grid is required for this type of migration rule.");
    }

    if(!this->mpVesselNetwork)
    {
        EXCEPTION("A vessel network is required for this type of migration rule.");
    }

    // Set up the output indices vector
    std::vector<int> indices(rNodes.size(), -1);

    // Get the point-node map from the regular grid
    std::vector<std::vector<boost::shared_ptr<VascularNode<DIM> > > > point_node_map = this->mpGrid->GetPointNodeMap();

    // Get the neighbour data from the regular grid
    std::vector<std::vector<unsigned> > neighbour_indices = this->mpGrid->GetNeighbourData();

    // Loop over all nodes, if they can move set the index
    for(unsigned idx = 0; idx < rNodes.size(); idx++)
    {
        // Get the grid index of the node
        unsigned grid_index = this->mpGrid->GetNearestGridIndex(rNodes[idx]->GetLocationVector());

        // Get the probability of moving into each of the neighbour sites
        std::vector<double> probability_of_moving = GetNeighbourMovementProbabilities(rNodes[idx], neighbour_indices[grid_index], grid_index);

        // Get the index of the neighbour to move into
        double sum = std::fabs(std::accumulate(probability_of_moving.begin(), probability_of_moving.end(), 0.0));
        if(sum > 0.0)
        {
            indices[idx] = GetNeighbourMovementIndex(probability_of_moving, neighbour_indices[grid_index]);
        }
    }
    return indices;
}

// Explicit instantiation
template class CaPopulationMigrationRule<2> ;
template class CaPopulationMigrationRule<3> ;
