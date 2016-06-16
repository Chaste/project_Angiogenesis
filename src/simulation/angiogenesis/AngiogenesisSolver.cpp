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

#include <boost/lexical_cast.hpp>
#include "UblasIncludes.hpp"
#include "RandomNumberGenerator.hpp"
#include "VascularNode.hpp"
#include "VesselSegment.hpp"
#include "VascularNode.hpp"
#include "AngiogenesisSolver.hpp"
#include "StalkCellMutationState.hpp"
#include "TipCellMutationState.hpp"

#include "Debug.hpp"

// helper method for creating a new tip cell, bypassing the need to query the sub-cellular model,
// which otherwise causes problems
CellPtr CreateNewTipCell(CellPtr parent_cell)
{

    // Create daughter cell with modified cell property collection
    CellPtr p_new_cell(new Cell(parent_cell->GetMutationState(), parent_cell->GetCellCycleModel()->CreateCellCycleModel()));

    // Initialise properties of daughter cell
    p_new_cell->GetCellCycleModel()->InitialiseDaughterCell();

    // Set the daughter cell to inherit the apoptosis time of the parent cell
    p_new_cell->SetApoptosisTime(parent_cell->GetApoptosisTime());

    return p_new_cell;
}

template<unsigned DIM>
AngiogenesisSolver<DIM>::AngiogenesisSolver() :
        mpNetwork(),
        mNodeAnastamosisRadius(5.0),
        mpMigrationRule(),
        mpSproutingRule(),
        mpBoundingDomain(),
        mpFileHandler(),
        mpVesselGrid(),
        mpCellPopulation(),
        mTipCells(),
        mCellNodeMap()
{

}

template<unsigned DIM>
AngiogenesisSolver<DIM>::~AngiogenesisSolver()
{

}

template<unsigned DIM>
boost::shared_ptr<AngiogenesisSolver<DIM> > AngiogenesisSolver<DIM>::Create()
{
    MAKE_PTR(AngiogenesisSolver<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
bool AngiogenesisSolver<DIM>::IsSproutingRuleSet()
{
    return bool(mpSproutingRule);
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::SetAnastamosisRadius(double radius)
{
    mNodeAnastamosisRadius = radius;
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::SetBoundingDomain(boost::shared_ptr<Part<DIM> > pDomain)
{
    mpBoundingDomain = pDomain;
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::SetCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> > cell_population)
{
    mpCellPopulation = cell_population;
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::SetMigrationRule(boost::shared_ptr<AbstractMigrationRule<DIM> > pMigrationRule)
{
    mpMigrationRule = pMigrationRule;
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::SetOutputFileHandler(boost::shared_ptr<OutputFileHandler> pHandler)
{
    mpFileHandler = pHandler;
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::SetSproutingRule(boost::shared_ptr<AbstractSproutingRule<DIM> > pSproutingRule)
{
    mpSproutingRule = pSproutingRule;
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::SetVesselGrid(boost::shared_ptr<RegularGrid<DIM> >pVesselGrid)
{
    mpVesselGrid = pVesselGrid;
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::SetVesselNetwork(boost::shared_ptr<VascularNetwork<DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::DoSprouting()
{
    // Get the candidate sprouts and set them as migrating
    std::vector<boost::shared_ptr<VascularNode<DIM> > > candidate_sprouts = mpSproutingRule->GetSprouts(mpNetwork->GetNodes());

    for(unsigned idx=0; idx<candidate_sprouts.size(); idx++)
    {
        candidate_sprouts[idx]->SetIsMigrating(true);
    }

    UpdateNodalPositions(true);

}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::UpdateNodalPositions(bool sprouting)
{
    // Move any nodes marked as migrating, either new sprouts or tips
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = mpNetwork->GetNodes();
    std::vector<boost::shared_ptr<VascularNode<DIM> > > tips;
    for(unsigned idx = 0; idx < nodes.size(); idx++)
    {
        if(sprouting)
        {
            if(nodes[idx]->IsMigrating() && nodes[idx]->GetNumberOfSegments() == 2)
            {
                tips.push_back(nodes[idx]);
            }
        }
        else
        {
            if(nodes[idx]->IsMigrating() && nodes[idx]->GetNumberOfSegments() == 1)
            {
                tips.push_back(nodes[idx]);
            }
        }
    }

    // Do lattice or off lattice movement
    if (mpVesselGrid)
    {
        std::vector<int> indices = mpMigrationRule->GetIndices(tips);

        for(unsigned idx=0; idx<tips.size();idx++)
        {
            if(indices[idx]>=0)
            {
                if(sprouting)
                {
                    mpNetwork->FormSprout(tips[idx]->GetLocation(), ChastePoint<DIM>(mpVesselGrid->GetLocationOf1dIndex(indices[idx])));
                    tips[idx]->SetIsMigrating(false);
                    mpNetwork->UpdateAll();
                }
                else
                {
                    boost::shared_ptr<VascularNode<DIM> > p_new_node = VascularNode<DIM>::Create(tips[idx]);
                    p_new_node->SetLocation(mpVesselGrid->GetLocationOf1dIndex(indices[idx]));
                    mpNetwork->ExtendVessel(tips[idx]->GetVesselSegment(0)->GetVessel(), tips[idx], p_new_node);
                    tips[idx]->SetIsMigrating(false);
                    p_new_node->SetIsMigrating(true);
                    mpNetwork->UpdateAll();
                }
            }
            else
            {
                if(sprouting && nodes[idx]->GetNumberOfSegments() == 2)
                {
                    tips[idx]->SetIsMigrating(false);
                }
            }
        }
    }
    else
    {

        mpMigrationRule->SetIsSprouting(sprouting);
        std::vector<c_vector<double, DIM> > movement_vectors = mpMigrationRule->GetDirections(tips);

        for(unsigned idx=0; idx<tips.size();idx++)
        {
            if(norm_2(movement_vectors[idx])>0.0)
            {
                bool do_move = true;
                if(mpBoundingDomain)
                {
                    if(!mpBoundingDomain->IsPointInPart(tips[idx]->GetLocationVector() + movement_vectors[idx]))
                    {
                        do_move = false;
                    }
                }

                if(do_move)
                {
                    if(sprouting)
                    {
                        mpNetwork->FormSprout(tips[idx]->GetLocation(), ChastePoint<DIM>(tips[idx]->GetLocationVector() + movement_vectors[idx]));
                        tips[idx]->SetIsMigrating(false);
                        mpNetwork->UpdateAll();
                    }
                    else
                    {
                        boost::shared_ptr<VascularNode<DIM> > p_new_node = VascularNode<DIM>::Create(tips[idx]);
                        p_new_node->SetLocation(tips[idx]->GetLocationVector() + movement_vectors[idx]);
                        mpNetwork->ExtendVessel(tips[idx]->GetVesselSegment(0)->GetVessel(), tips[idx], p_new_node);
                        tips[idx]->SetIsMigrating(false);
                        p_new_node->SetIsMigrating(true);
                    }
                }
            }
            else
            {
                if(sprouting && nodes[idx]->GetNumberOfSegments() == 2)
                {
                    tips[idx]->SetIsMigrating(false);
                }
            }
        }
    }

    mpNetwork->UpdateAll();
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::DoAnastamosis()
{
    mpNetwork->UpdateAll();
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = mpNetwork->GetNodes();

    for(unsigned idx = 0; idx < nodes.size(); idx++)
    {
        // If this is currently a tip
        if(nodes[idx]->IsMigrating() && nodes[idx]->GetNumberOfSegments() == 1)
        {

            if(mpVesselGrid)
            {
                std::vector<std::vector<boost::shared_ptr<VascularNode<DIM> > > > point_node_map = mpVesselGrid->GetPointNodeMap();
                unsigned grid_index = mpVesselGrid->GetNearestGridIndex(nodes[idx]->GetLocationVector());

                if(point_node_map[grid_index].size()>=2)
                {
                    boost::shared_ptr<VascularNode<DIM> > p_merge_node = VascularNode<DIM>::Create(nodes[idx]);

                    if(point_node_map[grid_index][0] == nodes[idx])
                    {
                        p_merge_node = mpNetwork->DivideVessel(point_node_map[grid_index][1]->GetVesselSegment(0)->GetVessel(), nodes[idx]->GetLocation());
                    }
                    else
                    {
                        p_merge_node = mpNetwork->DivideVessel(point_node_map[grid_index][0]->GetVesselSegment(0)->GetVessel(), nodes[idx]->GetLocation());
                    }

                    // Replace the tip node with the merge node
                    p_merge_node->SetIsMigrating(false);
                    if(nodes[idx]->GetVesselSegment(0)->GetNode(0) == nodes[idx])
                    {
                        nodes[idx]->GetVesselSegment(0)->ReplaceNode(0, p_merge_node);
                    }
                    else
                    {
                        nodes[idx]->GetVesselSegment(0)->ReplaceNode(1, p_merge_node);
                    }

                    mpNetwork->UpdateAll();
                }

            }
            else
            {
                // Get the nearest segment and check if it is close enough to the node for a merge
                std::pair<boost::shared_ptr<VesselSegment<DIM> >, double> segment_pair = mpNetwork->GetNearestSegment(nodes[idx], false);

                if(segment_pair.second <= mNodeAnastamosisRadius && nodes[idx]->GetVesselSegment(0)->GetLength() > segment_pair.second)
                {
                    // If there is a non-zero anastamosis radius move the tip onto the segment
                    c_vector<double, DIM> original_location = nodes[idx]->GetLocationVector();
                    if(mNodeAnastamosisRadius > 0.0)
                    {
                        c_vector<double, DIM> divide_location = segment_pair.first->GetPointProjection(original_location, true);
                        nodes[idx]->SetLocation(divide_location);
                    }
                    boost::shared_ptr<VascularNode<DIM> > p_merge_node = mpNetwork->DivideVessel(segment_pair.first->GetVessel(),
                                                                                                 nodes[idx]->GetLocation());

                    // Replace the node at the end of the migrating tip with the merge node
                    if((nodes[idx]->GetVesselSegment(0)->GetNode(0) == p_merge_node) ||
                            (nodes[idx]->GetVesselSegment(0)->GetNode(1) == p_merge_node))
                    {
                        nodes[idx]->SetLocation(original_location);
                    }
                    else
                    {
                        p_merge_node->SetIsMigrating(false);
                        if(nodes[idx]->GetVesselSegment(0)->GetNode(0) == nodes[idx])
                        {
                            nodes[idx]->GetVesselSegment(0)->ReplaceNode(0, p_merge_node);
                        }
                        else
                        {
                            nodes[idx]->GetVesselSegment(0)->ReplaceNode(1, p_merge_node);
                        }
                    }
                    mpNetwork->UpdateAll();
                }
            }
        }
    }
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::Increment()
{
    if(!mpNetwork)
    {
        EXCEPTION("The angiogenesis solver needs an initial vessel network");
    }

    if(mpCellPopulation)
    {
        // Update the tip cell collection and cell-node map
        mTipCells = std::vector<boost::shared_ptr<Cell> >();
        mCellNodeMap = std::map<boost::shared_ptr<Cell> , boost::shared_ptr<VascularNode<DIM> > >();
        MAKE_PTR(StalkCellMutationState, p_EC_state);
        MAKE_PTR(TipCellMutationState, p_EC_Tip_state);

        for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = mpCellPopulation->Begin(); cell_iter != mpCellPopulation->End(); ++cell_iter)
        {
            if ((*cell_iter)->GetMutationState()->IsSame(p_EC_state) || (*cell_iter)->GetMutationState()->IsSame(p_EC_Tip_state))
            {
                mCellNodeMap[*cell_iter] = mpNetwork->GetNearestNode(mpCellPopulation->GetLocationOfCellCentre((*cell_iter)));
                if(mCellNodeMap[*cell_iter]->GetNumberOfSegments() == 0)
                {
                    EXCEPTION("The node corresponding to this cell is not attached to any vessels.");
                }
                if((*cell_iter)->GetMutationState()->IsSame(p_EC_Tip_state))
                {
                    mTipCells.push_back(*cell_iter);
                }
            }
        }

        // Shuffle the tip cell collection
        RandomNumberGenerator::Instance()->Shuffle(mTipCells);

        // Do migration of existing tips
        double D = 1e-12/(2e-5*2e-5);
        double chi = 2e-5/(1e3*2e-5*2e-5);

        // Loop through active tips and move each one
        std::vector<boost::shared_ptr<Cell> > activeTips = mTipCells;

        for (unsigned tip_index = 0; tip_index < activeTips.size(); tip_index++)
        {
            // Get the vessel network node at the tip
            boost::shared_ptr<VascularNode<DIM> > p_node = mCellNodeMap[activeTips[tip_index]];
            c_vector<double,DIM> tip_location = p_node->GetLocationVector();

            // Make sure it is the correct type of cell and Do the move
            if(activeTips[tip_index]->GetMutationState()->IsSame(p_EC_Tip_state) && p_node->GetNumberOfSegments() == 1)
            {
                // Get the VEGF concentration at the tip
//                double tip_concentration = activeTips[tip_index]->GetCellData()->GetItem("VEGF");
                double tip_concentration = 0.0;

                // Get the neighbour data
                unsigned potts_index = mpCellPopulation->GetLocationIndexUsingCell(activeTips[tip_index]);

                // Get NBR Data
                std::map<std::string, std::vector<double> > nbr_data;
                nbr_data["Occupancy"] = std::vector<double> ();
                nbr_data["VEGF"] = std::vector<double> ();
                nbr_data["Index"] = std::vector<double> ();
                std::set<unsigned> neighbour_potts_indices = static_cast<PottsMesh<DIM>& >((mpCellPopulation->rGetMesh())).GetMooreNeighbouringNodeIndices(potts_index);
                for (std::set<unsigned>::iterator it=neighbour_potts_indices.begin(); it!=neighbour_potts_indices.end(); ++it)
                {
                    //nbr_data["Occupancy"].push_back(double(mpCellPopulation->IsSiteAvailable(*it, activeTips[tip_index])));
                    nbr_data["Occupancy"].push_back(0.0);
                    //double current_occupied_fraction = GetOccupiedVolumeFraction(index);
                    //double candidate_fraction = GetOccupyingVolumeFraction(pCell->GetMutationState());
                    //return(current_occupied_fraction + candidate_fraction <=1.0);

                    nbr_data["Index"].push_back(double(*it));
                    c_vector<double, DIM> neighbour_location = mpCellPopulation->rGetMesh().GetNode(*it)->rGetLocation();
                    nbr_data["VEGF"].push_back(0.0);
                }

                // check that there is space for the tip cell to move into
                bool space_for_tip_to_move_into = (std::fabs(std::accumulate(nbr_data["Occupancy"].begin(),nbr_data["Occupancy"].end(),0.0)) > 1e-6);

                if (space_for_tip_to_move_into)
                {
                    // try to move
                    std::vector<double> probability_of_moving(nbr_data["Occupancy"].size(),0.0);
                    for(unsigned idx=0; idx<nbr_data["Index"].size(); idx++)
                    {
                        // make sure that tip cell does not try to move into a location already occupied by the vessel that it
                        // comes from
                        c_vector<double, DIM> neighbour_location = mpCellPopulation->rGetMesh().GetNode(unsigned(nbr_data["Index"][idx]))->rGetLocation();
                        bool go_back_on_self = false;

                        for (unsigned seg_index = 0; seg_index < p_node->GetNumberOfSegments(); seg_index++)
                        {
                            if(p_node->GetVesselSegment(seg_index)->GetOppositeNode(p_node)->IsCoincident(ChastePoint<DIM>(neighbour_location)))
                            {
                                go_back_on_self = true;
                                break;
                            }
                        }

                        //ensure that the new sprout would not try to cross a vessel which is orientated diagonally
                        bool vessel_crosses_line_segment = mpNetwork->VesselCrossesLineSegment(neighbour_location,tip_location);
                        bool space_available = std::fabs(nbr_data["Occupancy"][idx] - 1.0) < 1e-6;
                        if (space_available && !vessel_crosses_line_segment && !go_back_on_self)
                        {
                            double k = 1;
                            double dij = norm_2(tip_location - neighbour_location);
                            double VEGF_diff = (nbr_data["VEGF"][idx] - tip_concentration);
                            // need to correct for eventuality that we are using a Moore neighbourhood
                            if (DIM == 2)
                            {
                                k = 2;
                            }
                            if (DIM == 3)
                            {
                                k = 26.0/6.0;
                            }

                            // probabilityOfMoving[idx] = (SimulationTime::Instance()->GetTimeStep()/(2*dij*dij))*(nbr_data["Occupancy"][idx]/(double)GetMaximumCarryingCapacity(mp_tip_mutation_state))*(D + gamma*VEGF_diff/(2*k));
                            probability_of_moving[idx] = (SimulationTime::Instance()->GetTimeStep()/(2*dij*dij))*(D + chi*VEGF_diff/(2*k));
                            if (probability_of_moving[idx] < 0.0)
                            {
                                probability_of_moving[idx] = 0.0;
                            }
                        }
                        else
                        {
                            probability_of_moving[idx] = 0.0;
                        }
                    }

                    if (std::fabs(std::accumulate(probability_of_moving.begin(),probability_of_moving.end(),0.0)) > 1e-16)
                    {
                        // use roulette-wheel style selection to select which location the tip will move into
                        unsigned location_index = nbr_data["Occupancy"].size()*2;
                        std::vector<double> cumulativeProbabilityVector(probability_of_moving.size());
                        std::partial_sum(probability_of_moving.begin(), probability_of_moving.end(), cumulativeProbabilityVector.begin());
                        if (cumulativeProbabilityVector.back() > 1.0)
                        {
                            std::string message;
                            message = "Cumulative probability of tip cell moving is greater than one.";
                            EXCEPTION(message);
                        }
                        assert(cumulativeProbabilityVector.size() == probability_of_moving.size());
                        double random_number = RandomNumberGenerator::Instance()->ranf();
                        for (unsigned ind = 0; ind < cumulativeProbabilityVector.size(); ind++)
                        {
                            if (random_number <= cumulativeProbabilityVector[ind])
                            {
                                location_index = ind;
                                break;
                            }
                        }
                        if (location_index < nbr_data["Occupancy"].size())
                        {
                            unsigned candidate_location_index = unsigned(nbr_data["Index"][location_index]);
                            c_vector<double, DIM> candidate_location = mpCellPopulation->rGetMesh().GetNode(candidate_location_index)->rGetLocation();
                            /*
                             * Note: the 'tip_location' here is a tip cell already located on the vessel, as decided by
                             * a sprouting rule. This actually creates the sprout, which moves to the candidate_location.
                             */
                            boost::shared_ptr<Vessel<DIM> > p_vessel = mpNetwork->FormSprout(tip_location, candidate_location);
                            boost::shared_ptr<VascularNode<DIM> > p_new_node = p_vessel->GetNodeAtOppositeEnd(p_node);

                            // Check for anastamosis
                            std::set<CellPtr> cells = mpCellPopulation->GetCellsUsingLocationIndex(candidate_location_index);
                            std::set<CellPtr>::iterator it;
                            for (it = cells.begin(); it != cells.end(); ++it)
                            {
                                // There is already an EC there
                                if((*it)->GetMutationState()->IsSame(p_EC_Tip_state) || (*it)->GetMutationState()->IsSame(p_EC_state))
                                {
                                    // Deselect self
                                    if (activeTips[tip_index]->GetMutationState()->IsSame(p_EC_Tip_state))
                                    {
                                        activeTips[tip_index]->SetMutationState(p_EC_state);
                                        typename std::vector<boost::shared_ptr<Cell> >::iterator it = std::find(mTipCells.begin(), mTipCells.end(), activeTips[tip_index]);
                                        if(it != mTipCells.end())
                                        {
                                            mTipCells.erase(it);
                                        }
                                        else
                                        {
                                            EXCEPTION("Tip cell not found.");
                                        }
                                    }
                                    else
                                    {
                                        EXCEPTION("Cell is not a tip cell.");
                                    }

                                    boost::shared_ptr<VascularNode<DIM> > p_other_node = mCellNodeMap[(*it)];
                                    mpNetwork->DivideVessel(p_other_node->GetVesselSegment(0)->GetVessel(), candidate_location);

                                    // Replace the tip node with the other node
                                    if(p_new_node->GetVesselSegment(0)->GetNode(0) == p_new_node)
                                    {
                                        p_new_node->GetVesselSegment(0)->ReplaceNode(0, p_other_node);
                                    }
                                    else
                                    {
                                        p_new_node->GetVesselSegment(0)->ReplaceNode(1, p_other_node);
                                    }

                                    mpNetwork->UpdateAll();

                                    // If we are a tip also de-select at neighbour location
                                    if((*it)->GetMutationState()->IsSame(p_EC_Tip_state))
                                    {
                                        if ((*it)->GetMutationState()->IsSame(p_EC_Tip_state))
                                        {
                                            (*it)->SetMutationState(p_EC_state);
                                            typename std::vector<boost::shared_ptr<Cell> >::iterator it = std::find(mTipCells.begin(), mTipCells.end(), (*it));
                                            if(it != mTipCells.end())
                                            {
                                                mTipCells.erase(it);
                                            }
                                            else
                                            {
                                                EXCEPTION("Tip cell not found.");
                                            }
                                        }
                                        else
                                        {
                                            EXCEPTION("Cell is not a tip cell.");
                                        }
                                    }
                                    break;
                                }
                            }

                            // If there has been no anastamosis make a new cell
                            if(activeTips[tip_index]->GetMutationState()->IsSame(p_EC_Tip_state))
                            {
                                // Create a new cell
                                CellPtr p_new_cell = CreateNewTipCell(activeTips[tip_index]);

                                // Add new cell to the cell population
                                mpCellPopulation->AddCellUsingLocationIndex(candidate_location_index, p_new_cell); // this doesn't actually add a cell!
                                //this->mCells.push_back(p_new_cell); // do it manually here...uh oh
                                mCellNodeMap[p_new_cell] = p_new_node;
                                assert(norm_2(mpCellPopulation->GetLocationOfCellCentre(p_new_cell) - mCellNodeMap[p_new_cell]->GetLocationVector())<1.e-4);
                                if (activeTips[tip_index]->GetMutationState()->IsSame(p_EC_Tip_state))
                                {
                                    activeTips[tip_index]->SetMutationState(p_EC_state);
                                    typename std::vector<boost::shared_ptr<Cell> >::iterator it = std::find(mTipCells.begin(), mTipCells.end(), activeTips[tip_index]);
                                    if(it != mTipCells.end())
                                    {
                                        mTipCells.erase(it);
                                    }
                                    else
                                    {
                                        EXCEPTION("Tip cell not found.");
                                    }
                                }
                                else
                                {
                                    EXCEPTION("Cell is not a tip cell.");
                                }
                                mTipCells.push_back(p_new_cell);
                            }
                        }
                    }
                }
            }
        }

        // Do Sprouting - Select candidate tips
        if(mpSproutingRule)
        {
            double p_sprout_max = 5e-1;
            //double half_max_vegf = 0.65e-3; // units of nano_molar
            double radius_of_exclusion = 2;
            std::vector<boost::shared_ptr<Cell> > candidate_tips = std::vector<boost::shared_ptr<Cell> >();

            for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = mpCellPopulation->Begin(); cell_iter != mpCellPopulation->End(); ++cell_iter)
            {
                if ((*cell_iter)->GetMutationState()->IsSame(p_EC_state) )
                {
                    unsigned cell_index = mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
                    c_vector<double, DIM> cell_location = mpCellPopulation->rGetMesh().GetNode(cell_index)->rGetLocation();

                    double prob_tip_selection = p_sprout_max*SimulationTime::Instance()->GetTimeStep()*1.0;
                    double distance_to_closest_tip_cell = radius_of_exclusion*2;
                    if (candidate_tips.size() == 0 && mTipCells.size() == 0)
                    {
                        distance_to_closest_tip_cell = radius_of_exclusion*2;
                    }
                    else
                    {
                        std::vector<boost::shared_ptr<Cell> >::iterator it;
                        for (it = mTipCells.begin(); it != mTipCells.end(); ++it)
                        {
                            if(mCellNodeMap[*cell_iter]->GetDistance(mCellNodeMap[*it]) < radius_of_exclusion)
                            {
                                distance_to_closest_tip_cell = mCellNodeMap[*cell_iter]->GetDistance(mCellNodeMap[*it]);
                                break;
                            }
                        }
                    }

                    if (RandomNumberGenerator::Instance()->ranf() < prob_tip_selection && mCellNodeMap[(*cell_iter)]->GetNumberOfSegments() > 1 && distance_to_closest_tip_cell > radius_of_exclusion)
                    {
                        if ((*cell_iter)->GetMutationState()->IsSame(p_EC_state))
                        {
                            (*cell_iter)->SetMutationState(p_EC_Tip_state);
                            mTipCells.push_back(*cell_iter);
                        }
                        else
                        {
                            EXCEPTION("Only stalk cells can be selected to be a tip cell.");
                        }
                        candidate_tips.push_back(*cell_iter);
                    }
                }
            }

            //DoSprouting(candidate_tips);
//            double D = 1e-12/(2e-5*2e-5);
//            double chi = 2e-5/(1e3*2e-5*2e-5);
//
//            // Loop through active tips and move each one
//            for (unsigned tip_index = 0; tip_index < activeTips.size(); tip_index++)
//            {
//                // Get the vessel network node at the tip
//                boost::shared_ptr<VascularNode<DIM> > p_node = mCellNodeMap[activeTips[tip_index]];
//                c_vector<double,DIM> tip_location = p_node->GetLocationVector();
//
//                // Do the move
//                if(activeTips[tip_index]->GetMutationState()->IsSame(mp_tip_mutation_state) && p_node->GetNumberOfSegments() > 1)
//                {
//                    // Get the VEGF concentration at the tip
//                    double tip_concentration = activeTips[tip_index]->GetCellData()->GetItem("VEGF");
//
//                    // Get the neighbour data
//                    unsigned potts_index = this->GetLocationIndexUsingCell(activeTips[tip_index]);
//                    std::map<std::string, std::vector<double> > nbr_data = GetNeighbourData(potts_index, activeTips[tip_index]);
//
//                    // check that there is space for the tip cell to move into
//                    bool space_for_tip_to_move_into = (std::fabs(std::accumulate(nbr_data["Occupancy"].begin(),nbr_data["Occupancy"].end(),0.0)) > 1e-6);
//
//                    if (!space_for_tip_to_move_into)
//                    {
//                        // sprouting fails ... de-select tip cell
//                        DeselectTipCell(activeTips[tip_index]);
//
//                    }
//                    else
//                    {
//
//                        // sprout
//
//
//                        std::vector<double> probability_of_moving(nbr_data["Occupancy"].size(),0.0);
//
//                        for(unsigned idx=0; idx<nbr_data["Index"].size(); idx++)
//                        {
//
//                            // make sure that tip cell does not try to move into a location already occupied by the vessel that it
//                            // comes from
//                            c_vector<double, DIM> neighbour_location = this->rGetMesh().GetNode(unsigned(nbr_data["Index"][idx]))->rGetLocation();
//                            bool sprout_already_attached_to_vessel_at_location = false;
//
//                            for (unsigned seg_index = 0; seg_index < p_node->GetNumberOfSegments(); seg_index++)
//                            {
//                                if(p_node->GetVesselSegment(seg_index)->GetOppositeNode(p_node)->IsCoincident(ChastePoint<DIM>(neighbour_location)))
//                                {
//                                    sprout_already_attached_to_vessel_at_location = true;
//                                    break;
//                                }
//                            }
//
//                            //ensure that the new sprout would not try to cross a vessel which is orientated diagonally
//                            bool vessel_crosses_line_segment = mpNetwork->VesselCrossesLineSegment(neighbour_location,tip_location);
//
//                            bool space_available = std::fabs(nbr_data["Occupancy"][idx] - 1.0) < 1e-6;
//
//                            if (space_available && !vessel_crosses_line_segment && !sprout_already_attached_to_vessel_at_location)
//                            {
//
//
//                                double k = 1;
//                                double dij = norm_2(tip_location - neighbour_location);
//                                double VEGF_diff = (nbr_data["VEGF"][idx] - tip_concentration);
//
//                                // need to correct for eventuality that we are using a Moore neighbourhood
//                                if (DIM == 2)
//                                {
//                                    k = 2;
//                                }
//                                if (DIM == 3)
//                                {
//                                    k = 26.0/6.0;
//                                }
//
//
//                                //                        probabilityOfMoving[idx] = (SimulationTime::Instance()->GetTimeStep()/(2*dij*dij))*(nbr_data["Occupancy"][idx]/(double)GetMaximumCarryingCapacity(mp_tip_mutation_state))*(D + gamma*VEGF_diff/(2*k));
//
//                                probability_of_moving[idx] = (SimulationTime::Instance()->GetTimeStep()/(2*dij*dij))*(D + chi*VEGF_diff/(2*k));
//
//                                if (probability_of_moving[idx] < 0.0)
//                                {
//                                    probability_of_moving[idx] = 0.0;
//                                }
//
//                            }
//                            else
//                            {
//                                probability_of_moving[idx] = 0.0;
//                            }
//
//                        }
//
//
//                        if (std::fabs(std::accumulate(probability_of_moving.begin(),probability_of_moving.end(),0.0)) > 1e-16)
//                        {
//
//                            // sprout
//
//                            // use roulette-wheel style selection to select which location the tip will move into
//                            unsigned location_index = nbr_data["Occupancy"].size()*2;
//
//                            std::vector<double> cumulativeProbabilityVector(probability_of_moving.size());
//
//                            std::partial_sum(probability_of_moving.begin(), probability_of_moving.end(), cumulativeProbabilityVector.begin());
//
//                            assert(cumulativeProbabilityVector.size() == probability_of_moving.size());
//
//                            double cumulativeProbability = cumulativeProbabilityVector.back();
//
//                            double random_number = RandomNumberGenerator::Instance()->ranf();
//
//                            for (unsigned ind = 0; ind < cumulativeProbabilityVector.size(); ind++)
//                            {
//                                // normalise probabilities so that cumulative probability of moving is 1
//                                if (random_number <= cumulativeProbabilityVector[ind]/cumulativeProbability)
//                                {
//                                    location_index = ind;
//                                    break;
//                                }
//
//                            }
//
//                            assert(location_index < nbr_data["Occupancy"].size());
//
//                            // make move into selected location
//
//                            unsigned candidate_location_index = unsigned(nbr_data["Index"][location_index]);
//
//                            c_vector<double, DIM> candidate_location = this->rGetMesh().GetNode(candidate_location_index)->rGetLocation();
//
//                            /*
//                             * Note: the 'tip_location' here is a tip cell already located on the vessel, as decided by
//                             * a sprouting rule. This actually creates the sprout, which moves to the candidate_location.
//                             */
//                            boost::shared_ptr<Vessel<DIM> > p_vessel = mpNetwork->FormSprout(tip_location, candidate_location);
//                            boost::shared_ptr<VascularNode<DIM> > p_new_node = p_vessel->GetNodeAtOppositeEnd(p_node);
//
//                            // Check for anastamosis
//                            std::set<CellPtr> cells = this->GetCellsUsingLocationIndex(candidate_location_index);
//                            std::set<CellPtr>::iterator it;
//                            for (it = cells.begin(); it != cells.end(); ++it)
//                            {
//                                // There is already an EC there
//                                if((*it)->GetMutationState()->IsSame(mp_tip_mutation_state) || (*it)->GetMutationState()->IsSame(mp_stalk_mutation_state))
//                                {
//                                    // Deselect self
//                                    DeselectTipCell(activeTips[tip_index]);
//
//                                    boost::shared_ptr<VascularNode<DIM> > p_other_node = mCellNodeMap[(*it)];
//                                    mpNetwork->DivideVessel(p_other_node->GetVesselSegment(0)->GetVessel(), candidate_location);
//
//                                    // Replace the tip node with the other node
//                                    if(p_new_node->GetVesselSegment(0)->GetNode(0) == p_new_node)
//                                    {
//                                        p_new_node->GetVesselSegment(0)->ReplaceNode(0, p_other_node);
//                                    }
//                                    else
//                                    {
//                                        p_new_node->GetVesselSegment(0)->ReplaceNode(1, p_other_node);
//                                    }
//
//                                    mpNetwork->UpdateAll();
//
//                                    // If we are a tip also de-select at neighbour location
//                                    if((*it)->GetMutationState()->IsSame(mp_tip_mutation_state))
//                                    {
//                                        DeselectTipCell((*it));
//                                    }
//                                    break;
//                                }
//                            }
//
//                            // If there has been no anastamosis make a new cell
//                            if(activeTips[tip_index]->GetMutationState()->IsSame(mp_tip_mutation_state))
//                            {
//                                // Create a new cell
//                                CellPtr p_new_cell = CreateNewTipCell(activeTips[tip_index]);
//
//                                // Add new cell to the cell population
//                                this->AddCellUsingLocationIndex(candidate_location_index, p_new_cell); // this doesn't actually add a cell!
//                                this->mCells.push_back(p_new_cell); // do it manually here
//                                mCellNodeMap[p_new_cell] = p_new_node;
//                                assert(norm_2(this->GetLocationOfCellCentre(p_new_cell) - mCellNodeMap[p_new_cell]->GetLocationVector())<1.e-4);
//                                DeselectTipCell(activeTips[tip_index]);
//                                mTipCells.push_back(p_new_cell);
//                            }
//
//                        }
//                        else
//                        {
//                            // sprouting fails ... de-select tip cell
//                            DeselectTipCell(activeTips[tip_index]);
//                        }
//
//                    }
//                }
//            }
        }

    }
    else
    {
        // If doing lattice based, add the network to the lattice
        if(mpVesselGrid)
        {
            mpVesselGrid->SetVesselNetwork(mpNetwork);
        }

        // Move any migrating nodes
        UpdateNodalPositions();

        // Check for anastamosis
        DoAnastamosis();

        // Do sprouting
        if(mpSproutingRule)
        {
            mpSproutingRule->SetVesselNetwork(mpNetwork);

            if(mpVesselGrid)
            {
                mpSproutingRule->SetGrid(mpVesselGrid);
            }

            DoSprouting();
            DoAnastamosis();
        }
    }
}

template<unsigned DIM>
void AngiogenesisSolver<DIM>::Run(bool writeOutput)
{
    // Loop for the duration of the simulation time
    while(!SimulationTime::Instance()->IsFinished())
    {
        // Write the vessel network if appropriate
        if(writeOutput && mpFileHandler && mpNetwork)
        {
            mpNetwork->Write(mpFileHandler->GetOutputDirectoryFullPath() + "/vessel_network_" +
                             boost::lexical_cast<std::string>(SimulationTime::Instance()->GetTimeStepsElapsed()) + ".vtp");
            if(mpCellPopulation)
            {
                mpCellPopulation->OpenWritersFiles(*mpFileHandler);
                mpCellPopulation->WriteResultsToFiles(mpFileHandler->GetRelativePath());
                mpCellPopulation->CloseWritersFiles();
            }
        }

        // Increment the solver and simulation time
        Increment();
        SimulationTime::Instance()->IncrementTimeOneStep();
    }
}

// Explicit instantiation
template class AngiogenesisSolver<2> ;
template class AngiogenesisSolver<3> ;
