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
#include "CellBasedPdeHandler.hpp"
#include "CaBasedCellPopulationWithVessels.hpp"
#include "RandomNumberGenerator.hpp"
#include "ReplicatableVector.hpp"
#include "VascularNode.hpp"
#include <numeric>

#include "Debug.hpp"

template<unsigned DIM>
CaBasedCellPopulationWithVessels<DIM>::CaBasedCellPopulationWithVessels(PottsMesh<DIM>& rMesh,
                                                                        std::vector<CellPtr>& rCells,
                                                                        const std::vector<unsigned> locationIndices,
                                                                        unsigned latticeCarryingCapacity,
                                                                        bool deleteMesh,
                                                                        bool validate)
                                                                        : CaBasedCellPopulation<DIM>(rMesh,
                                                                                                     rCells,
                                                                                                     locationIndices,
                                                                                                     latticeCarryingCapacity,
                                                                                                     deleteMesh,
                                                                                                     validate),
                                                                                                     mpNetwork(),
                                                                                                     mTipCells(),
                                                                                                     mp_tip_mutation_state(new TipCellMutationState),
                                                                                                     mp_stalk_mutation_state(new StalkCellMutationState),
                                                                                                     mp_pde_handler(),
                                                                                                     mCellNodeMap(),
                                                                                                     mVolumeFractionMap()
                                                                                                     {

                                                                                                     }

template<unsigned DIM>
CaBasedCellPopulationWithVessels<DIM>::~CaBasedCellPopulationWithVessels()
{

}

template<unsigned DIM>
void CaBasedCellPopulationWithVessels<DIM>::SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

template<unsigned DIM>
void CaBasedCellPopulationWithVessels<DIM>::DoSprouting(std::vector<boost::shared_ptr<Cell> > activeTips)
{
    double D = 1e-12/(2e-5*2e-5);
    double chi = 2e-5/(1e3*2e-5*2e-5);

    // Loop through active tips and move each one
    for (unsigned tip_index = 0; tip_index < activeTips.size(); tip_index++)
    {
        // Get the vessel network node at the tip
        boost::shared_ptr<VascularNode<DIM> > p_node = mCellNodeMap[activeTips[tip_index]];
        c_vector<double,DIM> tip_location = p_node->GetLocationVector();

        // Do the move
        if(activeTips[tip_index]->GetMutationState()->IsSame(mp_tip_mutation_state) && p_node->GetNumberOfSegments() > 1)
        {
            // Get the VEGF concentration at the tip
            double tip_concentration = activeTips[tip_index]->GetCellData()->GetItem("VEGF");

            // Get the neighbour data
            unsigned potts_index = this->GetLocationIndexUsingCell(activeTips[tip_index]);
            std::map<std::string, std::vector<double> > nbr_data = GetNeighbourData(potts_index, activeTips[tip_index]);

            // check that there is space for the tip cell to move into
            bool space_for_tip_to_move_into = (std::fabs(std::accumulate(nbr_data["Occupancy"].begin(),nbr_data["Occupancy"].end(),0.0)) > 1e-6);

            if (!space_for_tip_to_move_into)
            {
                // sprouting fails ... de-select tip cell
                DeselectTipCell(activeTips[tip_index]);

            }
            else
            {

                // sprout


                std::vector<double> probability_of_moving(nbr_data["Occupancy"].size(),0.0);

                for(unsigned idx=0; idx<nbr_data["Index"].size(); idx++)
                {

                    // make sure that tip cell does not try to move into a location already occupied by the vessel that it
                    // comes from
                    c_vector<double, DIM> neighbour_location = this->rGetMesh().GetNode(unsigned(nbr_data["Index"][idx]))->rGetLocation();
                    bool sprout_already_attached_to_vessel_at_location = false;

                    for (unsigned seg_index = 0; seg_index < p_node->GetNumberOfSegments(); seg_index++)
                    {
                        if(p_node->GetVesselSegment(seg_index)->GetOppositeNode(p_node)->IsCoincident(ChastePoint<DIM>(neighbour_location)))
                        {
                            sprout_already_attached_to_vessel_at_location = true;
                            break;
                        }
                    }

                    //ensure that the new sprout would not try to cross a vessel which is orientated diagonally
                    bool vessel_crosses_line_segment = mpNetwork->VesselCrossesLineSegment(neighbour_location,tip_location);

                    bool space_available = std::fabs(nbr_data["Occupancy"][idx] - 1.0) < 1e-6;

                    if (space_available && !vessel_crosses_line_segment && !sprout_already_attached_to_vessel_at_location)
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


                        //                        probabilityOfMoving[idx] = (SimulationTime::Instance()->GetTimeStep()/(2*dij*dij))*(nbr_data["Occupancy"][idx]/(double)GetMaximumCarryingCapacity(mp_tip_mutation_state))*(D + gamma*VEGF_diff/(2*k));

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

                    // sprout

                    // use roulette-wheel style selection to select which location the tip will move into
                    unsigned location_index = nbr_data["Occupancy"].size()*2;

                    std::vector<double> cumulativeProbabilityVector(probability_of_moving.size());

                    std::partial_sum(probability_of_moving.begin(), probability_of_moving.end(), cumulativeProbabilityVector.begin());

                    assert(cumulativeProbabilityVector.size() == probability_of_moving.size());

                    double cumulativeProbability = cumulativeProbabilityVector.back();

                    double random_number = RandomNumberGenerator::Instance()->ranf();

                    for (unsigned ind = 0; ind < cumulativeProbabilityVector.size(); ind++)
                    {
                        // normalise probabilities so that cumulative probability of moving is 1
                        if (random_number <= cumulativeProbabilityVector[ind]/cumulativeProbability)
                        {
                            location_index = ind;
                            break;
                        }

                    }

                    assert(location_index < nbr_data["Occupancy"].size());

                    // make move into selected location

                    unsigned candidate_location_index = unsigned(nbr_data["Index"][location_index]);

                    c_vector<double, DIM> candidate_location = this->rGetMesh().GetNode(candidate_location_index)->rGetLocation();

                    /*
                     * Note: the 'tip_location' here is a tip cell already located on the vessel, as decided by
                     * a sprouting rule. This actually creates the sprout, which moves to the candidate_location.
                     */
                    boost::shared_ptr<CaVessel<DIM> > p_vessel = mpNetwork->FormSprout(tip_location, candidate_location);
                    boost::shared_ptr<VascularNode<DIM> > p_new_node = p_vessel->GetNodeAtOppositeEnd(p_node);

                    // Check for anastamosis
                    std::set<CellPtr> cells = this->GetCellsUsingLocationIndex(candidate_location_index);
                    std::set<CellPtr>::iterator it;
                    for (it = cells.begin(); it != cells.end(); ++it)
                    {
                        // There is already an EC there
                        if((*it)->GetMutationState()->IsSame(mp_tip_mutation_state) || (*it)->GetMutationState()->IsSame(mp_stalk_mutation_state))
                        {
                            // Deselect self
                            DeselectTipCell(activeTips[tip_index]);

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
                            if((*it)->GetMutationState()->IsSame(mp_tip_mutation_state))
                            {
                                DeselectTipCell((*it));
                            }
                            break;
                        }
                    }

                    // If there has been no anastamosis make a new cell
                    if(activeTips[tip_index]->GetMutationState()->IsSame(mp_tip_mutation_state))
                    {
                        // Create a new cell
                        activeTips[tip_index]->ReadyToDivide();
                        CellPtr p_new_cell = activeTips[tip_index]->Divide();

                        // Add new cell to the cell population
                        this->AddCellUsingLocationIndex(candidate_location_index, p_new_cell); // this doesn't actually add a cell!
                        this->mCells.push_back(p_new_cell); // do it manually here
                        mCellNodeMap[p_new_cell] = p_new_node;
                        assert(norm_2(this->GetLocationOfCellCentre(p_new_cell) - mCellNodeMap[p_new_cell]->GetLocationVector())<1.e-4);
                        DeselectTipCell(activeTips[tip_index]);
                        mTipCells.push_back(p_new_cell);
                    }

                }
                else
                {
                    // sprouting fails ... de-select tip cell
                    DeselectTipCell(activeTips[tip_index]);
                }

            }
        }
    }
}

template<unsigned DIM>
void CaBasedCellPopulationWithVessels<DIM>::DoMigration()
{

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
        if(activeTips[tip_index]->GetMutationState()->IsSame(mp_tip_mutation_state) && p_node->GetNumberOfSegments() == 1)
        {
            // Get the VEGF concentration at the tip
            double tip_concentration = activeTips[tip_index]->GetCellData()->GetItem("VEGF");

            // Get the neighbour data
            unsigned potts_index = this->GetLocationIndexUsingCell(activeTips[tip_index]);
            std::map<std::string, std::vector<double> > nbr_data = GetNeighbourData(potts_index, activeTips[tip_index]);

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
                    c_vector<double, DIM> neighbour_location = this->rGetMesh().GetNode(unsigned(nbr_data["Index"][idx]))->rGetLocation();
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


                        //                        probabilityOfMoving[idx] = (SimulationTime::Instance()->GetTimeStep()/(2*dij*dij))*(nbr_data["Occupancy"][idx]/(double)GetMaximumCarryingCapacity(mp_tip_mutation_state))*(D + gamma*VEGF_diff/(2*k));

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

                    // move

                    // use roulette-wheel style selection to select which location the tip will move into
                    unsigned location_index = nbr_data["Occupancy"].size()*2;

                    std::vector<double> cumulativeProbabilityVector(probability_of_moving.size());

                    std::partial_sum(probability_of_moving.begin(), probability_of_moving.end(), cumulativeProbabilityVector.begin());

                    if (cumulativeProbabilityVector.back() > 1.0)
                    {
                        std::string message;
                        message = "Cumulative probability of tip cell moving is greater than one (" + boost::lexical_cast<std::string>(cumulativeProbabilityVector.back()) + "). Reduce time-step accordingly.";
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
                        // make move into selected location

                        unsigned candidate_location_index = unsigned(nbr_data["Index"][location_index]);

                        c_vector<double, DIM> candidate_location = this->rGetMesh().GetNode(candidate_location_index)->rGetLocation();

                        /*
                         * Note: the 'tip_location' here is a tip cell already located on the vessel, as decided by
                         * a sprouting rule. This actually creates the sprout, which moves to the candidate_location.
                         */
                        boost::shared_ptr<CaVessel<DIM> > p_vessel = mpNetwork->FormSprout(tip_location, candidate_location);
                        boost::shared_ptr<VascularNode<DIM> > p_new_node = p_vessel->GetNodeAtOppositeEnd(p_node);

                        // Check for anastamosis
                        std::set<CellPtr> cells = this->GetCellsUsingLocationIndex(candidate_location_index);
                        std::set<CellPtr>::iterator it;
                        for (it = cells.begin(); it != cells.end(); ++it)
                        {
                            // There is already an EC there
                            if((*it)->GetMutationState()->IsSame(mp_tip_mutation_state) || (*it)->GetMutationState()->IsSame(mp_stalk_mutation_state))
                            {
                                // Deselect self
                                DeselectTipCell(activeTips[tip_index]);

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
                                if((*it)->GetMutationState()->IsSame(mp_tip_mutation_state))
                                {
                                    DeselectTipCell((*it));
                                }
                                break;
                            }
                        }

                        // If there has been no anastamosis make a new cell
                        if(activeTips[tip_index]->GetMutationState()->IsSame(mp_tip_mutation_state))
                        {
                            // Create a new cell
                            activeTips[tip_index]->ReadyToDivide();
                            CellPtr p_new_cell = activeTips[tip_index]->Divide();

                            // Add new cell to the cell population
                            this->AddCellUsingLocationIndex(candidate_location_index, p_new_cell); // this doesn't actually add a cell!
                            this->mCells.push_back(p_new_cell); // do it manually here
                            mCellNodeMap[p_new_cell] = p_new_node;
                            assert(norm_2(this->GetLocationOfCellCentre(p_new_cell) - mCellNodeMap[p_new_cell]->GetLocationVector())<1.e-4);
                            DeselectTipCell(activeTips[tip_index]);
                            mTipCells.push_back(p_new_cell);
                        }
                    }
                }
            }
        }
    }
}

template<unsigned DIM>
bool CaBasedCellPopulationWithVessels<DIM>::IsSiteAvailable(unsigned index, CellPtr pCell)
{
    std::vector<unsigned> available_spaces = this->rGetAvailableSpaces();
    double current_occupied_fraction = GetOccupiedVolumeFraction(index);
    double candidate_fraction = GetOccupyingVolumeFraction(pCell->GetMutationState());
    return(current_occupied_fraction + candidate_fraction <=1.0);
    //
    //
    //    if (pCell->GetMutationState()->IsSame(mp_tip_mutation_state))
    //    {
    //        if (available_spaces[index] != 0)
    //        {
    //            return true;
    //        }
    //        else
    //        {
    //            // if any cell at location index is a stalk or tip cell then the
    //            // tip cell can move in to that location
    //            std::set<CellPtr> cells = this->GetCellsUsingLocationIndex(index);
    //            std::set<CellPtr>::iterator it;
    //
    //            for (it = cells.begin(); it != cells.end(); ++it)
    //            {
    //                if((*it)->GetMutationState()->IsSame(mp_tip_mutation_state) ||
    //                        (*it)->GetMutationState()->IsSame(mp_stalk_mutation_state))
    //                {
    //                    return true;
    //                }
    //            }
    //            return false;
    //        }
    //    }
    //    else
    //    {
    //        return (available_spaces[index] != 0);
    //    }
}

template<unsigned DIM>
void CaBasedCellPopulationWithVessels<DIM>::SetVolumeFraction(boost::shared_ptr<AbstractCellMutationState> mutation_state, double volume_fraction)
{
    if(volume_fraction >1.0)
    {
        EXCEPTION("Specified volume fractions should not be greater than 1.");
    }

    typedef std::map<boost::shared_ptr<AbstractCellMutationState> , double>::iterator it_type;
    it_type iterator;
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
double CaBasedCellPopulationWithVessels<DIM>::GetOccupyingVolumeFraction(boost::shared_ptr<AbstractCellMutationState> mutation_state)
{
    typedef std::map<boost::shared_ptr<AbstractCellMutationState> , double>::iterator it_type;
    it_type iterator;
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
unsigned CaBasedCellPopulationWithVessels<DIM>::GetMaximumCarryingCapacity(boost::shared_ptr<AbstractCellMutationState> mutation_state)
{
    return unsigned(1.0/GetOccupyingVolumeFraction(mutation_state));
}

template<unsigned DIM>
double CaBasedCellPopulationWithVessels<DIM>::GetOccupiedVolumeFraction(unsigned index)
{
    std::set<CellPtr> cells = this->GetCellsUsingLocationIndex(index);
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
void CaBasedCellPopulationWithVessels<DIM>::SelectTipCell(boost::shared_ptr<Cell> pCell)
{
    if (pCell->GetMutationState()->IsSame(mp_stalk_mutation_state))
    {
        pCell->SetMutationState(mp_tip_mutation_state);
        mTipCells.push_back(pCell);
    }
    else
    {
        EXCEPTION("Only stalk cells can be selected to be a tip cell.");
    }
}

template<unsigned DIM>
void CaBasedCellPopulationWithVessels<DIM>::DeselectTipCell(boost::shared_ptr<Cell> pCell)
{
    if (pCell->GetMutationState()->IsSame(mp_tip_mutation_state))
    {
        pCell->SetMutationState(mp_stalk_mutation_state);
        typename std::vector<boost::shared_ptr<Cell> >::iterator it = std::find(mTipCells.begin(), mTipCells.end(), pCell);
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

template<unsigned DIM>
unsigned CaBasedCellPopulationWithVessels<DIM>::GetNumberOfTipCells()
{
    return mTipCells.size();
}

template<unsigned DIM>
std::vector<boost::shared_ptr<Cell> > CaBasedCellPopulationWithVessels<DIM>::GetTipCells()
{
    return mTipCells;
}

template<unsigned DIM>
void CaBasedCellPopulationWithVessels<DIM>::AddPdeHandler(boost::shared_ptr<CellBasedPdeHandler<DIM> > pde_handler)
{
    mp_pde_handler = pde_handler;
}

template<unsigned DIM>
std::map<std::string, std::vector<double> > CaBasedCellPopulationWithVessels<DIM>::GetNeighbourData(unsigned meshIndex, CellPtr pCell)
{
    std::map<std::string, std::vector<double> > data_map;
    data_map["Occupancy"] = std::vector<double> ();
    data_map["VEGF"] = std::vector<double> ();
    data_map["Index"] = std::vector<double> ();

    std::set<unsigned> neighbour_potts_indices = static_cast<PottsMesh<DIM>& >((this->mrMesh)).GetMooreNeighbouringNodeIndices(meshIndex);
    for (std::set<unsigned>::iterator it=neighbour_potts_indices.begin(); it!=neighbour_potts_indices.end(); ++it)
    {
        data_map["Occupancy"].push_back(double(IsSiteAvailable(*it, pCell)));
        data_map["Index"].push_back(double(*it));
        c_vector<double, DIM> neighbour_location = this->rGetMesh().GetNode(*it)->rGetLocation();
        data_map["VEGF"].push_back(mp_pde_handler->GetPdeSolutionAtPoint(neighbour_location, "VEGF"));
    }
    return data_map;
}

template<unsigned DIM>
void CaBasedCellPopulationWithVessels<DIM>::UpdateVascularCellPopulation()
{
    // Update the tip cell collection and cell-node map
    mTipCells = std::vector<boost::shared_ptr<Cell> >();
    mCellNodeMap = std::map<boost::shared_ptr<Cell> , boost::shared_ptr<VascularNode<DIM> > >();

    for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = this->Begin(); cell_iter != this->End(); ++cell_iter)
    {
        if ((*cell_iter)->GetMutationState()->IsSame(mp_stalk_mutation_state) || (*cell_iter)->GetMutationState()->IsSame(mp_tip_mutation_state))
        {
            mCellNodeMap[*cell_iter] = mpNetwork->GetNearestNode(this->GetLocationOfCellCentre((*cell_iter)));

            if(mCellNodeMap[*cell_iter]->GetNumberOfSegments() == 0)
            {
                EXCEPTION("The node corresponding to this cell is not attached to any vessels.");
            }

            if((*cell_iter)->GetMutationState()->IsSame(mp_tip_mutation_state))
            {
                mTipCells.push_back(*cell_iter);
            }
        }
    }

    // Shuffle the tip cell collection
    RandomNumberGenerator::Instance()->Shuffle(mTipCells);

    // Do migration of existing tips
    DoMigration();

    // Do Sprouting - Select candidate tips
    //
    double p_sprout_max = 5e-1;
    double half_max_vegf = 0.65; // units of nano_molar
    double radius_of_exclusion = 2;
    std::vector<boost::shared_ptr<Cell> > candidate_tips = std::vector<boost::shared_ptr<Cell> >();
    if(mp_pde_handler)
    {
        for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = this->Begin(); cell_iter != this->End(); ++cell_iter)
        {
            if ((*cell_iter)->GetMutationState()->IsSame(mp_stalk_mutation_state) )
            {
                double vegf_conc = (*cell_iter)->GetCellData()->GetItem("VEGF");
                double prob_tip_selection = p_sprout_max*SimulationTime::Instance()->GetTimeStep()*vegf_conc/(vegf_conc + half_max_vegf);

                // check to see that there are no other tip cells within the prescribed radius of exclusion
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

                if (RandomNumberGenerator::Instance()->ranf() < prob_tip_selection && mCellNodeMap[*cell_iter]->GetNumberOfSegments() > 1 && distance_to_closest_tip_cell > radius_of_exclusion)
                {
                    SelectTipCell(*cell_iter);
                    candidate_tips.push_back(*cell_iter);
                }
            }
        }
    }

    DoSprouting(candidate_tips);
}



// Explicit instantiation
template class CaBasedCellPopulationWithVessels<2>;
template class CaBasedCellPopulationWithVessels<3>;
