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
    // Loop through active tips and move each one
    for (unsigned tip_index = 0; tip_index < activeTips.size(); tip_index++)
    {
        // Get the vessel network node at the tip
        c_vector<double,DIM> tip_location = this->GetLocationOfCellCentre(activeTips[tip_index]);
        boost::shared_ptr<VascularNode<DIM> > p_node = mCellNodeMap[activeTips[tip_index]];

        // Make sure it is the correct type of cell
        bool is_tip_type = activeTips[tip_index]->GetMutationState()->IsSame(mp_tip_mutation_state);
        bool has_correct_node_type = false;
        if(p_node->GetNumberOfSegments() > 1)
        {
            has_correct_node_type = true;
        }
        // Do the move
        if(is_tip_type && has_correct_node_type)
        {
            // Get the VEGF concentration at the tip
            double tip_concentration = activeTips[tip_index]->GetCellData()->GetItem("VEGF");

            // Get the neighbour data
            unsigned potts_index = this->GetLocationIndexUsingCell(activeTips[tip_index]);
            std::map<std::string, std::vector<double> > nbr_data = GetNeighbourData(potts_index, activeTips[tip_index]);

            // Identify the potts mesh location index corresponding to the direction of highest gradient
            double max_gradient = 0.0; // Assume only movement due to positive gradients
            int max_grandient_index = -1;
            for(unsigned idx=0; idx<nbr_data["Index"].size(); idx++)
            {
                // make sure that tip cell does not go back on itself
                c_vector<double, DIM> neighbour_location = this->rGetMesh().GetNode(unsigned(nbr_data["Index"][idx]))->rGetLocation();
                bool back_on_self = false;

                for (unsigned seg_index = 0; seg_index < p_node->GetNumberOfSegments(); seg_index++)
                {
                    MARK;
                    if(p_node->GetVesselSegment(seg_index)->GetOppositeNode(p_node)->IsCoincident(ChastePoint<DIM>(neighbour_location)))
                    {
                        back_on_self = true;
                        break;
                    }
                    MARK;
                }

                double gradient = (nbr_data["VEGF"][idx] - tip_concentration) / norm_2(tip_location - neighbour_location);
                if(gradient > max_gradient && bool(nbr_data["Occupancy"][idx]) && !back_on_self)
                {
                    max_gradient = gradient;
                    max_grandient_index = unsigned(nbr_data["Index"][idx]);
                }
            }

            // If there is a free location try to move into it
            if(max_grandient_index>-1)
            {
                c_vector<double, DIM> candidate_location = this->rGetMesh().GetNode(max_grandient_index)->rGetLocation();

                // Make a node at the location
                boost::shared_ptr<VascularNode<DIM> > p_new_node;

                /*
                 * Note: the 'tip_location' here is a tip cell already located on the vessel, as decided by
                 * a sprouting rule. This actually creates the sprout, which moves to the candidate_location.
                 */
                MARK;
                boost::shared_ptr<CaVessel<DIM> > p_vessel = mpNetwork->FormSprout(tip_location, candidate_location);
                p_new_node = p_vessel->GetNodeAtOppositeEnd(p_node);
                MARK;

                // Check for anastamosis
                std::set<CellPtr> cells = this->GetCellsUsingLocationIndex(max_grandient_index);
                std::set<CellPtr>::iterator it;
                for (it = cells.begin(); it != cells.end(); ++it)
                {
                    // There is already an EC there
                    if((*it)->GetMutationState()->IsSame(mp_tip_mutation_state) || (*it)->GetMutationState()->IsSame(mp_stalk_mutation_state))
                    {
                        // Deselect self
                        DeselectTipCell(activeTips[tip_index]);

                        // Update vessel network
                        mpNetwork->UpdateNodes();
                        mpNetwork->UpdateVesselNodes();
                        MARK;
                        boost::shared_ptr<VascularNode<DIM> > p_other_node = mpNetwork->DivideVessel(mCellNodeMap[(*it)]->GetVesselSegment(0)->GetVessel(),
                                                                                                     candidate_location);
                        MARK;
                        std::vector<boost::shared_ptr<VascularNode<DIM> > > merge_nodes;
                        merge_nodes.push_back(p_other_node);
                        merge_nodes.push_back(p_new_node);
                        mpNetwork->MergeCoincidentNodes(merge_nodes);
                        MARK;
                        // If we are a tip also de-select at neighbour location
                        if((*it)->GetMutationState()->IsSame(mp_tip_mutation_state))
                        {
                            DeselectTipCell((*it));
                        }

                        // todo need to merge vessels at nodes with connectivity 2 - can do this at the end of the whole method
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
                    this->AddCellUsingLocationIndex(max_grandient_index, p_new_cell); // this doesn't actually add a cell!
                    this->mCells.push_back(p_new_cell); // do it manually here
                    mCellNodeMap[p_new_cell] = mpNetwork->GetNearestNode(this->GetLocationOfCellCentre(p_new_cell));
                    assert(norm_2(this->GetLocationOfCellCentre(p_new_cell) - mCellNodeMap[p_new_cell]->GetLocationVector())<1.e-4);

                    DeselectTipCell(activeTips[tip_index]);
                    mTipCells.push_back(p_new_cell);
                }
            }
        }
    }
}

template<unsigned DIM>
void CaBasedCellPopulationWithVessels<DIM>::DoMigration()
{
    // Loop through active tips and move each one
    std::vector<boost::shared_ptr<Cell> > activeTips = mTipCells;

    for (unsigned tip_index = 0; tip_index < activeTips.size(); tip_index++)
    {
        // Get the vessel network node at the tip
        c_vector<double,DIM> tip_location = this->GetLocationOfCellCentre(activeTips[tip_index]);
        boost::shared_ptr<VascularNode<DIM> > p_node = mCellNodeMap[activeTips[tip_index]];

        // Make sure it is the correct type of cell and Do the move
        if(activeTips[tip_index]->GetMutationState()->IsSame(mp_tip_mutation_state) && p_node->GetNumberOfSegments() == 1)
        {
            // Get the VEGF concentration at the tip
            double tip_concentration = activeTips[tip_index]->GetCellData()->GetItem("VEGF");

            // Get the neighbour data
            unsigned potts_index = this->GetLocationIndexUsingCell(activeTips[tip_index]);
            std::map<std::string, std::vector<double> > nbr_data = GetNeighbourData(potts_index, activeTips[tip_index]);

            // Identify the potts mesh location index corresponding to the direction of highest gradient
            double max_gradient = 0.0; // Assume only movement due to positive gradients
            int max_grandient_index = -1;
            for(unsigned idx=0; idx<nbr_data["Index"].size(); idx++)
            {
                // make sure that tip cell does not go back on itself
                c_vector<double, DIM> neighbour_location = this->rGetMesh().GetNode(unsigned(nbr_data["Index"][idx]))->rGetLocation();
                bool back_on_self = false;
                MARK;
                if(p_node->GetVesselSegment(0)->GetOppositeNode(p_node)->IsCoincident(ChastePoint<DIM>(neighbour_location)))
                {
                    back_on_self = true;
                }
                MARK;
                double gradient = (nbr_data["VEGF"][idx] - tip_concentration) / norm_2(tip_location - neighbour_location);
                if(gradient > max_gradient && bool(nbr_data["Occupancy"][idx]) && !back_on_self)
                {
                    max_gradient = gradient;
                    max_grandient_index = unsigned(nbr_data["Index"][idx]);
                }
            }

            // If there is a free location try to move into it
            if(max_grandient_index>-1)
            {
                c_vector<double, DIM> candidate_location = this->rGetMesh().GetNode(max_grandient_index)->rGetLocation();

                MARK;
                // Make a node at the location
                boost::shared_ptr<VascularNode<DIM> > p_new_node;
                boost::shared_ptr<CaVessel<DIM> > p_vessel = p_node->GetVesselSegment(0)->GetVessel();
                p_new_node = boost::shared_ptr<VascularNode<DIM> > (new VascularNode<DIM>(candidate_location));
                p_node->SetIsMigrating(false);
                MARK;
                p_new_node->SetIsMigrating(true);
                p_vessel->AddSegment(CaVesselSegment<DIM>::Create(p_node, p_new_node));
                mpNetwork->UpdateNodes();
                mpNetwork->UpdateVesselNodes();
                MARK;
                // Check for anastamosis
                std::set<CellPtr> cells = this->GetCellsUsingLocationIndex(max_grandient_index);
                std::set<CellPtr>::iterator it;
                for (it = cells.begin(); it != cells.end(); ++it)
                {
                    // There is already an EC there
                    if((*it)->GetMutationState()->IsSame(mp_tip_mutation_state) || (*it)->GetMutationState()->IsSame(mp_stalk_mutation_state))
                    {
                        // Deselect self
                        DeselectTipCell(activeTips[tip_index]);

                        // Update vessel network
                        mpNetwork->UpdateNodes();
                        mpNetwork->UpdateVesselNodes();
                        MARK;
                        // todo problem line here ... something may be going wrong with the node map
                        boost::shared_ptr<VascularNode<DIM> > p_other_node = mpNetwork->DivideVessel(mCellNodeMap[(*it)]->GetVesselSegment(0)->GetVessel(),
                                                                                                     candidate_location);
                        MARK;
                        std::vector<boost::shared_ptr<VascularNode<DIM> > > merge_nodes;
                        merge_nodes.push_back(p_other_node);
                        merge_nodes.push_back(p_new_node);
                        mpNetwork->MergeCoincidentNodes(merge_nodes);
                        MARK;
                        // If we are a tip also de-select at neighbour location
                        if((*it)->GetMutationState()->IsSame(mp_tip_mutation_state))
                        {
                            DeselectTipCell((*it));
                        }

                        // todo need to merge vessels at nodes with connectivity 2 - can do this at the end of the whole method
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
                    this->AddCellUsingLocationIndex(max_grandient_index, p_new_cell); // this doesn't actually add a cell!
                    this->mCells.push_back(p_new_cell); // do it manually here
                    mCellNodeMap[p_new_cell] = mpNetwork->GetNearestNode(this->GetLocationOfCellCentre(p_new_cell));
                    assert(norm_2(this->GetLocationOfCellCentre(p_new_cell) - mCellNodeMap[p_new_cell]->GetLocationVector())<1.e-4);

                    DeselectTipCell(activeTips[tip_index]);
                    mTipCells.push_back(p_new_cell);
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
    double prob_max = 0.1;
    double half_max_vegf = 0.5;
    std::vector<boost::shared_ptr<Cell> > candidate_tips = std::vector<boost::shared_ptr<Cell> >();
    if(mp_pde_handler)
    {
        for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = this->Begin(); cell_iter != this->End(); ++cell_iter)
        {
            if ((*cell_iter)->GetMutationState()->IsSame(mp_stalk_mutation_state) )
            {
                double vegf_conc = (*cell_iter)->GetCellData()->GetItem("VEGF");
                double prob_tip_selection = prob_max*SimulationTime::Instance()->GetTimeStep()*vegf_conc/(vegf_conc + half_max_vegf);
                if (RandomNumberGenerator::Instance()->ranf() < prob_tip_selection && mCellNodeMap[*cell_iter]->GetNumberOfSegments() > 1)
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
