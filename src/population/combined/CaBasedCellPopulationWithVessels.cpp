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
                                                                                                     mp_pde_handler()
                                                                                                     {
                                                                                                     }

template<unsigned DIM>
void CaBasedCellPopulationWithVessels<DIM>::SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

template<unsigned DIM>
bool CaBasedCellPopulationWithVessels<DIM>::IsSiteAvailable(unsigned index, CellPtr pCell)
{
    std::vector<unsigned> available_spaces = this->rGetAvailableSpaces();
    if (pCell->GetMutationState()->IsSame(mp_tip_mutation_state))
    {
        if (available_spaces[index] != 0)
        {
            return true;
        }
        else
        {
            // if any cell at location index is a stalk or tip cell then the
            // tip cell can move in to that location
            std::set<CellPtr> cells = this->GetCellsUsingLocationIndex(index);
            std::set<CellPtr>::iterator it;

            for (it = cells.begin(); it != cells.end(); ++it)
            {
                if((*it)->GetMutationState()->IsSame(mp_tip_mutation_state) ||
                        (*it)->GetMutationState()->IsSame(mp_stalk_mutation_state))
                {
                    return true;
                }
            }
            return false;
        }
    }
    else
    {
        return (available_spaces[index] != 0);
    }
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
void CaBasedCellPopulationWithVessels<DIM>::UpdateVascularCellPopulation()
{

    // move already existing tip cells

    // Shuffle the tip cell collection
    RandomNumberGenerator::Instance()->Shuffle(mTipCells);
    std::vector<boost::shared_ptr<Cell> > candidate_tips = mTipCells;

    for (unsigned tip_index = 0; tip_index < candidate_tips.size(); tip_index++)
    {
        // Get the vessel network node at the tip

        c_vector<double,DIM> tip_location = this->GetLocationOfCellCentre(candidate_tips[tip_index]);
        boost::shared_ptr<VascularNode<DIM> > p_node = mpNetwork->GetNearestNode(tip_location);

        // If the type has been reset to stalk then skip it
        // also if the tip cell is still in the middle of a vessel then skip it
        // (tip cells can be selected without forming a sprout immediately)
        if(candidate_tips[tip_index]->GetMutationState()->IsSame(mp_tip_mutation_state) &&
                p_node->GetNumberOfSegments() == 1)
        {
            // Get VEGF concentration at the tip
            double tip_concentration =  candidate_tips[tip_index]->GetCellData()->GetItem("VEGF");

            // Get the location of neighbour points on the potts mesh
            unsigned potts_index = this->GetLocationIndexUsingCell(candidate_tips[tip_index]);
            std::set<unsigned> neighbour_potts_indices = static_cast<PottsMesh<DIM>& >((this->mrMesh)).GetMooreNeighbouringNodeIndices(potts_index);
            std::vector<unsigned> vec_neighbour_potts_indices;
            std::vector<c_vector<double,DIM> > neighbour_locations;
            std::vector<bool> neighbour_occupancy;

            // Get the spatial locations and occupancy of the potts neighbours
            for (std::set<unsigned>::iterator it=neighbour_potts_indices.begin(); it!=neighbour_potts_indices.end(); ++it)
            {
                neighbour_locations.push_back(this->rGetMesh().GetNode(*it)->rGetLocation());
                neighbour_occupancy.push_back(IsSiteAvailable(*it, candidate_tips[tip_index])); //true if site is free
                vec_neighbour_potts_indices.push_back(*it);
            }

            // Get the VEGF concentrations at the neighbours
            std::vector<double> neighbour_concentrations;
            for(unsigned idx=0; idx<neighbour_locations.size(); idx++)
            {
                neighbour_concentrations.push_back(mp_pde_handler->GetPdeSolutionAtPoint(neighbour_locations[idx], "VEGF"));
            }

            // Identify the potts mesh location index corresponding to the direction of highest gradient
            double max_gradient = 0.0; // Assume only migration due to positive gradients
            int max_grandient_index = -1;
            for(unsigned idx=0; idx<neighbour_locations.size(); idx++)
            {
                // make sure that tip cell does not go back on itself
                bool candidateLocationCoincidentWithAttachedVessel = false;
                assert(p_node->GetNumberOfSegments() == 1);
                if(p_node->GetVesselSegment(0)->GetOppositeNode(p_node)->IsCoincident(ChastePoint<DIM>(neighbour_locations[idx])))
                {
                    candidateLocationCoincidentWithAttachedVessel = true;
                }

                double gradient = (neighbour_concentrations[idx] - tip_concentration) / norm_2(tip_location - neighbour_locations[idx]);
                if(gradient > max_gradient && neighbour_occupancy[idx] && !candidateLocationCoincidentWithAttachedVessel)
                {
                    max_gradient = gradient;
                    max_grandient_index = vec_neighbour_potts_indices[idx];
                }
            }

                if(max_grandient_index>-1)
                {
                    c_vector<double, DIM> candidate_location = this->rGetMesh().GetNode(max_grandient_index)->rGetLocation();

                    boost::shared_ptr<CaVessel<DIM> > p_vessel = p_node->GetVesselSegment(0)->GetVessel();
                    // create new node at candidate location
                    boost::shared_ptr<VascularNode<DIM> > p_new_node(new VascularNode<DIM>(candidate_location));
                    // create new segment and add it to the vessel
                    boost::shared_ptr<CaVesselSegment<DIM> > p_new_segment =
                              CaVesselSegment<DIM>::Create(p_node, p_new_node);
                    p_node->SetIsMigrating(false);
                    p_new_node->SetIsMigrating(true);
                    p_vessel->AddSegment(p_new_segment);

                    // if there is a tip or sprout already here do anastamosis
                    std::set<CellPtr> cells = this->GetCellsUsingLocationIndex(max_grandient_index);
                    std::set<CellPtr>::iterator it;

                    for (it = cells.begin(); it != cells.end(); ++it)
                    {
                        if((*it)->GetMutationState()->IsSame(mp_tip_mutation_state))
                        {
                            // deselect self
                            DeselectTipCell(candidate_tips[tip_index]);

                            // todo although not important for this particular cell-based type of model we should set the is migrating flag for
                            // nodes to false when they form anastomoses

                            // deselect at neighbour location
                            DeselectTipCell((*it));

                            // update vessel network
                            boost::shared_ptr<VascularNode<DIM> > p_other_node = mpNetwork->DivideVessel(mpNetwork->GetNearestVessel(candidate_location), candidate_location);
                            std::vector<boost::shared_ptr<VascularNode<DIM> > > merge_nodes;
                            merge_nodes.push_back(p_other_node);
                            merge_nodes.push_back(p_new_node);
                            mpNetwork->MergeCoincidentNodes(merge_nodes);

                            // todo need to merge vessels at nodes with connectivity 2 - can do this at the end of the whole method

                            break;

                        }
                        else if((*it)->GetMutationState()->IsSame(mp_stalk_mutation_state))
                        {
                            // deselect self
                            DeselectTipCell(candidate_tips[tip_index]);

                            // update vessel network
                            boost::shared_ptr<VascularNode<DIM> > p_other_node = mpNetwork->DivideVessel(mpNetwork->GetNearestVessel(candidate_location), candidate_location);
                            std::vector<boost::shared_ptr<VascularNode<DIM> > > merge_nodes;
                            merge_nodes.push_back(p_other_node);
                            merge_nodes.push_back(p_new_node);
                            mpNetwork->MergeCoincidentNodes(merge_nodes);
                            break;
                        }
                    }

                    if(candidate_tips[tip_index]->GetMutationState()->IsSame(mp_tip_mutation_state))
                    {
                        // Create a new cell - N/B: tip cells do not really divide but this is the most convenient way
                        // to implement a discrete snail-trail-type model
                        candidate_tips[tip_index]->ReadyToDivide();
                        CellPtr p_new_cell = candidate_tips[tip_index]->Divide();

                        // Add new cell to the cell population
                        this->AddCellUsingLocationIndex(max_grandient_index, p_new_cell); // this doesn't actually add a cell?!
                        this->mCells.push_back(p_new_cell); // do it manually

                        DeselectTipCell(candidate_tips[tip_index]);
                        mTipCells.push_back(p_new_cell);
                    }

            }
        }

    }

    // implement vessel sprouting

    // If there is a VEGF field use it to select new tip cells
    double prob_max = 0.1;
    double half_max_vegf = 0.5;
    if(mp_pde_handler)
    {
        for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = this->Begin(); cell_iter != this->End(); ++cell_iter)
        {
            if ((*cell_iter)->GetMutationState()->IsSame(mp_stalk_mutation_state) )
            {
                double vegf_conc = (*cell_iter)->GetCellData()->GetItem("VEGF");
                double prob_tip_selection = prob_max*SimulationTime::Instance()->GetTimeStep()*vegf_conc/(vegf_conc + half_max_vegf);
                if (RandomNumberGenerator::Instance()->ranf() < prob_tip_selection
                        && mpNetwork->GetNearestNode(this->GetLocationOfCellCentre((*cell_iter)))->GetNumberOfSegments() > 1)
                {
                    SelectTipCell(*cell_iter);
                }
            }
        }
    }

    // Shuffle the tip cell collection
    RandomNumberGenerator::Instance()->Shuffle(mTipCells);
    candidate_tips = mTipCells;

    // Do sprouting and migration - needs a VEGF field
    if(mp_pde_handler)
    {
        ReplicatableVector pde_solution_repl(mp_pde_handler->GetPdeSolution("VEGF"));
        for (unsigned tip_index = 0; tip_index < candidate_tips.size(); tip_index++)
        {
            // If the type has been reset to stalk then skip it
            if(candidate_tips[tip_index]->GetMutationState()->IsSame(mp_tip_mutation_state))
            {
                // Get VEGF concentration at the tip
                c_vector<double,DIM> tip_location = this->GetLocationOfCellCentre(candidate_tips[tip_index]);
                double tip_concentration =  candidate_tips[tip_index]->GetCellData()->GetItem("VEGF");

                // Get the location of neighbour points on the potts mesh
                unsigned potts_index = this->GetLocationIndexUsingCell(candidate_tips[tip_index]);
                std::set<unsigned> neighbour_potts_indices = static_cast<PottsMesh<DIM>& >((this->mrMesh)).GetMooreNeighbouringNodeIndices(potts_index);
                std::vector<unsigned> vec_neighbour_potts_indices;
                std::vector<c_vector<double,DIM> > neighbour_locations;
                std::vector<bool> neighbour_occupancy;

                // Get the vessel network node at the tip
                boost::shared_ptr<VascularNode<DIM> > p_node = mpNetwork->GetNearestNode(tip_location);

                // Get the spatial locations and occupancy of the potts neighbours
                for (std::set<unsigned>::iterator it=neighbour_potts_indices.begin(); it!=neighbour_potts_indices.end(); ++it)
                {
                    neighbour_locations.push_back(this->rGetMesh().GetNode(*it)->rGetLocation());
                    neighbour_occupancy.push_back(IsSiteAvailable(*it, candidate_tips[tip_index])); //true if site is free
                    vec_neighbour_potts_indices.push_back(*it);
                }

                // Get the VEGF concentrations at the neighbours
                std::vector<double> neighbour_concentrations;
                for(unsigned idx=0; idx<neighbour_locations.size(); idx++)
                {
                    neighbour_concentrations.push_back(mp_pde_handler->GetPdeSolutionAtPoint(neighbour_locations[idx], "VEGF"));
                }

                // Identify the potts mesh location index corresponding to the direction of highest gradient
                double max_gradient = 0.0; // Assume only migration due to positive gradients
                int max_grandient_index = -1;
                for(unsigned idx=0; idx<neighbour_locations.size(); idx++)
                {
                    bool candidateLocationCoincidentWithAttachedVessel = false;
                    for (unsigned seg_index = 0; seg_index < p_node->GetNumberOfSegments(); seg_index++)
                    {
                        if(p_node->GetVesselSegment(seg_index)->GetOppositeNode(p_node)->IsCoincident(ChastePoint<DIM>(neighbour_locations[idx])))
                        {
                            candidateLocationCoincidentWithAttachedVessel = true;
                            break;
                        }
                    }
                    double gradient = (neighbour_concentrations[idx] - tip_concentration) / norm_2(tip_location - neighbour_locations[idx]);
                    if(gradient > max_gradient && neighbour_occupancy[idx] && !candidateLocationCoincidentWithAttachedVessel)
                    {
                        max_gradient = gradient;
                        max_grandient_index = vec_neighbour_potts_indices[idx];
                    }
                }

                // Can form a sprout
                if (p_node->GetNumberOfSegments() > 1)
                {
                    if(max_grandient_index>-1)
                    {
                        c_vector<double, DIM> candidate_location = this->rGetMesh().GetNode(max_grandient_index)->rGetLocation();

                        boost::shared_ptr<CaVessel<DIM> > p_vessel = mpNetwork->FormSprout(tip_location, candidate_location);
                        boost::shared_ptr<VascularNode<DIM> > p_sprout_node = p_vessel->GetNodeAtOppositeEnd(p_node);

                        // if there is a tip or sprout already here do anastamosis
                        std::set<CellPtr> cells = this->GetCellsUsingLocationIndex(max_grandient_index);
                        std::set<CellPtr>::iterator it;

                        for (it = cells.begin(); it != cells.end(); ++it)
                        {
                            if((*it)->GetMutationState()->IsSame(mp_tip_mutation_state))
                            {
                                // deselect self
                                DeselectTipCell(candidate_tips[tip_index]);

                                // todo although not important for this particular cell-based type of model we should set the is migrating flag for
                                // nodes to false when they form anastomoses

                                // deselect at neighbour location
                                DeselectTipCell((*it));

                                // update vessel network
                                boost::shared_ptr<VascularNode<DIM> > p_other_node = mpNetwork->DivideVessel(mpNetwork->GetNearestVessel(candidate_location), candidate_location);
                                std::vector<boost::shared_ptr<VascularNode<DIM> > > merge_nodes;
                                merge_nodes.push_back(p_other_node);
                                merge_nodes.push_back(p_sprout_node);
                                mpNetwork->MergeCoincidentNodes(merge_nodes);

                                // todo need to merge vessels at nodes with connectivity 2 - can do this at the end of the whole method

                                break;

                            }
                            else if((*it)->GetMutationState()->IsSame(mp_stalk_mutation_state))
                            {
                                // deselect self
                                DeselectTipCell(candidate_tips[tip_index]);

                                // update vessel network
                                boost::shared_ptr<VascularNode<DIM> > p_other_node = mpNetwork->DivideVessel(mpNetwork->GetNearestVessel(candidate_location), candidate_location);
                                std::vector<boost::shared_ptr<VascularNode<DIM> > > merge_nodes;
                                merge_nodes.push_back(p_other_node);
                                merge_nodes.push_back(p_sprout_node);
                                mpNetwork->MergeCoincidentNodes(merge_nodes);
                                break;
                            }
                        }

                        if(candidate_tips[tip_index]->GetMutationState()->IsSame(mp_tip_mutation_state))
                        {
                            // Create a new cell - N/B: tip cells do not really divide but this is the most convenient way
                            // to implement a discrete snail-trail-type model
                            candidate_tips[tip_index]->ReadyToDivide();
                            CellPtr p_new_cell = candidate_tips[tip_index]->Divide();

                            // Add new cell to the cell population
                            this->AddCellUsingLocationIndex(max_grandient_index, p_new_cell); // this doesn't actually add a cell?!
                            this->mCells.push_back(p_new_cell); // do it manually

                            DeselectTipCell(candidate_tips[tip_index]);
                            mTipCells.push_back(p_new_cell);
                        }
                    }
                }
            }
        }
    }





}

// Explicit instantiation
template class CaBasedCellPopulationWithVessels<2>;
template class CaBasedCellPopulationWithVessels<3>;
