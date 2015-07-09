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
            EXCEPTION("Tip cell is not contained inside tip cell container.");
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

    double pMax = 0.1;
    double halfMaxVEGF = 0.5;

    // iterate through cells and select tip cells based on local concentration of VEGF
    for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = this->Begin();
             cell_iter != this->End();
             ++cell_iter)
    {

        if ((*cell_iter)->GetMutationState()->IsSame(mp_stalk_mutation_state))
        {
            double localVEGFConcentration = 0;

            if (mp_pde_handler != NULL)
            {
//                localVEGFConcentration = mp_pde_handler->GetPdeSolutionAtPoint(GetLocationOfCellCentre(*cell_iter), "VEGF");
                localVEGFConcentration = (*cell_iter)->GetCellData()->GetItem("VEGF");
            }

            double probabilityOfTipSelection = pMax*SimulationTime::Instance()->GetTimeStep()*localVEGFConcentration/(localVEGFConcentration + halfMaxVEGF);

            if (RandomNumberGenerator::Instance()->ranf() < probabilityOfTipSelection)
            {
                SelectTipCell(*cell_iter);
            }

        }

    }
    RandomNumberGenerator::Instance()->Shuffle(mTipCells);

    // decide whether each tip cell migrates
    ReplicatableVector pde_solution_repl(mp_pde_handler->GetPdeSolution("VEGF"));
    for (unsigned tip_index = 0; tip_index < GetNumberOfTipCells(); tip_index++)
    {

        // Get tip location
        c_vector<double,DIM> tip_location = this->GetLocationOfCellCentre(mTipCells[tip_index]);
        double tip_concentration =  mp_pde_handler->GetPdeSolutionAtPoint(tip_location, "VEGF");

        // Get neighbouring Potts mesh node locations
        std::set<unsigned> neighbour_location_indices = this->GetNeighbouringLocationIndices(mTipCells[tip_index]);
        std::vector<unsigned> vec_neighbour_locations;
        std::vector<c_vector<double,DIM> > neighbour_locations;
        std::vector<bool> neighbour_occupancy;

        for (std::set<unsigned>::iterator it=neighbour_location_indices.begin(); it!=neighbour_location_indices.end(); ++it)
        {
            neighbour_locations.push_back(this->rGetMesh().GetNode(*it)->rGetLocation());
            neighbour_occupancy.push_back(this->IsSiteAvailable(*it, mTipCells[tip_index]));
            vec_neighbour_locations.push_back(*it);
        }

        // Get the VEGF concentrations at the neighbouring mesh node locations
        std::vector<double> neighbour_concentrations(2* DIM);
        for(unsigned idx=0; idx<neighbour_locations.size(); idx++)
        {
            neighbour_concentrations[idx] = mp_pde_handler->GetPdeSolutionAtPoint(neighbour_locations[idx], "VEGF");
        }

        // Identify the direction with highest gradient that is also free
        double max_gradient = 0.0; // what about negative grads
        int max_grandient_index = -1;
        for(unsigned idx=0; idx<neighbour_locations.size(); idx++)
        {
            double gradient = (neighbour_concentrations[idx] - tip_concentration) / norm_2(tip_location - neighbour_locations[idx]);
            if(gradient > max_gradient && neighbour_occupancy[idx])
            {
                max_gradient = gradient;
                max_grandient_index = vec_neighbour_locations[idx];
            }
        }

        // Sprout into the neighbour
        if(max_gradient > 1.0 && max_grandient_index>-1)
        {
            // Sprout...but we need a sprout functionality first
        }
    }
}


//template<unsigned DIM>
//void CaBasedCellPopulationWithVessels<DIM>::AsscoiateVesselNetworkWithCells(
//                      boost::shared_ptr<AbstractCellMutationState> pStalkCellMutatationState,
//                      boost::shared_ptr<AbstractCellMutationState> pTipCellMutatationState)
//{
//    if(!mpNetwork)
//    {
//        EXCEPTION("A vessel network has not been assigned to the population.");
//    }
//
//    // TODO the vessel network should be fully snapped to the potts mesh. This should be tested for here.
//    // For now create nodes on the vessel network wherever potts mesh nodes intersect segments.
//    for (unsigned index=0; index < this->rGetMesh().GetNumNodes(); index++)
//    {
//        c_vector<double, DIM> node_location = this->rGetMesh().GetNode(index)->rGetLocation();
//        std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double> segment_distance_pair =
//                mpNetwork->GetNearestSegment(node_location);
//
//        // If the mesh node is on the segment, create a vessel node at the location (if there isn't one
//        // already there) and assign a cell to it.
//        if (segment_distance_pair.second < 1e-6)
//        {
//            boost::shared_ptr<CaVessel<DIM> > pVessel = segment_distance_pair.first->GetVessel();
//            boost::shared_ptr<VascularNode<DIM> > pNode = pVessel->DivideSegment(this->rGetMesh().GetNode(index)->GetPoint());
//
//            // If there is a cell at this location assign it stalk or tip mutation states, if they have been defined.
//            if(this->IsCellAttachedToLocationIndex(index))
//            {
//                CellPtr pCell = this->GetCellUsingLocationIndex(index);
//                pNode->SetCell(pCell);
//
//                if(pNode->GetNumberOfSegments()==1 && pTipCellMutatationState)
//                {
//                    pCell->SetMutationState(pTipCellMutatationState);
//                }
//                else
//                {
//                    pCell->SetMutationState(pStalkCellMutatationState);
//                }
//            }
//            else
//            {
//                // For now, throw an exception, in future add cells at vessel node locations
//                EXCEPTION("No cell defined at the location of a vessel node.");
//            }
//            mpNetwork->UpdateNodes();
//            mpNetwork->UpdateSegments();
//            mpNetwork->UpdateVesselNodes();
//        }
//    }
//}

// Explicit instantiation
template class CaBasedCellPopulationWithVessels<2>;
template class CaBasedCellPopulationWithVessels<3>;
