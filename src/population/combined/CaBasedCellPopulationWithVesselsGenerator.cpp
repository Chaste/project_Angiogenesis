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

#include "CaBasedCellPopulationWithVesselsGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "DefaultCellProliferativeType.hpp"
#include "StalkCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "Owen2011OxygenBasedCellCycleModel.hpp"


template<unsigned DIM>
CaBasedCellPopulationWithVesselsGenerator<DIM>::CaBasedCellPopulationWithVesselsGenerator() :
    mIncludeNormalPopulation(false)
{

}

template<unsigned DIM>
CaBasedCellPopulationWithVesselsGenerator<DIM>::~CaBasedCellPopulationWithVesselsGenerator()
{

}

template<unsigned DIM>
boost::shared_ptr<CaBasedCellPopulationWithVessels<DIM> > CaBasedCellPopulationWithVesselsGenerator<DIM>::CreateCellPopulation(PottsMesh<DIM>& rMesh,
                                                                                                                               boost::shared_ptr<CaVascularNetwork<DIM> > pVascularNetwork)
{
    // create endothelial cell population
    std::vector<unsigned> location_indices;
    for (unsigned index=0; index < rMesh.GetNumNodes(); index++)
    {
        location_indices.push_back(index);
    }

    std::vector<CellPtr> cells;
    MAKE_PTR(DefaultCellProliferativeType, p_diff_type);
    MAKE_PTR(StalkCellMutationState, p_EC_state);
    MAKE_PTR(WildTypeCellMutationState, p_normal_state);
    CellsGenerator<Owen2011OxygenBasedCellCycleModel, DIM> cells_generator;
    cells_generator.GenerateBasicRandom(cells, rMesh.GetNumNodes(), p_diff_type);

    boost::shared_ptr<CaBasedCellPopulationWithVessels<DIM> > cell_population(new CaBasedCellPopulationWithVessels<DIM>(rMesh, cells, location_indices));
    cell_population->SetVesselNetwork(pVascularNetwork);

    // create normal cell population
    for (unsigned index=0; index < rMesh.GetNumNodes(); index++)
    {
        std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double> segment_distance_pair =
                pVascularNetwork->GetNearestSegment(rMesh.GetNode(index)->rGetLocation());

        if (segment_distance_pair.second < 1e-6 || pVascularNetwork->GetDistanceToNearestNode(rMesh.GetNode(index)->rGetLocation()) < 1e-6)
        {
            if (pVascularNetwork->GetDistanceToNearestNode(rMesh.GetNode(index)->rGetLocation()) >= 1e-3)
            {
                boost::shared_ptr<CaVessel<DIM> > pVessel = segment_distance_pair.first->GetVessel();
                pVessel->DivideSegment(rMesh.GetNode(index)->GetPoint());
                pVessel->UpdateNodes();
            }

            pVascularNetwork->UpdateNodes();
            pVascularNetwork->UpdateSegments();
            pVascularNetwork->UpdateVesselNodes();

            // cell is a stalk cell
            cell_population->GetCellUsingLocationIndex(index)->SetMutationState(p_EC_state);
            // associate stalk cell with node
            pVascularNetwork->GetNearestNode(rMesh.GetNode(index)->rGetLocation())->SetCell(cell_population->GetCellUsingLocationIndex(index));
        }
        else
        {
            if (mIncludeNormalPopulation)
            {
                cell_population->GetCellUsingLocationIndex(index)->SetMutationState(p_normal_state);
            }
            else
            {
                cell_population->GetCellUsingLocationIndex(index)->Kill();
            }
        }
    }
    cell_population->RemoveDeadCells();
    return cell_population;
}

template<unsigned DIM>
void CaBasedCellPopulationWithVesselsGenerator<DIM>::SetIncludeNormalCellPopulation(bool include)
{
    mIncludeNormalPopulation = include;
}

//explicit instantiation
template class CaBasedCellPopulationWithVesselsGenerator<2>;
template class CaBasedCellPopulationWithVesselsGenerator<3>;
