/*
 * OnLatticeVascularTumourCellPopulationGenerator.cpp
 *
 *  Created on: 11 Jun 2015
 *      Author: chaste
 */

#include "OnLatticeVascularTumourCellPopulationGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "DefaultCellProliferativeType.hpp"
#include "StalkCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"

#include "Debug.hpp"

template<unsigned DIM>
OnLatticeVascularTumourCellPopulationGenerator<DIM>::OnLatticeVascularTumourCellPopulationGenerator()
{

}

template<unsigned DIM>
OnLatticeVascularTumourCellPopulationGenerator<DIM>::~OnLatticeVascularTumourCellPopulationGenerator()
{

}

template<unsigned DIM>
boost::shared_ptr<CaBasedCellPopulation<DIM> > OnLatticeVascularTumourCellPopulationGenerator<DIM>::CreateCellPopulation(PottsMesh<DIM>& rMesh,
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
    CellsGenerator<FixedDurationGenerationBasedCellCycleModel, DIM> cells_generator;
    cells_generator.GenerateBasicRandom(cells, rMesh.GetNumNodes(), p_diff_type);

    boost::shared_ptr<CaBasedCellPopulation<DIM> > cell_population(new CaBasedCellPopulation<DIM>(rMesh, cells, location_indices));


    // create normal cell population
    for (unsigned index=0; index < rMesh.GetNumNodes(); index++)
    {

        std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double> segment_distance_pair =
                pVascularNetwork->GetNearestSegment(rMesh.GetNode(index)->rGetLocation());

        if (segment_distance_pair.second < 1e-3 || pVascularNetwork->GetDistanceToNearestNode(rMesh.GetNode(index)->rGetLocation()) < 1e-3)
        {
            if (pVascularNetwork->GetDistanceToNearestNode(rMesh.GetNode(index)->rGetLocation()) >= 1e-3)
            {
                boost::shared_ptr<CaVessel<DIM> > pVessel = segment_distance_pair.first->GetVessel();
                pVessel->DivideSegment(rMesh.GetNode(index)->GetPoint());
            }

            pVascularNetwork->UpdateNodes();
            pVascularNetwork->UpdateSegments();
            pVascularNetwork->UpdateVesselNodes();

            // cell is a stalk cell
            cell_population->GetCellUsingLocationIndex(index)->SetMutationState(p_EC_state);
            // associate stalk cell with node
            pVascularNetwork->GetNearestNode(rMesh.GetNode(index)->rGetLocation())->SetCellPopulation(cell_population);
            pVascularNetwork->GetNearestNode(rMesh.GetNode(index)->rGetLocation())->SetCell(cell_population->GetCellUsingLocationIndex(index));

        }
        else
        {
            cell_population->GetCellUsingLocationIndex(index)->SetMutationState(p_normal_state);
        }
    }



    return cell_population;

        }

//explicit instantiation
//______________________
template class OnLatticeVascularTumourCellPopulationGenerator<2>;
template class OnLatticeVascularTumourCellPopulationGenerator<3>;
