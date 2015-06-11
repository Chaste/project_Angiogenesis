/*
 * OnLatticeVascularTumourCellPopulationGenerator.cpp
 *
 *  Created on: 11 Jun 2015
 *      Author: chaste
 */

#include "OnLatticeVascularTumourCellPopulationGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"

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

    std::vector<unsigned> location_indices;
    unsigned numEndothelialCells = 0;

    for (unsigned index=0; index < rMesh.GetNumNodes(); index++)
    {
        std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double> segment_distance_pair =
                                    pVascularNetwork->GetNearestSegment(rMesh.GetNode(index)->rGetLocation());

        if (segment_distance_pair.second < 1e-6 || pVascularNetwork->GetDistanceToNearestNode(rMesh.GetNode(index)->rGetLocation()) < 1e-6)
        {
            location_indices.push_back(index);
            numEndothelialCells++;
        }
    }

    std::vector<CellPtr> cells;
    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
    CellsGenerator<FixedDurationGenerationBasedCellCycleModel, DIM> cells_generator;
    cells_generator.GenerateBasicRandom(cells, numEndothelialCells, p_diff_type);

    boost::shared_ptr<CaBasedCellPopulation<DIM> > cell_population(new CaBasedCellPopulation<DIM>(rMesh, cells, location_indices));

    return cell_population;

}

//explicit instantiation
//______________________
template class OnLatticeVascularTumourCellPopulationGenerator<2>;
template class OnLatticeVascularTumourCellPopulationGenerator<3>;
