/*
 * OnLatticeVascularTumourCellPopulationGenerator.hpp
 *
 *  Created on: 11 Jun 2015
 *      Author: chaste
 */

#ifndef ONLATTICEVASCULARTUMOURCELLPOPULATIONGENERATOR_HPP_
#define ONLATTICEVASCULARTUMOURCELLPOPULATIONGENERATOR_HPP_

#include "CaBasedCellPopulation.hpp"
#include "CaVascularNetwork.hpp"
#include "SmartPointers.hpp"

/**
 * Class which sets up an on lattice cell population consisting of
 * normal cells, tumour cells and vessels. Endothelial cells are
 * associated with the vascular network provided as an argument
 * to the CreateCellPopulation() method.
 */
template<unsigned DIM>
class OnLatticeVascularTumourCellPopulationGenerator
{

public:

    /*
     * Constructor
     */
    OnLatticeVascularTumourCellPopulationGenerator();

    /*
     * Destructor
     */
    virtual ~OnLatticeVascularTumourCellPopulationGenerator();

    /*
     * Returns cell population
     */
    boost::shared_ptr<CaBasedCellPopulation<DIM> > CreateCellPopulation(PottsMesh<DIM>& rMesh,
                                  boost::shared_ptr<CaVascularNetwork<DIM> > pVascularNetwork);

};

#endif /* ONLATTICEVASCULARTUMOURCELLPOPULATIONGENERATOR_HPP_ */
