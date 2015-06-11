/*
 * OnLatticeVascularTumourGrowthSimulation.hpp
 *
 *  Created on: 11 Jun 2015
 *      Author: chaste
 */

#ifndef ONLATTICEVASCULARTUMOURGROWTHSIMULATION_HPP_
#define ONLATTICEVASCULARTUMOURGROWTHSIMULATION_HPP_


#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellPopulation.hpp"


template<unsigned DIM>
class OnLatticeVascularTumourGrowthSimulation
{

public:

    OnLatticeVascularTumourGrowthSimulation(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                            bool deleteCellPopulationInDestructor,
                                            bool initialiseCells);

    virtual ~OnLatticeVascularTumourGrowthSimulation();

};

#endif /* ONLATTICEVASCULARTUMOURGROWTHSIMULATION_HPP_ */
