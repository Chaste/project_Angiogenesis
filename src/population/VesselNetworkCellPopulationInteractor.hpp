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

#ifndef VesselNetworkCellPopulationInteractor_hpp
#define VesselNetworkCellPopulationInteractor_hpp

#include <vector>
#include <string>
#include "VascularNetwork.hpp"
#include "AbstractCellPopulation.hpp"
#include "AbstractCellMutationState.hpp"
#include "UnitCollections.hpp"

/**
 * The class manages interactions between vessel networks and cell populations. It helps to 'decompartmentalize'
 * each framework by describe interactions in this class rather than individual population classes.
 */
template<unsigned DIM>
class VesselNetworkCellPopulationInteractor
{

protected:

    /**
     * The vessel network can be stored here. The cell population is not stored.
     */
    boost::shared_ptr<VascularNetwork<DIM> > mpNetwork;

public:

    /**
     * Constructor.
     */
    VesselNetworkCellPopulationInteractor();

    /**
     * Destructor.
     */
    virtual ~VesselNetworkCellPopulationInteractor();

    /**
     * Set any cells in a population that overlaps with a vessel to either the tip or stalk cell mutation state.
     * @param cellPopulation the cell population
     * @param pTipMutationState the mutation state for vessel tips
     * @param pStalkState the mutation state for non tip regions
     * @param threshold the max distance from a cell location to vessel centre for labelling
     */
    void LabelVesselsInCellPopulation(AbstractCellPopulation<DIM>& cellPopulation,
                                      boost::shared_ptr<AbstractCellMutationState> pTipMutationState,
                                      boost::shared_ptr<AbstractCellMutationState> pStalkState,
                                      units::quantity<unit::length> threshold = 1.25e-9 * unit::metres);

    /**
     * Divide vessels in the network between cell locations in the cell population by adding new nodes at locations
     * coincident with cells.
     * @param rCellPopulation the cell population
     */
    void PartitionNetworkOverCells(AbstractCellPopulation<DIM>& rCellPopulation, units::quantity<unit::length> threshold = 1.25e-9 * unit::metres);

    /**
     * Remove any cells not overlapping with the vessel network. Does not label the cells.
     * @param rCellPopulation the cell population
     * @param threshold the max distance from a cell location to vessel centre for killing.
     */
    void KillNonVesselOverlappingCells(AbstractCellPopulation<DIM>& rCellPopulation, units::quantity<unit::length> threshold = 1.25e-9 * unit::metres);

    /**
     * Remove any cells overlapping with the vessel network. Does not label the cells.
     * @param rCellPopulation the cell population
     * @param threshold the max distance from a cell location to vessel centre for killing.
     */
    void KillOverlappingVesselCells(AbstractCellPopulation<DIM>& rCellPopulation, units::quantity<unit::length> threshold = 1.25e-9 * unit::metres);

    /**
     * Set the vessel network
     * @param pNetwork the vessel network
     */
    void SetVesselNetwork(boost::shared_ptr<VascularNetwork<DIM> > pNetwork);

};

#endif /* VesselNetworkCellPopulationInteractor_hpp */
