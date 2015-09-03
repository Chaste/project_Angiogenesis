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

#ifndef CABASEDCELLPOPULATIONWITHVESSELS_HPP_
#define CABASEDCELLPOPULATIONWITHVESSELS_HPP_

#include "CaBasedCellPopulation.hpp"
#include "PottsMesh.hpp"
#include "VertexMesh.hpp"
#include "AbstractCaUpdateRule.hpp"
#include "StalkCellMutationState.hpp"
#include "TipCellMutationState.hpp"

template<unsigned DIM>
class CellBasedPdeHandler; // circular definition

#include "CaVascularNetwork.hpp"

template<unsigned DIM>
class AbstractCaUpdateRule; // Circular definition

template<unsigned DIM>
class CaBasedCellPopulationWithVessels : public CaBasedCellPopulation<DIM>
{
    friend class TestCaBasedCellPopulationWithVessels;

private:

    boost::shared_ptr<CaVascularNetwork<DIM> > mpNetwork;

    std::vector<boost::shared_ptr<Cell> > mTipCells;

    boost::shared_ptr<TipCellMutationState> mp_tip_mutation_state;

    boost::shared_ptr<StalkCellMutationState> mp_stalk_mutation_state;

    boost::shared_ptr<CellBasedPdeHandler<DIM> > mp_pde_handler;

public:

    /**
     * Create a new cell population facade from a mesh, a vector of location indices
     * and a collection of cells.
     *
     * There must be precisely one CellPtr for each entry of the locationIndices vector.
     *
     * @param rMesh reference to a PottsMesh
     * @param rCells reference to a vector of CellPtrs
     * @param locationIndices a vector of location indices that correspond to real cells
     * @param latticeCarryingCapacity an optional parameter to allow more than one cell per site
     * @param deleteMesh set to true if you want the cell population to free the mesh memory on destruction
     *                   (defaults to false)
     * @param validate whether to validate the cell population when it is created (defaults to false as not used in CA simulations)
     */
    CaBasedCellPopulationWithVessels(PottsMesh<DIM>& rMesh,
                                     std::vector<CellPtr>& rCells,
                                     const std::vector<unsigned> locationIndices,
                                     unsigned latticeCarryingCapacity=1u,
                                     bool deleteMesh=false,
                                     bool validate=false);

    /**
     * Add a vessel network to the population
     *
     * @param pNetwork
     */
    void SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> >pNetwork);

    /**
     * Select a cell to take on tip cell mutation
     */
    void SelectTipCell(boost::shared_ptr<Cell> pCell);

    /**
     * Select a cell to take on stalk cell mutation
     */
    void DeselectTipCell(boost::shared_ptr<Cell> pCell);

    /**
     * Return number of tip cells in population.
     */
    unsigned GetNumberOfTipCells();

    /**
     * Find if a given node has space available.
     *
     * @param index  The global index of a specified node.
     * @param pCell  The cell wanting to divide into the lattice site.
     *
     * @return whether the node is an empty site
     */
    virtual bool IsSiteAvailable(unsigned index, CellPtr pCell);

    /**
     * Return vector of tip cells.
     */
    std::vector<boost::shared_ptr<Cell> > GetTipCells();

    void AddPdeHandler(boost::shared_ptr<CellBasedPdeHandler<DIM> > pde_handler);

    /**
     * Method to update endothelial cell population
     * \\TODO break this up into smaller pieces and make sure those smaller pieces are in-fitting
     * with methods in CaBasedCellPopulation
     */
    void UpdateVascularCellPopulation();

};

#endif /*CABASEDCELLPOPULATIONWITHVESSELS_HPP_*/
