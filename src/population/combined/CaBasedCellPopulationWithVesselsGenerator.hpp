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

#ifndef CABASEDCELLPOPULATIONWITHVESSELSGENERATOR_HPP_
#define CABASEDCELLPOPULATIONWITHVESSELSGENERATOR_HPP_

#include "SmartPointers.hpp"
#include "CaBasedCellPopulationWithVessels.hpp"
#include "VascularNetwork.hpp"

/**
 * Class which sets up an on lattice cell population consisting of
 * normal cells, tumour cells and vessels. Endothelial cells are
 * associated with the vascular network provided as an argument
 * to the CreateCellPopulation() method.
 */
template<unsigned DIM>
class CaBasedCellPopulationWithVesselsGenerator
{
    /**
     * Whether to fill spaces with wild type cells
     */
    bool mIncludeNormalPopulation;

public:

    /**
     * Constructor
     */
    CaBasedCellPopulationWithVesselsGenerator();

    /**
     * Destructor
     */
    virtual ~CaBasedCellPopulationWithVesselsGenerator();

    /**
     * Return the cell population
     * @param rMesh reference to a potts mesh
     * @param pVascularNetwork the vessel network
     * @return a pointer to the cell population
     */
    boost::shared_ptr<CaBasedCellPopulationWithVessels<DIM> > CreateCellPopulation(PottsMesh<DIM>& rMesh,
                                                                                   boost::shared_ptr<VascularNetwork<DIM> > pVascularNetwork);

    /**
     * Whether to include a normal cell population around the endothelial cell population.
     * @param include include normal cells
     */
    void SetIncludeNormalCellPopulation(bool include);

};

#endif /* CABASEDCELLPOPULATIONWITHVESSELSGENERATOR_HPP_ */
