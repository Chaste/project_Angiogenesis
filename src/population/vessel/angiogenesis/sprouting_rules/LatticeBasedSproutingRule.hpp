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

#ifndef LATTICEBASEDSPROUTINGRULE_HPP_
#define LATTICEBASEDSPROUTINGRULE_HPP_

#include <vector>
#include <string>
#include "VascularNode.hpp"
#include "SmartPointers.hpp"
#include "AbstractSproutingRule.hpp"
#include "RegularGrid.hpp"
#include "AbstractRegularGridHybridSolver.hpp"

/**
 * A concrete sprouting rule for lattice based simulations based on
 * a model described in Owen et al. 2011. Default Parameter values are taken
 * from that study.
 */
template<unsigned DIM>
class LatticeBasedSproutingRule : public AbstractSproutingRule<DIM>
{

protected:

    /**
     * The grid for the vessel simulation
     */
    boost::shared_ptr<RegularGrid<DIM> > mpGrid;

    /**
     * A hybrid solver for calculating the VEGF field
     */
    boost::shared_ptr<AbstractRegularGridHybridSolver<DIM> > mpSolver;

    /**
     * Cell motility for random walks
     */
    double mCellMotility;

    /**
     * Cell chemotactic sensitivity
     */
    double mCellChemotacticParameter;

    /**
     * Cell sprouting probability in high vegf regions
     */
    double mMaxSproutingProbability;

    /**
     * VEGF concentration for half max sprouting probability
     */
    double mHalfMaxSproutingProbability;

public:

    /**
     * Constructor.
     */
    LatticeBasedSproutingRule();

    /**
     * Destructor.
     */
    virtual ~LatticeBasedSproutingRule();

    /**
     * Construct a new instance of the class and return a shared pointer to it.
     * @return a pointer to a new instance of the class
     */
    static boost::shared_ptr<LatticeBasedSproutingRule<DIM> > Create();

    /**
     * Set the grid for the vessel network
     * @param pGrid the grid for the vessel network
     */
    void SetGrid(boost::shared_ptr<RegularGrid<DIM> > pGrid);

    /**
     * Set the grid for the vessel network
     * @param pGrid the grid for the vessel network
     */
    void SetHybridSolver(boost::shared_ptr<AbstractRegularGridHybridSolver<DIM> > pSolver);


    void SetCellMotilityParameter(double cellMotility);

    void SetCellChemotacticParameter(double cellChemotacticParameter);

    void SetMaxSproutingProbability(double maxSproutingProbability);

    void SetHalfMaxSproutingProbability(double halfMaxSproutingProbability);

    /**
     * Calculate direction for each node that will sprout
     * @param rNodes nodes to calculate directions for
     * @return a vector of bools which true for nodes that will sprout
     */
    virtual std::vector<c_vector<double, DIM> > GetSproutDirections(const std::vector<boost::shared_ptr<VascularNode<DIM> > >& rNodes);

};

#endif /* RANDOMLATTICEBASEDSPROUTINGRULE_HPP_ */
