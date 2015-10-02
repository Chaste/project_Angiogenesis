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

#ifndef ABSTRACTANGIOGENESISSOLVER_HPP_
#define ABSTRACTANGIOGENESISSOLVER_HPP_

#include <vector>
#include <string>
#include "SmartPointers.hpp"
#include "CaVascularNetwork.hpp"
#include "AbstractHybridSolver.hpp"
#include "AbstractGrowthDirectionModifier.hpp"
#include "AbstractSproutingRule.hpp"
#include "SimpleFlowSolver.hpp"
#include "SimpleStructuralAdaptationSolver.hpp"
#include "Part.hpp"

template<unsigned DIM>
class AbstractAngiogenesisSolver
{
    boost::shared_ptr<CaVascularNetwork<DIM> > mpNetwork;

    double mGrowthVelocity;

    double mEndTime;

    unsigned mOutputFrequency;

    std::string mOutputDirectory;

    double mNodeAnastamosisRadius;

    std::vector<boost::shared_ptr<AbstractHybridSolver<DIM> > > mPdeSolvers;

    std::vector<boost::shared_ptr<AbstractGrowthDirectionModifier<DIM> > > mGrowthDirectionModifiers;

    boost::shared_ptr<AbstractSproutingRule<DIM> > mpSproutingRule;

    boost::shared_ptr<SimpleFlowSolver<DIM> > mpFlowSolver;

    boost::shared_ptr<SimpleStructuralAdaptationSolver<DIM> > mpStructuralAdaptationSolver;

    boost::shared_ptr<Part<DIM> > mpBoundingDomain;

public:

    /**
     * Constructor.
     * @param pNetwork the network to perform angiogenesis on
     */
    AbstractAngiogenesisSolver(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork);

    /**
     * Destructor.
     */
    virtual ~AbstractAngiogenesisSolver();

    /* Factory constructor method
     * @return a shared pointer to a new solver
     */
    static boost::shared_ptr<AbstractAngiogenesisSolver> Create(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork);

    /**
     * Add a growth direction modifier to the collection
     */
    void AddGrowthDirectionModifier(boost::shared_ptr<AbstractGrowthDirectionModifier<DIM> > pModifier);

    /**
     * Add a PDE solver to the collection
     */
    void AddPdeSolver(boost::shared_ptr<AbstractHybridSolver<DIM> > pPdeSolver);

    void SetAnastamosisRadius(double radius);

    void SetBoundingDomain(boost::shared_ptr<Part<DIM> > pDomain);

    /**
     * Return the current PDE solvers
     */
    std::vector<boost::shared_ptr<AbstractHybridSolver<DIM> > > GetPdeSolvers();

    /**
     * Set the simulation end time
     */
    void SetEndTime(double time);

    /**
     * Set the flow solver for the network
     */
    void SetFlowSolver(boost::shared_ptr<SimpleFlowSolver<DIM> > pFlowSolver);

    /**
     * Set the base growth velocity for migrating tips
     */
    void SetGrowthVelocity(double velocity);

    /**
     * Set the outputdirectory for results
     */
    void SetOutputDirectory(const std::string& rDirectory);

    /**
     * Set the results output frequency
     */
    void SetOutputFrequency(unsigned frequency);

    /**
     * Set the rule for managing sprouting
     */
    void SetSproutingRule(boost::shared_ptr<AbstractSproutingRule<DIM> > pSproutingRule);

    /**
     * Set the structural adaptation solver for the network
     */
    void SetStructuralAdaptationSolver(boost::shared_ptr<SimpleStructuralAdaptationSolver<DIM> > pStructuralAdaptationSolver);

    /**
     * Increment one step in time
     */
    void Increment();

    /**
     * Run until the specified end time
     */
    void Run();


protected:

    /**
     * Identify and grow sprouts
     */
    void DoSprouting();

    /**
     * Do the anastamosis step
     */
    void DoAnastamosis();

    /**
     * Get the growth direction for moving tips
     */
    virtual c_vector<double, DIM> GetGrowthDirection(c_vector<double, DIM> currentDirection, boost::shared_ptr<VascularNode<DIM> > pNode);

    /**
     * Update the position of all nodes
     */
    void UpdateNodalPositions(const std::string& speciesLabel = "Default");
};

#endif /* ABSTRACTANGIOGENESISSOLVER_HPP_ */
