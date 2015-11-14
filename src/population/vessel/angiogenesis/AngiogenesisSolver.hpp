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

#ifndef ANGIOGENESISSOLVER_HPP_
#define ANGIOGENESISSOLVER_HPP_

#include <vector>
#include <string>
#include "SmartPointers.hpp"
#include "CaVascularNetwork.hpp"
#include "Part.hpp"
#include "AbstractGrowthDirectionModifier.hpp"
#include "AbstractSproutingRule.hpp"

/**
 * This class is for simulating modifications to the vessel network due to angiogenesis.
 * The current implementation is based on sprouting angiogenesis.
 */
template<unsigned DIM>
class AngiogenesisSolver
{
    /**
     * The vessel network
     */
    boost::shared_ptr<CaVascularNetwork<DIM> > mpNetwork;

    /**
     * The growth velocity of vessels in angiogenesis simulations
     */
    double mGrowthVelocity;

    /**
     * The end time for solves if used in standalone
     */
    double mEndTime;

    /**
     * The radius in which anastamosis is allowed in angiogenesis simulations
     */
    double mNodeAnastamosisRadius;

    /**
     * The collection of growth direction modifiers
     */
    std::vector<boost::shared_ptr<AbstractGrowthDirectionModifier<DIM> > > mGrowthDirectionModifiers;

    /**
     * The sprouting rule for angiogenesis
     */
    boost::shared_ptr<AbstractSproutingRule<DIM> > mpSproutingRule;

    /**
     * The bounding domain for the vessel network
     */
    boost::shared_ptr<Part<DIM> > mpBoundingDomain;

public:

    /**
     * Constructor.
     */
    AngiogenesisSolver();

    /**
     * Destructor.
     */
    virtual ~AngiogenesisSolver();

    /**
     * Factory constructor method
     * @return a shared pointer to a new solver
     */
    static boost::shared_ptr<AngiogenesisSolver> Create();

    /**
     * Add a growth direction modifier to the collection
     * @param pModifier a growth direction modifier
     */
    void AddGrowthDirectionModifier(boost::shared_ptr<AbstractGrowthDirectionModifier<DIM> > pModifier);

    bool IsSproutingRuleSet();

    /**
     * Set the radius within which anastamosis of vessels is allowed
     * @param radius the radius within which anastamosis of vessels is allowed
     */
    void SetAnastamosisRadius(double radius);

    /**
     * A domain which vessels a not permitted to leave
     * @param pDomain the domain which vessels a not permitted to leave
     */
    void SetBoundingDomain(boost::shared_ptr<Part<DIM> > pDomain);

    /**
     * Set the simulation end time
     * @param time the simulation end time
     */
    void SetEndTime(double time);

    /**
     * Set the base growth velocity for migrating tips
     * @param velocity the velocity of node growth
     */
    void SetGrowthVelocity(double velocity);

    /**
     * Set the rule for managing sprouting
     * @param pSproutingRule the rule for vessel sprouting
     */
    void SetSproutingRule(boost::shared_ptr<AbstractSproutingRule<DIM> > pSproutingRule);

    /**
     * Set the vessel network
     * @param pNetwork the vessel network
     */
    void SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork);

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
     * @param currentDirection the current growth direction
     * @param pNode the moving node on the tip
     */
    virtual c_vector<double, DIM> GetGrowthDirection(c_vector<double, DIM> currentDirection, boost::shared_ptr<VascularNode<DIM> > pNode);

    /**
     * Update the position of all nodes
     * @param speciesLabel the name of the species from which to sample concentration values.
     */
    void UpdateNodalPositions(const std::string& speciesLabel = "Default");
};

#endif /* ANGIOGENESISSOLVER_HPP_ */
