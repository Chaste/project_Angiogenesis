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

#ifndef SIMPLESTRCUTURALADAPATATIONSOLVER_HPP
#define SIMPLESTRCUTURALADAPATATIONSOLVER_HPP

#include "SmartPointers.hpp"
#include "CaVascularNetwork.hpp"
#include "Alarcon03MechanicalStimulusCalculator.hpp"
#include "Alarcon03MetabolicStimulusCalculator.hpp"
#include "Alarcon03RadiusCalculator.hpp"
#include "Alarcon03ViscosityCalculator.hpp"
#include "Alarcon03WallShearStressCalculator.hpp"
#include "PoiseuilleImpedanceCalculator.hpp"
#include "FlowSolver.hpp"
#include "AbstractHaematocritSolver.hpp"
#include "AbstractStructuralAdaptationSolver.hpp"

/**
 * This is a concrete implementation of a structural adaptation solver. It iteratively changes
 * vessel radii in response to a collection of flow based stimuli until the rate of change of
 * the radius is below a specified tolerance.
 */
template<unsigned DIM>
class StructuralAdaptationSolver : public AbstractStructuralAdaptationSolver<DIM>
{

private:

    /**
     * A calculator to determine haematocrit in the network
     */
    boost::shared_ptr<AbstractHaematocritSolver<DIM> > mpHaematocritCalculator;

    /**
     * A solver to calculate flow rates and pressures in the network
     */
    boost::shared_ptr<FlowSolver<DIM> > mpFlowSolver;

    /**
     * A calculator to determine radius changes
     */
    boost::shared_ptr<Alarcon03RadiusCalculator<DIM> > mpRadiusCalculator;

    /**
     * A calculator to determine metabolic stimuli
     */
    boost::shared_ptr<Alarcon03MetabolicStimulusCalculator<DIM> > mpMetabolicStimulusCalculator;

    /**
     * A calculator to determine impedance in vessels
     */
    boost::shared_ptr<PoiseuilleImpedanceCalculator<DIM> > mpImpedanceCalculator;

    /**
     * A calculator to determine mechanical stimuli
     */
    boost::shared_ptr<Alarcon03MechanicalStimulusCalculator<DIM> > mpMechanicalStimulusCalculator;

    /**
     * A calculator to determine fluid effective viscosity in vessels
     */
    boost::shared_ptr<Alarcon03ViscosityCalculator<DIM> > mpViscosityCalculator;

    /**
     * A calculator to determine wall shear stress in vessels
     */
    boost::shared_ptr<Alarcon03WallShearStressCalculator<DIM> > mpWallShearStressCalculator;


public:

    /**
     * Constructor.
     */
    StructuralAdaptationSolver();

    /**
     * Destructor.
     */
    virtual ~StructuralAdaptationSolver();

    /**
     * Factor constructor. Construct a new instance of the class and return a shared pointer to it.
     * @return a pointer to a new instance of the class.
     */
    static boost::shared_ptr<StructuralAdaptationSolver<DIM> > Create();

    /**
     * Perform a single iteration to update the radius and calculators
     */
    virtual void Iterate();

    /**
     * Set the flow calculator
     * @param pSolver the flow solver.
     */
    void SetFlowSolver(boost::shared_ptr<FlowSolver<DIM> > pSolver);

    /**
     * Set the haematocirt solver
     * @param pCalculator the haematocirt calculator.
     */
    void SetHaematocritCalculator(boost::shared_ptr<AbstractHaematocritSolver<DIM> > pCalculator);

    /**
     * Set the impedance calculator
     * @param pCalculator the impedance calculator.
     */
    void SetImpedanceCalculator(boost::shared_ptr<PoiseuilleImpedanceCalculator<DIM> > pCalculator);

    /**
     * Set the metabolic stimulus calculator
     * @param pCalculator the metabolic stimulus calculator.
     */
    void SetMetabolicStimulusCalculator(boost::shared_ptr<Alarcon03MetabolicStimulusCalculator<DIM> > pCalculator);

    /**
     * Set the mechanical stimulus calculator
     * @param pCalculator the mechanical stimulus calculator.
     */
    void SetMechanicalStimulusCalculator(boost::shared_ptr<Alarcon03MechanicalStimulusCalculator<DIM> > pCalculator);

    /**
     * Set the radius calculator
     * @param pCalculator the radius calculator.
     */
    void SetRadiusCalculator(boost::shared_ptr<Alarcon03RadiusCalculator<DIM> > pCalculator);

    /**
     * Set the viscosity calculator
     * @param pCalculator the viscosity calculator.
     */
    void SetViscosityCalculator(boost::shared_ptr<Alarcon03ViscosityCalculator<DIM> > pCalculator);

    /**
     * Set the wall shear stress calculator
     * @param pCalculator the wall shear stress calculator.
     */
    void SetWallShearStressCalculator(boost::shared_ptr<Alarcon03WallShearStressCalculator<DIM> > pCalculator);

};

#endif //SIMPLESTRCUTURALADAPATATIONSOLVER_HPP
