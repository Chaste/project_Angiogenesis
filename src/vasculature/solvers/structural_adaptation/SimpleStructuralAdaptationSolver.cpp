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

#include "SimpleStructuralAdaptationSolver.hpp"

template<unsigned DIM>
SimpleStructuralAdaptationSolver<DIM>::SimpleStructuralAdaptationSolver()
	: StructuralAdaptationSolver<DIM>(),
	  haematocritCalculator(new Alarcon03HaematocritCalculator<DIM>()),
	  nodePressureCalculator(new Alarcon03NodePressureCalculator<DIM>()),
	  radiusCalculator(new Alarcon03RadiusCalculator<DIM>()),
	  metabolicStimulusCalculator(new Alarcon03MetabolicStimulusCalculator<DIM>()),
	  flowRateCalculator(new Alarcon03FlowRateCalculator<DIM>()),
	  flowVelocityCalculator(new Alarcon03FlowVelocityCalculator<DIM>()),
	  impedanceCalculator(new Alarcon03ImpedanceCalculator<DIM>()),
	  mechanicalStimulusCalculator(new Alarcon03MechanicalStimulusCalculator<DIM>()),
	  viscosityCalculator(new Alarcon03ViscosityCalculator<DIM>()),
	  wallShearStressCalculator(new Alarcon03WallShearStressCalculator<DIM>()),
	  downstreamConductedStimulusCalculator(new NullVesselPropertyLocalCalculator<DIM>()),
	  upstreamConductedStimulusCalculator(new NullVesselPropertyLocalCalculator<DIM>())
{
    
}

template<unsigned DIM>
SimpleStructuralAdaptationSolver<DIM>::~SimpleStructuralAdaptationSolver()
{
    
}

// setters
template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetRadiusCalculator(boost::shared_ptr<Alarcon03RadiusCalculator<DIM> > pCalulcator)
{
    radiusCalculator = pCalulcator;
}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetMetabolicStimulusCalculator(boost::shared_ptr<Alarcon03MetabolicStimulusCalculator<DIM> > pCalulcator)
{
    metabolicStimulusCalculator = pCalulcator;
}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetFlowRateCalculator(boost::shared_ptr<Alarcon03FlowRateCalculator<DIM> > pCalulcator)
{
    flowRateCalculator = pCalulcator;
}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetFlowVelocityCalculator(boost::shared_ptr<Alarcon03FlowVelocityCalculator<DIM> > pCalulcator)
{
    flowVelocityCalculator = pCalulcator;
}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetImpedanceCalculator(boost::shared_ptr<PoiseuilleImpedanceCalculator<DIM> > pCalulcator)
{
    impedanceCalculator = pCalulcator;
}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetMechanicalStimulusCalculator(boost::shared_ptr<Alarcon03MechanicalStimulusCalculator<DIM> > pCalulcator)
{
    mechanicalStimulusCalculator = pCalulcator;
}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetViscosityCalculator(boost::shared_ptr<Alarcon03ViscosityCalculator<DIM> > pCalulcator)
{
    viscosityCalculator = pCalulcator;
}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetWallShearStressCalculator(boost::shared_ptr<Alarcon03WallShearStressCalculator<DIM> > pCalulcator)
{
    wallShearStressCalculator = pCalulcator;
}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetDownstreamConductedStimulusCalculator(boost::shared_ptr<AbstractVesselPropertyLocalCalculator<DIM> > pCalulcator)
{
    downstreamConductedStimulusCalculator = pCalulcator;
}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetUpstreamConductedStimulusCalculator(boost::shared_ptr<AbstractVesselPropertyLocalCalculator<DIM> > pCalulcator)
{
    upstreamConductedStimulusCalculator = pCalulcator;
}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetNodePressureCalculator(boost::shared_ptr<SimpleFlowSolver<DIM> > pCalulcator)
{
    nodePressureCalculator = pCalulcator;
}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetHaematocritCalculator(boost::shared_ptr<AbstractHaematocritCalculator<DIM> > pCalulcator)
{
    haematocritCalculator = pCalulcator;
}

// method for performing the Calculator
template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::Iterate(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork)
{
    
	viscosityCalculator->Calculate(vascularNetwork);
	impedanceCalculator->Calculate(vascularNetwork);

    nodePressureCalculator->Calculate(vascularNetwork);
    
    for (int i = 0; i < VN->GetNumberOfVesselsInNetwork(); i++)
    {
        flowRateCalculator->Calculate(VN->GetVessel(i));
        flowVelocityCalculator->Calculate(VN->GetVessel(i));
    }
    
    haematocritCalculator->Calculate(vascularNetwork);
    
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = vascularNetwork->GetVesselSegments();

    for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)
    {
    	segments[segment_index]->SetData("ShrinkingStimulus", 1.79);
    }

	wallShearStressCalculator->Calculate(vascularNetwork);
	metabolicStimulusCalculator->Calculate(vascularNetwork);
	mechanicalStimulusCalculator->Calculate(vascularNetwork);
	upstreamConductedStimulusCalculator->Calculate(vascularNetwork);
	downstreamConductedStimulusCalculator->Calculate(vascularNetwork);
	radiusCalculator->Calculate(vascularNetwork);
    
}

// Explicit instantiation
template class SimpleStructuralAdaptationSolver<2>;
template class SimpleStructuralAdaptationSolver<3>;
