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
	  mHaematocritCalculator(new ConstantHaematocritSolver<DIM>()),
	  mFlowSolver(new SimpleFlowSolver<DIM>()),
	  mRadiusCalculator(new Alarcon03RadiusCalculator<DIM>()),
	  mMetabolicStimulusCalculator(new Alarcon03MetabolicStimulusCalculator<DIM>()),
	  mImpedanceCalculator(new PoiseuilleImpedanceCalculator<DIM>()),
	  mMechanicalStimulusCalculator(new Alarcon03MechanicalStimulusCalculator<DIM>()),
	  mViscosityCalculator(new Alarcon03ViscosityCalculator<DIM>()),
	  mWallShearStressCalculator(new Alarcon03WallShearStressCalculator<DIM>())
//	  ,
//	  downstreamConductedStimulusCalculator(new NullVesselPropertyLocalCalculator<DIM>()),
//	  upstreamConductedStimulusCalculator(new NullVesselPropertyLocalCalculator<DIM>())
{

}

template<unsigned DIM>
SimpleStructuralAdaptationSolver<DIM>::~SimpleStructuralAdaptationSolver()
{

}

// setters
template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetRadiusCalculator(boost::shared_ptr<Alarcon03RadiusCalculator<DIM> > pCalculator)
{
    mRadiusCalculator = pCalculator;
}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetMetabolicStimulusCalculator(boost::shared_ptr<Alarcon03MetabolicStimulusCalculator<DIM> > pCalculator)
{
    mMetabolicStimulusCalculator = pCalculator;
}


template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetImpedanceCalculator(boost::shared_ptr<PoiseuilleImpedanceCalculator<DIM> > pCalculator)
{
    mImpedanceCalculator = pCalculator;
}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetMechanicalStimulusCalculator(boost::shared_ptr<Alarcon03MechanicalStimulusCalculator<DIM> > pCalculator)
{
    mMechanicalStimulusCalculator = pCalculator;
}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetViscosityCalculator(boost::shared_ptr<Alarcon03ViscosityCalculator<DIM> > pCalculator)
{
    mViscosityCalculator = pCalculator;
}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetWallShearStressCalculator(boost::shared_ptr<Alarcon03WallShearStressCalculator<DIM> > pCalculator)
{
    mWallShearStressCalculator = pCalculator;
}

//template<unsigned DIM>
//void SimpleStructuralAdaptationSolver<DIM>::SetDownstreamConductedStimulusCalculator(boost::shared_ptr<AbstractVesselPropertyLocalCalculator<DIM> > pCalculator)
//{
//    downstreamConductedStimulusCalculator = pCalculator;
//}
//
//template<unsigned DIM>
//void SimpleStructuralAdaptationSolver<DIM>::SetUpstreamConductedStimulusCalculator(boost::shared_ptr<AbstractVesselPropertyLocalCalculator<DIM> > pCalculator)
//{
//    upstreamConductedStimulusCalculator = pCalculator;
//}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetNodePressureCalculator(boost::shared_ptr<SimpleFlowSolver<DIM> > pCalculator)
{
    mFlowSolver = pCalculator;
}

template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::SetHaematocritCalculator(boost::shared_ptr<ConstantHaematocritSolver<DIM> > pCalculator)
{
    mHaematocritCalculator = pCalculator;
}

// method for performing the Calculator
template<unsigned DIM>
void SimpleStructuralAdaptationSolver<DIM>::Iterate(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork)
{

	mRadiusCalculator->SetTimestep(this->GetTimeIncrement());

	mViscosityCalculator->Calculate(vascularNetwork);
	mImpedanceCalculator->Calculate(vascularNetwork);

	mFlowSolver->Implement(vascularNetwork);

    mHaematocritCalculator->Calculate(vascularNetwork);

    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = vascularNetwork->GetVesselSegments();

    for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)
    {
    	segments[segment_index]->SetData("Shrinking Stimulus", 1.79);
    }

	mWallShearStressCalculator->Calculate(vascularNetwork);
	mMetabolicStimulusCalculator->Calculate(vascularNetwork);
	mMechanicalStimulusCalculator->Calculate(vascularNetwork);
//	upstreamConductedStimulusCalculator->Calculate(vascularNetwork);
//	downstreamConductedStimulusCalculator->Calculate(vascularNetwork);
	mRadiusCalculator->Calculate(vascularNetwork);

}

// Explicit instantiation
template class SimpleStructuralAdaptationSolver<2>;
template class SimpleStructuralAdaptationSolver<3>;
