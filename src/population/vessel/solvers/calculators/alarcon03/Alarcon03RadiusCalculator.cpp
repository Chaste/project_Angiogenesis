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

#include "Alarcon03RadiusCalculator.hpp"

template<unsigned DIM>
Alarcon03RadiusCalculator<DIM>::Alarcon03RadiusCalculator()
	: mMinRadius(1e-6),
	  mMaxRadius(50e-6),
	  mTimeStep(0.0001)
{
    
}


template<unsigned DIM>
Alarcon03RadiusCalculator<DIM>::~Alarcon03RadiusCalculator()
{
    
}

template<unsigned DIM>
void Alarcon03RadiusCalculator<DIM>::SetMinRadius(double minRadius)
{
	mMinRadius = minRadius;
}

template<unsigned DIM>
void Alarcon03RadiusCalculator<DIM>::SetMaxRadius(double maxRadius)
{
    mMaxRadius = maxRadius;
}

template<unsigned DIM>
void Alarcon03RadiusCalculator<DIM>::SetTimestep(double dt)
{
    mTimeStep = dt;
}

template<unsigned DIM>
void Alarcon03RadiusCalculator<DIM>::Calculate(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork)

{
	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = vascularNetwork->GetVesselSegments();

	for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)
	{

		double mechanical_stimulus = segments[segment_index]->GetFlowProperties()->GetMechanicalStimulus();
		double metabolic_stimulus = segments[segment_index]->GetFlowProperties()->GetMetabolicStimulus();
		double upstream_conducted_stimulus = segments[segment_index]->GetFlowProperties()->GetUpstreamConductedStimulus();
		double downstream_stimulus = segments[segment_index]->GetFlowProperties()->GetDownstreamConductedStimulus();
		double shrinking_stimulus = segments[segment_index]->GetFlowProperties()->GetShrinkingStimulus();

		double total_stimulus = mechanical_stimulus + metabolic_stimulus + upstream_conducted_stimulus + downstream_stimulus - shrinking_stimulus;

		double radius = segments[segment_index]->GetRadius();

//		std::cout << "step "<< mTimeStep << " stim "<< mechanical_stimulus << std::endl;
		radius *= 1.0 + mTimeStep*total_stimulus;

		if (radius > mMaxRadius)
		{
			radius = mMaxRadius;
		}
		if (radius < mMinRadius)
		{
			radius = mMinRadius;
		}

		segments[segment_index]->SetRadius(radius);
	}
}

// Explicit instantiation
template class Alarcon03RadiusCalculator<2>;
template class Alarcon03RadiusCalculator<3>;