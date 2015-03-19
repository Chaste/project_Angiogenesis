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

#include "Alarcon03MetabolicStimulusCalculator.hpp"

template<unsigned DIM>
Alarcon03MetabolicStimulusCalculator<DIM>::Alarcon03MetabolicStimulusCalculator()
	:mQRef(6.666*pow(10.0, -10)),
	mKm(0.83),
	mMaxStimulus(1*pow(10.0, 10.0))
{

}

template<unsigned DIM>
Alarcon03MetabolicStimulusCalculator<DIM>::~Alarcon03MetabolicStimulusCalculator()
{

}

template<unsigned DIM>
double Alarcon03MetabolicStimulusCalculator<DIM>::GetQRef()
{
	return mQRef;
}

template<unsigned DIM>
double Alarcon03MetabolicStimulusCalculator<DIM>::GetKm()
{
	return mKm;
}

template<unsigned DIM>
double Alarcon03MetabolicStimulusCalculator<DIM>::GetMaxStimulus()
{
	return mMaxStimulus;
}

template<unsigned DIM>
void Alarcon03MetabolicStimulusCalculator<DIM>::SetQRef(double qRef)
{
	if(mQRef <= 0.0)
	{
		EXPCETION("Reference flow rate must be positive.");
	}
	mQRef = qRef;

}

template<unsigned DIM>
void Alarcon03MetabolicStimulusCalculator<DIM>::SetKm(double km)
{
	if(km <= 0.0)
	{
		EXPCETION("Parameter km must be positive.");
	}
	mKm = km;
}

template<unsigned DIM>
void Alarcon03MetabolicStimulusCalculator<DIM>::SetMaxStimulus(double maxStimulus)
{
	if(mMaxStimulus <= 0.0)
	{
		EXPCETION("Max Stimulus parameter must be positive");
	}

	mMaxStimulus = maxStimulus;
}

template<unsigned DIM>
void Alarcon03MetabolicStimulusCalculator<DIM>::Calculate(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork)
{

	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = vascularNetwork->GetVesselSegments();
	for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)
	{

		double metabolic_stimulus;
		double haematocrit = segments[segment_index]->template GetData<double>("Haematocrit");
		double flow_rate = segments[segment_index]->template GetData<double>("Absolute Flow Rate");

		if (flow_rate > 0.0)
		{
			if(haematocrit > 0.0)
			{
				metabolic_stimulus = mKm*log10(mQRef/(flow_rate*haematocrit) + 1.0);
			}
			else
			{
				// todo explore alternatives to sticking in an arbitary large number.
				// some large number - unphysical if Metabolic stimulus term becomes infinite
				metabolic_stimulus = mMaxStimulus;
			}
		}
		else
		{
			metabolic_stimulus = 0;
		}

		segments[segment_index]->SetData("Metabolic Stimulus", metabolic_stimulus);
	}
}

// Explicit instantiation
template class Alarcon03MetabolicStimulusCalculator<2>;
template class Alarcon03MetabolicStimulusCalculator<3>;
