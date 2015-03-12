//
//  Alarcon03MetabolicStimulusCalculator.cpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 07/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#include <stdio.h>
#include "Alarcon03MetabolicStimulusCalculator.hpp"

template<unsigned DIM>
Alarcon03MetabolicStimulusCalculator<DIM>::Alarcon03MetabolicStimulusCalculator() :
Q_ref(6.666*pow(10.0,(-10))),
k_m(0.83),
MaxStimulus(1*pow(10.0,10.0))
{

}

template<unsigned DIM>
Alarcon03MetabolicStimulusCalculator<DIM>::~Alarcon03MetabolicStimulusCalculator()
{

}

template<unsigned DIM>
double Alarcon03MetabolicStimulusCalculator<DIM>::GetQRef()
{
	return Q_ref;
}

template<unsigned DIM>
double Alarcon03MetabolicStimulusCalculator<DIM>::GetKm()
{
	return k_m;
}

template<unsigned DIM>
double Alarcon03MetabolicStimulusCalculator<DIM>::GetMaxStimulus()
{
	return MaxStimulus;
}

template<unsigned DIM>
void Alarcon03MetabolicStimulusCalculator<DIM>::SetQRef(double qref)
{
	assert(Q_ref > 0);
	Q_ref = qref;
}

template<unsigned DIM>
void Alarcon03MetabolicStimulusCalculator<DIM>::SetKm(double km)
{
	assert(k_m >= 0);
	k_m = km;
}

template<unsigned DIM>
void Alarcon03MetabolicStimulusCalculator<DIM>::SetMaxStimulus(double maxstimulus)
{
	assert(MaxStimulus > 1000);
	MaxStimulus = maxstimulus;
}

template<unsigned DIM>
void Alarcon03MetabolicStimulusCalculator<DIM>::Calculate(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork)
{

	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = vascularNetwork->GetVesselSegments();

	for (unsigned segIndex = 0; segIndex < segments.size(); segIndex++)
	{

		double metabolic_stimulus;

		if (fabs(segments[segIndex]->template GetData<double>("Flow Rate")) > 0)
		{
			if(segments[segIndex]->template GetData<double>("Haematocrit Level") > 0)
			{
				metabolic_stimulus = k_m*log10(((Q_ref)/(fabs(segments[segIndex]->template GetData<double>("Flow Rate"))*(segments[segIndex]->template GetData<double>("Haematocrit Level")))) + 1.0);
			}
			else
			{
				// some large number - unphysical if Metabolic stimulus term becomes infinite
				metabolic_stimulus = MaxStimulus;
			}
		}
		else
		{
			metabolic_stimulus = 0;
		}

		segments[segIndex]->SetData("Metabolic Stimulus", metabolic_stimulus);

	}

}

// Explicit instantiation
template class Alarcon03MetabolicStimulusCalculator<2>;
template class Alarcon03MetabolicStimulusCalculator<3>;
