//
//  Alarcon03MechanicalStimulusCalculator.cpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 07/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#include <stdio.h>
#include "Alarcon03MechanicalStimulusCalculator.hpp"
#include "Debug.hpp"

template<unsigned DIM>
Alarcon03MechanicalStimulusCalculator<DIM>::Alarcon03MechanicalStimulusCalculator()
	:Tau_ref(0.05),
	Tau_P(0.0)
{
    
}

template<unsigned DIM>
Alarcon03MechanicalStimulusCalculator<DIM>::~Alarcon03MechanicalStimulusCalculator()
{
    
}

template<unsigned DIM>
double Alarcon03MechanicalStimulusCalculator<DIM>::GetTauP()
{
    return Tau_P;
}

template<unsigned DIM>
double Alarcon03MechanicalStimulusCalculator<DIM>::GetTauRef()
{
    return Tau_ref;
}

template<unsigned DIM>
void Alarcon03MechanicalStimulusCalculator<DIM>::SetTauRef(double tau_ref)
{
    assert(tau_ref > 0);
    Tau_ref = tau_ref;
}

template<unsigned DIM>
void Alarcon03MechanicalStimulusCalculator<DIM>::Calculate(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork)
{
    
	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = vascularNetwork->GetVesselSegments();

	for (unsigned segIndex = 0; segIndex < segments.size(); segIndex++)
	{
	    double mechanical_stimulus;

	    // get average pressure on segment in Pa
	    double average_pressure = (segments[segIndex]->GetNode(0)->template GetData<double>("Pressure") + segments[segIndex]->GetNode(1)->template GetData<double>("Pressure"))/2;

	    // convert to mmHg
	    average_pressure *= 760/(1.01*pow(10.0,5));

	    if (log10(average_pressure) < 1)
	    {
	        Tau_P = 1.4; //
	    }
	    else
	    {
		    // tau_p calculated in pascals
	    	// factor of 0.1 introduced in order to convert original expression
	    	// (calculated in units of dyne/cm^2) to pascals
	    	Tau_P = 0.1*(100.0 - 86.0*pow(exp(-5.0*log10(log10(average_pressure))), 5.4));
	    }

	    mechanical_stimulus = log10((segments[segIndex]->template GetData<double>("Wall Shear Stress") + Tau_ref)/Tau_P);
	    segments[segIndex]->SetData("Mechanical Stimulus", mechanical_stimulus);

	}

}

// Explicit instantiation
template class Alarcon03MechanicalStimulusCalculator<2>;
template class Alarcon03MechanicalStimulusCalculator<3>;
