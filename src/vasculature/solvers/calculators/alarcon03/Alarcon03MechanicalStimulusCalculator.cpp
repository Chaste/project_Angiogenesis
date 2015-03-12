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

#include "Alarcon03MechanicalStimulusCalculator.hpp"

template<unsigned DIM>
Alarcon03MechanicalStimulusCalculator<DIM>::Alarcon03MechanicalStimulusCalculator()
	:mTauRef(0.05),
	 mTauP(0.0)
{
    
}

template<unsigned DIM>
Alarcon03MechanicalStimulusCalculator<DIM>::~Alarcon03MechanicalStimulusCalculator()
{
    
}

template<unsigned DIM>
double Alarcon03MechanicalStimulusCalculator<DIM>::GetTauP()
{
    return mTauP;
}

template<unsigned DIM>
double Alarcon03MechanicalStimulusCalculator<DIM>::GetTauRef()
{
    return mTauRef;
}

template<unsigned DIM>
void Alarcon03MechanicalStimulusCalculator<DIM>::SetTauRef(double TauRef)
{
	if(TauRef <= 0.0)
	{
		EXCEPTION("Reference Wall Shear Stress must be positive.");
	}

    mTauRef = TauRef;
}

template<unsigned DIM>
void Alarcon03MechanicalStimulusCalculator<DIM>::Calculate(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork)
{
    
	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = vascularNetwork->GetVesselSegments();

	for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)
	{
	    // get average pressure in segment. It is stored in pascal, so is converted to mmHg for the calculation.
		double node0_pressure = segments[segment_index]->GetNode(0)->template GetData<double>("Pressure");
		double node1_pressure = segments[segment_index]->GetNode(1)->template GetData<double>("Pressure");

	    double average_pressure = (node0_pressure + node1_pressure)*760.0/(2.0*1.01*pow(10.0,5));

	    // The calculation does not work for pressures less than 1 mmHg, so we specify a cut-off value of TauP for lower
	    // pressures.
	    if (log10(average_pressure) < 1.0)
	    {
	    	mTauP = 1.4;
	    }
	    else
	    {
		    // tau_p calculated in pascals
	    	// factor of 0.1 introduced in order to convert original expression
	    	// (calculated in units of dyne/cm^2) to pascals
	    	mTauP = 0.1*(100.0 - 86.0*pow(exp(-5.0*log10(log10(average_pressure))), 5.4));
	    }

	    double wall_shear_stress = segments[segment_index]->template GetData<double>("Wall Shear Stress");
	    double mechanical_stimulus = log10((wall_shear_stress + mTauRef)/mTauP);
	    segments[segment_index]->SetData("Mechanical Stimulus", mechanical_stimulus);
	}
}

// Explicit instantiation
template class Alarcon03MechanicalStimulusCalculator<2>;
template class Alarcon03MechanicalStimulusCalculator<3>;
