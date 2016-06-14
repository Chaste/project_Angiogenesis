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

#include "Exception.hpp"
#include "MechanicalStimulusCalculator.hpp"

template<unsigned DIM>
MechanicalStimulusCalculator<DIM>::MechanicalStimulusCalculator() : AbstractVesselNetworkCalculator<DIM>(),
    mTauRef(0.05 * unit::pascals),
    mTauP(0.0 * unit::pascals)
{

}

template<unsigned DIM>
MechanicalStimulusCalculator<DIM>::~MechanicalStimulusCalculator()
{

}

template<unsigned DIM>
units::quantity<unit::pressure> MechanicalStimulusCalculator<DIM>::GetTauP()
{
    return mTauP;
}

template<unsigned DIM>
void MechanicalStimulusCalculator<DIM>::SetTauRef(units::quantity<unit::pressure> TauRef)
{
    mTauRef = TauRef;
}

template<unsigned DIM>
void MechanicalStimulusCalculator<DIM>::Calculate()
{

    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = this->mpNetwork->GetVesselSegments();

    for (unsigned idx = 0; idx < segments.size(); idx++)
    {
        // get average pressure in a segment. It is stored in pascal, so is converted to mmHg for the calculation.
        double node0_pressure = segments[idx]->GetNode(0)->GetFlowProperties()->GetPressure();
        double node1_pressure = segments[idx]->GetNode(1)->GetFlowProperties()->GetPressure();

        // Conversion to mmHg
        double average_pressure = (node0_pressure + node1_pressure) * 760.0 / (2.0 * 1.01 * pow(10.0, 5));

        // The calculation does not work for pressures less than 1 mmHg, so we specify a cut-off value of TauP for lower
        // pressures.
        if (log10(average_pressure) < 1.0)
        {
            mTauP = 1.4 * unit::pascals;
        }
        else
        {
            // tau_p calculated in pascals
            // factor of 0.1 introduced in order to convert original expression (calculated in units of dyne/cm^2) to pascals
            mTauP = 0.1 * (100.0 - 86.0 * pow(exp(-5.0 * log10(log10(average_pressure))), 5.4)) * unit::pascals;
        }

        units::quantity<unit::rate> mechanical_stimulus = log10((segments[idx]->GetFlowProperties()->GetWallShearStress() + mTauRef) / mTauP) * unit::reciprocal_seconds;
        segments[idx]->GetFlowProperties()->SetStimulus(segments[idx]->GetFlowProperties()->GetStimulus() + mechanical_stimulus);
    }
}

// Explicit instantiation
template class MechanicalStimulusCalculator<2> ;
template class MechanicalStimulusCalculator<3> ;