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

#include "PoiseuilleImpedanceCalculator.hpp"
#include "VesselSegment.hpp"
#include "VasculatureData.hpp"
#include "LinearSystem.hpp"
#include "Exception.hpp"
#include "ReplicatableVector.hpp"
#include "MathsCustomFunctions.hpp"
#include "Debug.hpp"

template<unsigned DIM>
PoiseuilleImpedanceCalculator<DIM>::PoiseuilleImpedanceCalculator()
{

}

template<unsigned DIM>
PoiseuilleImpedanceCalculator<DIM>::~PoiseuilleImpedanceCalculator()
{

}

template<unsigned DIM>
void PoiseuilleImpedanceCalculator<DIM>::Calculate(boost::shared_ptr<VascularNetwork<DIM> > vascularNetwork)
{
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > segments = vascularNetwork->GetVesselSegments();
    for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)
    {
        double length = segments[segment_index]->GetLength();
        double radius = segments[segment_index]->GetRadius();
        double viscosity = segments[segment_index]->GetFlowProperties()->GetViscosity();
        if (radius <= 0.0)
        {
            EXCEPTION("Radius should be a positive number.");
        }
        if (viscosity <= 0.0)
        {
            EXCEPTION("Viscosity should be a positive number.");
        }
        double impedance = 8.0 * viscosity * length / (M_PI * SmallPow(radius, 4u));
        segments[segment_index]->GetFlowProperties()->SetImpedance(impedance);
    }
}

// Explicit instantiation
template class PoiseuilleImpedanceCalculator<2> ;
template class PoiseuilleImpedanceCalculator<3> ;

