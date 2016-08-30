/*

Copyright (c) 2005-2016, University of Oxford.
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

#include <cassert>
#include <cmath>

#include "DimensionalSimulationTime.hpp"

/** Pointer to the single instance */
boost::shared_ptr<DimensionalSimulationTime> DimensionalSimulationTime::mpInstance = boost::shared_ptr<DimensionalSimulationTime>();
boost::shared_ptr<SimulationTime> DimensionalSimulationTime::mpSimulationTimeInstance = boost::shared_ptr<SimulationTime>();

boost::shared_ptr<DimensionalSimulationTime> DimensionalSimulationTime::Instance()
{
    if (mpInstance)
    {
        mpInstance = boost::shared_ptr<DimensionalSimulationTime>(new DimensionalSimulationTime);
        std::atexit(Destroy);
    }
    return mpInstance;
}

DimensionalSimulationTime::DimensionalSimulationTime()
    : mReferenceTime(60.0 * unit::seconds)
{
    // Make sure there's only one instance - enforces correct serialization
    assert(bool(mpInstance) == false);

    // Careful, we are grabbing a raw pointer with a shared one. Have to hope no-one uses it to make
    // more shared pointers to SimulationTime and carefully manage deletion in Destroy.
    mpSimulationTimeInstance = boost::shared_ptr<SimulationTime>(SimulationTime::Instance());
}

units::quantity<unit::time> DimensionalSimulationTime::GetReferenceTimeScale()
{
    return mReferenceTime;
}

void DimensionalSimulationTime::SetReferenceTimeScale(units::quantity<unit::time> referenceTimeScale)
{
    mReferenceTime = referenceTimeScale;
}

boost::shared_ptr<SimulationTime> DimensionalSimulationTime::GetSimulationTime()
{
    return mpSimulationTimeInstance;
}

void DimensionalSimulationTime::Destroy()
{
    if(mpSimulationTimeInstance)
    {
        mpSimulationTimeInstance = boost::shared_ptr<SimulationTime>();
        SimulationTime::Destroy();
    }

    if (mpInstance)
    {
        mpInstance = boost::shared_ptr<DimensionalSimulationTime>();
    }
}
