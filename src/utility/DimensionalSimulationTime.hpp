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

#ifndef DIMENSIONALSIMULATIONTIME_HPP_
#define DIMENSIONALSIMULATIONTIME_HPP_

#include <boost/shared_ptr.hpp>
#include "ChasteSerialization.hpp"
#include <boost/serialization/shared_ptr.hpp>
#include "SerializableSingleton.hpp"
#include "TimeStepper.hpp"
#include "UnitCollection.hpp"
#include "SimulationTime.hpp"

/**
 * Simulation time object with extra information regarding the units of the time increment
 */
class DimensionalSimulationTime : public SimulationTime
{

    units::quantity<unit::time> mReferenceTime;

public:

    /**
     * @return a pointer to the simulation time object.
     * The first time this is called the simulation time object is created.
     */
    static DimensionalSimulationTime* Instance();

    /**
     * @return the reference time scale
     */
    units::quantity<unit::time> GetReferenceTimeScale();

    /**
     * Sets reference time scale
     */
    void SetReferenceTimeScale(units::quantity<unit::time> referenceTimeScale);


protected:
    /**
     * Default simulation time constructor
     *
     * Sets up time, you must set the start time,
     * end time and number of time steps before using the object.
     */
    DimensionalSimulationTime();

private:
    /**
     * A pointer to the singleton instance of this class.
     */
    static DimensionalSimulationTime* mpInstance;

    static boost::shared_ptr<TimeStepper> mpTimeStepper;

};

#endif /*SIMULATIONTIME_HPP_*/
