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

#ifndef BaseUnits_HPP_
#define BaseUnits_HPP_

#include <boost/shared_ptr.hpp>
#include "ChasteSerialization.hpp"
#include <boost/serialization/shared_ptr.hpp>
#include "SerializableSingleton.hpp"
#include "UnitCollection.hpp"
#include "SimulationTime.hpp"

/**
 * This singleton manages base units of time and length. Its use needs some care.
 *
 * The default values are Length: microns, Time: minutes, Mass: kg
 *
 * After a length unit is set any features (Dimensional Chaste Point, Readers, Writers, Solvers) created
 * after that point will use the length base unit as their reference length scale. Any already
 * existing features won't be changed. It is better to explicitly set reference length scales in a feature
 * rather than changing this value often in a simulation.
 *
 * After a time unit is set any timers depending on SimulationTime will use it as the reference time scale. Extra care should be taken if
 * this value is changed during a simulation looping over solvers, to make sure it is changed back when needed.
 */
class BaseUnits : public SerializableSingleton<BaseUnits>
{
    /**
     * A pointer to the singleton instance of this class.
     */
    static boost::shared_ptr<BaseUnits> mpInstance;

    /**
     * The time unit for increments
     */
    units::quantity<unit::time> mTime;

    /**
     * The length unit
     */
    units::quantity<unit::length> mLength;

    /**
     * The mass unit
     */
    units::quantity<unit::mass> mMass;

public:

    /**
     * @return a pointer to the dimensional simulation time object.
     * The first time this is called the simulation time object is created.
     */
    static boost::shared_ptr<BaseUnits> Instance();

    /**
     * @return the reference time scale
     */
    units::quantity<unit::time> GetReferenceTimeScale();

    /**
     * @return the reference time scale
     */
    units::quantity<unit::length> GetReferenceLengthScale();

    /**
     * @return the reference time scale
     */
    units::quantity<unit::mass> GetReferenceMassScale();

    /**
     * Sets reference time scale
     */
    void SetReferenceTimeScale(units::quantity<unit::time> referenceTimeScale);

    /**
     * Sets reference length scale
     */
    void SetReferenceLengthScale(units::quantity<unit::length> referenceLengthScale);

    /**
     * Sets reference mass scale
     */
    void SetReferenceMassScale(units::quantity<unit::mass> referenceMassScale);

    /**
     * Destroy the current BaseUnits instance
     */
    static void Destroy();

protected:

    /**
     * Default  constructor, force instantiation through factory method.
     *
     */
    BaseUnits();

};

#endif /*BaseUnits_HPP_*/