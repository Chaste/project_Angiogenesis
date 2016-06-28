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

#ifndef AbstractVesselNetworkComponentProperties_HPP_
#define AbstractVesselNetworkComponentProperties_HPP_

#include <string>
#include <map>
#include <boost/enable_shared_from_this.hpp>
#include "UnitCollection.hpp"
/**
 * This class contains common functionality for property containers for all vessel network components.
 *
 * Note: It is named 'Abstract' to discourage instantiation, but is not strictly an abstract class.
 * A pure virtual destructor is avoided as it prevents Python wrapping.
 */
template<unsigned DIM>
class AbstractVesselNetworkComponentProperties : public boost::enable_shared_from_this<AbstractVesselNetworkComponentProperties<DIM> >
{

protected:

    /**
     * Reference length scale for dimensionalizing units
     */
    units::quantity<unit::length> mReferenceLength;

    /**
     * Reference time scale for dimensionalizing units
     */
    units::quantity<unit::time> mReferenceTime;

    /**
     * Reference mass scale for dimensionalizing units
     */
    units::quantity<unit::mass> mReferenceMass;

public:

    /**
     * Constructor
     */
    AbstractVesselNetworkComponentProperties();

    /**
     * Destructor
     */
    virtual ~AbstractVesselNetworkComponentProperties();

    /**
     * Return the reference length for the component
     *
     * @return the reference length for the component
     */
    virtual units::quantity<unit::length> GetReferenceLength() const;

    /**
     * Return the reference time for the component
     *
     * @return the reference time for the component
     */
    virtual units::quantity<unit::time> GetReferenceTime() const;

    /**
     * Return the reference mass for the component
     *
     * @return the reference mass for the component
     */
    virtual units::quantity<unit::mass> GetReferenceMass() const;

    /**
     * Return a map of output data for writing to file
     *
     * @return a map of output data for use by writers
     */
    virtual std::map<std::string, double> GetOutputData() const;

    /**
     * Return the reference length for the component in SI units
     * @return the reference length in SI units
     */
    virtual double GetReferenceLengthSI() const;

    /**
     * Return the reference time for the component in SI units
     * @return a pair containing the reference time value and unit as a string.
     */
    virtual double GetReferenceTimeSI() const;

    /**
     * Return the reference mass for the component in SI units
     * @return a pair containing the reference mass value and unit as a string.
     */
    virtual double GetReferenceMassSI() const;

    /**
     * Set the reference length for the component
     *
     * @param referenceLength the reference length
     */
    virtual void SetReferenceLength(units::quantity<unit::length> referenceLength);

    /**
     * Set the reference length for the component in SI units. This is used by the Python interface.
     * @param referenceLength the reference length
     */
    virtual void SetReferenceLengthSI(double referenceLength);

    /**
     * Set the reference time for the component
     *
     * @param referenceTimethe reference time
     */
    virtual void SetReferenceTime(units::quantity<unit::time> referenceTime);

    /**
     * Set the reference time for the component in SI Units. This is used by the Python interface.
     * @param referenceTime the reference time
     */
    virtual void SetReferenceTimeSI(double referenceTime);

    /**
     * Set the reference mass for the component
     *
     * @param referenceMass the reference mass
     */
    virtual void SetReferenceMass(units::quantity<unit::mass> referenceMass);

    /**
     * Set the reference mass for the component in SI Units. This is used by the Python interface.
     * @param referenceMass the reference mass
     */
    virtual void SetReferenceMassSI(double referenceMass);

};

#endif /* AbstractVesselNetworkComponentProperties_HPP_ */
