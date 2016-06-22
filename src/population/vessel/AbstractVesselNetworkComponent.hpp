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

#ifndef AbstractVesselNetworkComponent_HPP_
#define AbstractVesselNetworkComponent_HPP_

#include <vector>
#include <string>
#include <map>
#include <boost/enable_shared_from_this.hpp>
#include "UblasVectorInclude.hpp"
#include "UnitCollection.hpp"
#include "SmartPointers.hpp"
#include "AbstractVesselNetworkComponentProperties.hpp"

/**
 * This class contains functionality common to all components of a vessel network
 */
template<unsigned DIM>
class AbstractVesselNetworkComponent
{

protected:

    /**
     * Container for generic component data.
     */
    std::map<std::string, double> mOutputData;

    /**
     * Id tag, useful for post-processing
     */
    unsigned mId;

    /**
     * The component radius
     */
    units::quantity<unit::length> mRadius;

    /**
     * Reference length scale for non-dimensionalizing units
     */
    units::quantity<unit::length> mReferenceLength;

    /**
     * Reference time scale for non-dimensionalizing units
     */
    units::quantity<unit::time> mReferenceTime;

    /**
     * Reference mass scale for non-dimensionalizing units
     */
    units::quantity<unit::mass> mReferenceMass;

public:

    /**
     * Constructor.
     */
    AbstractVesselNetworkComponent();

    /**
     * Destructor
     */
    virtual ~AbstractVesselNetworkComponent() = 0;

    /**
     * Return the component Id
     *
     * @return the component id
     */
    virtual unsigned GetId() const;

    /**
     * Return the output data for the given key. This method is relatively slow compared to GetOutputData as the
     * data map is reconstructed each time. It is used by the Python framework.
     *
     * @param rKey the key to be queried
     * @return the component data for the input key
     */
    virtual double GetOutputDataValue(const std::string& rKey);

    /**
     * Return a map of output data for writers
     * @return a map of component data for use by the vtk writer
     */
    virtual std::map<std::string, double> GetOutputData();

    /**
     * Return the keys of the output data map
     * @param verbose include all flow data
     * @return a map of component data for use by the vtk writer
     */
    virtual std::vector<std::string> GetOutputDataKeys();

    /**
     * Return the dimensional radius of the component
     *
     * @return the dimensional radius of the component
     */
    virtual units::quantity<unit::length> GetDimensionalRadius() const;

    /**
     * Return the non-dimensional radius of the component
     *
     * @return the non-dimensional radius of the component
     */
    virtual double GetRadius() const;

    /**
     * Return the radius of the component in SI units
     *
     * @return the dimensional radius of the component
     */
    virtual double GetRadiusSI() const;
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
     * Assign the Id
     * @param id the id for the component
     */
    virtual void SetId(unsigned id);

    /**
     * Add output data to the component using the identifying key
     * @param rKey the key for the data being assigned to the node
     * @param value the value to be stored
     */
    virtual void SetOutputData(const std::string& rKey, double value);

    /**
     * Set the dimensionless component radius
     * @param radius the component radius
     */
    virtual void SetRadius(double radius);

    /**
     * Set the component radius
     * @param radius the component radius
     */
    virtual void SetDimensionalRadius(units::quantity<unit::length> radius);

    /**
     * Set the component radius
     * @param radius the component radius
     */
    virtual void SetRadiusSI(double radius);

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
     * @param referenceMass the reference mass
     */
    virtual void SetReferenceMass(units::quantity<unit::mass> referenceMass);

    /**
     * Set the reference mass for the component in SI Units. This is used by the Python interface.
     * @param referenceMass the reference mass
     */
    virtual void SetReferenceMassSI(double referenceMass);

};

#endif /* AbstractVesselNetworkComponent_HPP_ */
