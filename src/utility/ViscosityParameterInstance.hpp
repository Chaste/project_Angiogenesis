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

#ifndef ViscosityParameterInstance_HPP_
#define ViscosityParameterInstance_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "SmartPointers.hpp"
#include "UnitCollection.hpp"
#include "BaseParameterInstance.hpp"

/**
 * This is a class for storing often used for length type parameters. Note, templating of the
 * individual parameter types is avoided to ease Python wrapping.
 */

class ViscosityParameterInstance : public BaseParameterInstance
{
    /**
     * Archiving
     */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<BaseParameterInstance>(*this);
        ar & mValue;
    }

    /**
     * The value of the parameter
     */
    units::quantity<unit::dynamic_viscosity> mValue;

public:

    /**
     * Constructor
     */
    ViscosityParameterInstance();

    /**
     * Constructor
     */
    ViscosityParameterInstance(units::quantity<unit::dynamic_viscosity> value,
                              const std::string& rName,
                              const std::string& rShortDescription,
                              const std::string& rSymbol,
                              const std::string& rBibliographicInfromation);

    /**
     * Destructor
     */
    virtual ~ViscosityParameterInstance();

    /**
     * Set the default value
     */
    void SetValue(units::quantity<unit::dynamic_viscosity> value);


    units::quantity<unit::dynamic_viscosity> GetValue(const std::string& rCallingClass = "User");

    std::string GetValueAsString();

};

#endif /*ViscosityParameterInstance_HPP_*/
