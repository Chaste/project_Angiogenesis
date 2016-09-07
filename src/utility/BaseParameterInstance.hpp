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

#ifndef BaseParameterInstance_HPP_
#define BaseParameterInstance_HPP_

#include <ostream>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include "ChasteSerialization.hpp"
#include "SmartPointers.hpp"
#include "UnitCollection.hpp"

/**
 * This is a class for storing often used parameters and their metadata, such as a short textual description, bib information
 * from a paper and units. It is templated over unit type. Specific parameter values can inherit from this class, with
 * hard-coded metadata.
 */
class BaseParameterInstance : public boost::enable_shared_from_this<BaseParameterInstance>
{
    /**
     * Archiving
     */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & mName;
        ar & mShortDescription;
        ar & mSourceInformation;
    }

    /**
     * A name to distinguish it in a ParameterCollection
     */
    std::string mName;

    /**
     * A short string based description
     */
    std::string mShortDescription;

    /**
     * Information on the parameter source, such as bibliographic information
     */
    std::string mSourceInformation;

    std::string mSymbol;

    units::quantity<unit::dimensionless> mBaseValue;

public:

    /**
     * Constructor
     */
    BaseParameterInstance();

    /**
     * Constructor
     */
    BaseParameterInstance(const std::string& rName,
                          const std::string& rShortDescription,
                          const std::string& rSymbol,
                          const std::string& rBibliographicInfromation);

    /**
     * Destructor
     */
    virtual ~BaseParameterInstance();

    /**
     * Return the parameter name
     * @return the parameter name as a string
     */
    std::string GetName();

    /**
     * Return the bibliographic information
     * @return the bibliographic information
     */
    std::string GetBibliographicInformation();

    /**
     * Return the short description
     * @return the short description
     */
    std::string GetShortDescription();

    virtual std::string GetValueAsString();

    std::string GetSymbol();

    void RegisterWithCollection(const std::string& rCallingClass);

    /**
     * Set the bibliographic information
     * @param rSourceInformation the bibliographic information
     */
    void SetBibliographicInformation(const std::string& rSourceInformation);

    /**
     * Set the parameter name, this is used to distinguish it in the ParameterCollection map
     * @param rShortDescription the short description
     */
    void SetName(const std::string& rName);

    /**
     * Set the short description of the parameter
     * @param rShortDescription the short description
     */
    void SetShortDescription(const std::string& rShortDescription);

    void SetSymbol(const std::string& rSymbol);

    friend std::ostream& operator<< (std::ostream& stream, const boost::shared_ptr<BaseParameterInstance>& rParameter);

};

#endif /*BaseParameterInstance_HPP_*/
