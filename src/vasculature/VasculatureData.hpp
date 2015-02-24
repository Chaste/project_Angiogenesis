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

#ifndef VASCULATUREDATA_HPP_
#define VASCULATUREDATA_HPP_

#include <map>
#include <vector>
#include <string>
#include <boost/any.hpp>
#include <boost/enable_shared_from_this.hpp>
#include "Exception.hpp"
#include "SmartPointers.hpp"

class VasculatureData : public boost::enable_shared_from_this<VasculatureData>
{
    std::map<std::string, boost::any> mDataMap;

public:

    VasculatureData();

    ~VasculatureData();

    std::map<std::string, boost::any> & GetMap();

    template<typename T> T GetData(const std::string& variableName);

    // These classes help with python wrapping as conversion to-from boost::any type is difficult.
    // The generic GetData method is preferred for C++ usage.
    // todo remove when a python dict->vasculature data converter is implemented
    double GetDoubleData(const std::string& variableName);

    unsigned GetUnsignedData(const std::string& variableName);

    std::vector<double> GetVectorDoubleData(const std::string& variableName);

    template<typename T> bool IsType(const std::string& variableName);

    void SetMap(std::map<std::string, boost::any> map);

    void SetData(const std::string& variableName, const boost::any& value);

    // These classes help with python wrapping as conversion to-from boost::any type is difficult.
    // The generic SetData method is preferred for C++ usage.
    // todo remove when a python dict->vasculature data converter is implemented
    void SetDoubleData(const std::string& variableName, double value);

    void SetUnsignedData(const std::string& variableName, unsigned value);

    void SetVectorDoubleData(const std::string& variableName, std::vector<double> value);
};

// Templated methods are defined here as explicit instantiation would limit the types that can
// be stored in the data maps.
template<typename T> T VasculatureData::GetData(const std::string& variableName)
{
	// Check if the key is in the map
	std::map<std::string, boost::any>::const_iterator it = mDataMap.find(variableName);
	if (it == mDataMap.end())
	{
		EXCEPTION("No key: '" << variableName << "' found in property register.");
	}

	// Try to return the data in the form of the requested type
	try
	{
		return boost::any_cast<T>(it->second);
	}
	catch(const boost::bad_any_cast&)
	{
		EXCEPTION("Invalid type specified for the requested key: " << variableName);
	}
}

template<typename T> bool VasculatureData::IsType(const std::string& variableName)
{
	try
	{
		boost::any_cast<T>(mDataMap[variableName]);
		return true;
	}
	catch(const boost::bad_any_cast&)
	{
		return false;
	}
}

#endif /* VASCULATUREDATA_HPP_ */
