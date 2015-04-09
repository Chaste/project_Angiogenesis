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

#include "VasculatureData.hpp"
#include <stdexcept>

VasculatureData::VasculatureData()
{
}

VasculatureData::~VasculatureData()
{
}

template<typename T> T VasculatureData::GetData(const std::string& variableName)
{
	// Try to return the data in the form of the requested type
	try
	{
		return boost::any_cast<T>(mDataMap.at(variableName));
	}
	catch(const boost::bad_any_cast&)
	{
		EXCEPTION("Invalid type specified for the requested key: " << variableName);
	}
	catch(const std::out_of_range&)
	{
		EXCEPTION("No key: '" << variableName << "' found in property register.");
	}
}

std::vector<std::string> VasculatureData::GetKeys(bool castable_to_double) const
{
	std::vector<std::string> keys;
	std::map<std::string, boost::any>::const_iterator it;
	for(it = mDataMap.begin(); it != mDataMap.end(); it++)
	{
		if(!castable_to_double)
		{
			keys.push_back(it->first);
		}
		else
		{
			if(it->second.type()== typeid(int) || it->second.type()== typeid(unsigned)
					|| it->second.type()== typeid(bool) || it->second.type()== typeid(double))
			{
				keys.push_back(it->first);
			}
		}
	}
	return keys;
}

std::map<std::string, boost::any> VasculatureData::GetMap() const
{
	return mDataMap;
}

bool VasculatureData::HasKey(const std::string& rKey) const
{
	// Check if the key is in the map
	std::map<std::string, boost::any>::const_iterator it = mDataMap.find(rKey);
	if (it == mDataMap.end())
	{
		return false;
	}
	return true;
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

void VasculatureData::SetMap(std::map<std::string, boost::any> map)
{
	mDataMap = map;
}

template<typename T> void VasculatureData::SetData(const std::string& variableName, T value)
{
	mDataMap[variableName] = value;
}

// Explicit instantiation

template bool VasculatureData::GetData<bool>(const std::string& variableName);
template double VasculatureData::GetData<double>(const std::string& variableName);
template unsigned VasculatureData::GetData<unsigned>(const std::string& variableName);
template std::vector<double> VasculatureData::GetData<std::vector<double> >(const std::string& variableName);

template bool VasculatureData::IsType<bool>(const std::string& variableName);
template bool VasculatureData::IsType<double>(const std::string& variableName);
template bool VasculatureData::IsType<unsigned>(const std::string& variableName);
template bool VasculatureData::IsType<std::vector<double> >(const std::string& variableName);

template void VasculatureData::SetData(const std::string& variableName, bool value);
template void VasculatureData::SetData(const std::string& variableName, double value);
template void VasculatureData::SetData(const std::string& variableName, unsigned value);
template void VasculatureData::SetData(const std::string& variableName, std::vector<double> value);
