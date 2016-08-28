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
#include <iostream>

#include "ParameterCollection.hpp"

/** Pointer to the single instance */
ParameterCollection* ParameterCollection::mpInstance = NULL;

ParameterCollection* ParameterCollection::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new ParameterCollection;
        std::atexit(Destroy);
    }
    return mpInstance;
}

ParameterCollection::ParameterCollection()
    : mParameters()
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == NULL);

}

void ParameterCollection::DumpToFile(const std::string& rFilename)
{

    typedef std::map<std::string, boost::shared_ptr<BaseParameterInstance> >::iterator it_type;
    for(it_type iterator = mParameters.begin(); iterator != mParameters.end(); iterator++)
    {
        std::cout << "Name: " << iterator->first << "Description:" << (iterator->second)->GetName();
    }
}

template<class UNIT>
boost::shared_ptr<ParameterInstance<UNIT> > ParameterCollection::GetParameter(const std::string& rName)
{
    boost::shared_ptr<BaseParameterInstance> p_base_parameter = mParameters[rName];

    // Try to cast to the templated type
    boost::shared_ptr<ParameterInstance<UNIT> > p_parameter = boost::dynamic_pointer_cast<ParameterInstance<UNIT> >(p_base_parameter);

    if(p_parameter)
    {
        return p_parameter;
    }
    else
    {
        EXCEPTION("Named parameter does not match (can not be cast to) template type");
    }

}

void ParameterCollection::AddParameter(boost::shared_ptr<BaseParameterInstance> pParameter)
{
    mParameters[pParameter->GetName()] = pParameter;
}

void ParameterCollection::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = NULL;
    }
}

template boost::shared_ptr<ParameterInstance<unit::length> > ParameterCollection::GetParameter(const std::string&);
template boost::shared_ptr<ParameterInstance<unit::time> > ParameterCollection::GetParameter(const std::string&);
template boost::shared_ptr<ParameterInstance<unit::dimensionless> > ParameterCollection::GetParameter(const std::string&);
template boost::shared_ptr<ParameterInstance<unit::mass> > ParameterCollection::GetParameter(const std::string&);
template boost::shared_ptr<ParameterInstance<unit::pressure> > ParameterCollection::GetParameter(const std::string&);
template boost::shared_ptr<ParameterInstance<unit::rate> > ParameterCollection::GetParameter(const std::string&);
template boost::shared_ptr<ParameterInstance<unit::flow_impedance> > ParameterCollection::GetParameter(const std::string&);
template boost::shared_ptr<ParameterInstance<unit::flow_rate> > ParameterCollection::GetParameter(const std::string&);
