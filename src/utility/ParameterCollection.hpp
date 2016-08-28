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

#ifndef ParameterCollection_HPP_
#define ParameterCollection_HPP_

#include <boost/shared_ptr.hpp>
#include "ChasteSerialization.hpp"
#include <boost/serialization/shared_ptr.hpp>
#include "SerializableSingleton.hpp"
#include "UnitCollection.hpp"
#include "ParameterInstance.hpp"
#include "BaseParameterInstance.hpp"
#include <map>

/**
 * This singleton holds parameter values used in a simulation. It allows the value of parameters
 * used in a simulation to be dumped to file on completion. New parameters can be added at run-time,
 * which can hold meta-data on parameter sources.
 */
class ParameterCollection : public SerializableSingleton<ParameterCollection>
{
    /**
     * A pointer to the singleton instance of this class.
     */
    static ParameterCollection* mpInstance;

    /**
     * Parameter Collection
     */
    std::map<std::string, boost::shared_ptr<BaseParameterInstance> > mParameters;

public:

    /**
     * @return a pointer to the singleton instance
     * The first time this is called the simulation object is created.
     */
    static ParameterCollection* Instance();

    /**
     * Dump the parameters to file
     */
    void DumpToFile(const std::string& rFilename);

    /**
     * Add a parameter
     */
    void AddParameter(boost::shared_ptr<BaseParameterInstance> pParameter);

    /**
     * Destroy the current ParameterCollection instance.
     *
     * This method *must* be called before program exit, to avoid a memory
     * leak.
     */
    static void Destroy();

    template<class UNIT>
    boost::shared_ptr<ParameterInstance<UNIT> > GetParameter(const std::string& rName);

protected:

    /**
     * Default constructor
     *
     * Sets up an initial parameter collection on construction
     */
    ParameterCollection();

};

#endif /*ParameterCollection_HPP_*/