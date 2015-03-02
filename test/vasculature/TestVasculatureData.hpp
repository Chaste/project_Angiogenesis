/*

Copyright (c) 2005-2014, University of Oxford.
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

#ifndef TESTVASCULARNETWORKDATA_HPP_
#define TESTVASCULARNETWORKDATA_HPP_

#include <cxxtest/TestSuite.h>
#include "FakePetscSetup.hpp"
#include "VasculatureData.hpp"

class TestVasculatureData : public CxxTest::TestSuite
{
public:

    /*
     */
    void TestConstructor() throw(Exception)
    {
    	// Make a data instance
    	VasculatureData data;

    	// Pass in some data
    	double value = 12.0;
    	data.SetData("radius", value);

    	// Check that the correct type is distinguished
    	TS_ASSERT(data.IsType<double>("radius"));
    	TS_ASSERT(!data.IsType<std::vector<double> >("radius"));

    	// Check that the data is correctly returned
    	TS_ASSERT_DELTA(data.GetData<double>("radius"), 12.0, 1.e-6);

    	// Check that a suitable Exception is thrown if an incorrect type is
    	// specified.
    	TS_ASSERT_THROWS_THIS(data.GetData<std::vector<double> >("radius"), "Invalid type specified for the requested key: radius");

    	// Check that a suitable Exception is thrown if a non-existing key is requested
    	TS_ASSERT_THROWS_THIS(data.GetData<double>("nonsense"), "No key: 'nonsense' found in property register.");
    }
};

#endif /*TESTVASCULARNETWORKDATA_HPP_*/
