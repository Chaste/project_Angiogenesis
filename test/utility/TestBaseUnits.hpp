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

#ifndef TESTBASEUNITS_HPP
#define TESTBASEUNITS_HPP

#include <cxxtest/TestSuite.h>
#include <SmartPointers.hpp>
#include "UnitCollection.hpp"
#include "BaseUnits.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * Check setting up and destroying the singleton
 */
class TestUnitCollection : public CxxTest::TestSuite
{

public:

    void TestSetupAndDestroy()
    {
        BaseUnits::Instance()->SetReferenceLengthScale(10.0*unit::metres);
        BaseUnits::Instance()->SetReferenceMassScale(15.0*unit::kg);
        BaseUnits::Instance()->SetReferenceTimeScale(20.0*unit::seconds);
        TS_ASSERT_DELTA(BaseUnits::Instance()->GetReferenceLengthScale().value(), 10.0, 1.e-6);
        TS_ASSERT_DELTA(BaseUnits::Instance()->GetReferenceMassScale().value(), 15.0, 1.e-6);
        TS_ASSERT_DELTA(BaseUnits::Instance()->GetReferenceTimeScale().value(), 20.0, 1.e-6);

        BaseUnits::Instance()->Destroy();
        TS_ASSERT_DELTA(BaseUnits::Instance()->GetReferenceLengthScale().value(), 1.e-6, 1.e-8);
        TS_ASSERT_DELTA(BaseUnits::Instance()->GetReferenceMassScale().value(), 1.0, 1.e-6);
        TS_ASSERT_DELTA(BaseUnits::Instance()->GetReferenceTimeScale().value(), 60.0, 1.e-6);
    }

};

#endif // TESTBASEUNITS_HPP
