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

#ifndef TESTCAVESSEL_HPP_
#define TESTCAVESSEL__HPP_

#include <cxxtest/TestSuite.h>
#include "SmartPointers.hpp"
#include "CaVascularNetworkNode.hpp"
#include "CaVessel.hpp"
#include "ChastePoint.hpp"
#include "FakePetscSetup.hpp"

class TestCaVessel : public CxxTest::TestSuite
{
public:

	void TestConstructor() throw(Exception)
	{
		// Make a vessel
		CaVessel<2> vessel;

		// Test Getters and Setters
		vessel.SetRadius(10.0);
		TS_ASSERT_DELTA(vessel.GetRadius(), 10.0, 1.e-6);
		vessel.SetPreviousRadius(15.0);
		TS_ASSERT_DELTA(vessel.GetPreviousRadius(), 15.0, 1.e-6);
		vessel.SetHaematocritLevel(0.3);
		TS_ASSERT_DELTA(vessel.GetHaematocritLevel(), 0.3, 1.e-6);
		vessel.SetFlowVelocity(100.0);
		TS_ASSERT_DELTA(vessel.GetFlowVelocity(), 100.0, 1.e-6);
		vessel.SetFlowRate(0.001);
		TS_ASSERT_DELTA(vessel.GetFlowRate(), 0.001, 1.e-6);
		vessel.SetImpedance(300.0);
		TS_ASSERT_DELTA(vessel.GetImpedance(), 300.0, 1.e-6);
		vessel.SetWallShearStress(100.0);
		TS_ASSERT_DELTA(vessel.GetWallShearStress(), 100.0, 1.e-6);
		vessel.SetViscosity(0.001);
		TS_ASSERT_DELTA(vessel.GetViscosity(), 0.001, 1.e-6);
		vessel.SetMechanicalStimulus(90.0);
		TS_ASSERT_DELTA(vessel.GetMechanicalStimulus(), 90.0, 1.e-6);
		vessel.SetMetabolicStimulus(20.0);
		TS_ASSERT_DELTA(vessel.GetMetabolicStimulus(), 20.0, 1.e-6);
		vessel.SetShrinkingStimulus(30.0);
		TS_ASSERT_DELTA(vessel.GetShrinkingStimulus(), 30.0, 1.e-6);
		vessel.SetDownstreamConductedStimulus(50.0);
		TS_ASSERT_DELTA(vessel.GetDownstreamConductedStimulus(), 50.0, 1.e-6);
		vessel.SetUpstreamConductedStimulus(90.0);
		TS_ASSERT_DELTA(vessel.GetUpstreamConductedStimulus(), 90.0, 1.e-6);
	}
};

#endif /*TESTCAVESSEL_HPP_*/
