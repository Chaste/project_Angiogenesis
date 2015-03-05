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

#ifndef TestPoiseuilleImpedanceCalculator_HPP_
#define TestPoiseuilleImpedanceCalculator_HPP_

#include <cxxtest/TestSuite.h>
#include "SmartPointers.hpp"
#include "PoiseuilleImpedanceCalculator.hpp"
#include "VasculatureData.hpp"
#include "MathsCustomFunctions.hpp"

#include "FakePetscSetup.hpp"

class TestPoiseuilleImpedanceCalculator : public CxxTest::TestSuite
{

	typedef boost::shared_ptr<VascularNode<2> > NodePtr2;
	typedef boost::shared_ptr<VascularNode<3> > NodePtr3;
	typedef boost::shared_ptr<CaVesselSegment<2> > SegmentPtr2;
	typedef boost::shared_ptr<CaVesselSegment<3> > SegmentPtr3;
	typedef boost::shared_ptr<CaVessel<2> > VesselPtr2;
	typedef boost::shared_ptr<CaVessel<3> > VesselPtr3;

public:



	void TestFlowThroughSingleSegment() throw(Exception)
	{

		// Make some nodes
		std::vector<ChastePoint<3> > points;
		points.push_back(ChastePoint<3>(0, 0, 0));
		points.push_back(ChastePoint<3>(5, 0, 0));

		std::vector<NodePtr3> nodes;
		for(unsigned i=0; i < points.size(); i++)
		{
			nodes.push_back(NodePtr3 (VascularNode<3>::Create(points[i])));
		}

		SegmentPtr3 p_segment(CaVesselSegment<3>::Create(nodes[0], nodes[1]));
		VesselPtr3 p_vessel(CaVessel<3>::Create(p_segment));

		// Generate the network
		boost::shared_ptr<CaVascularNetwork<3> > p_vascular_network(new CaVascularNetwork<3>());

		p_vascular_network->AddVessel(p_vessel);

		VasculatureData data;
		double viscosity = 2e-3;
		double radius = 5e-6;
		data.SetData("Radius", radius);
		data.SetData("Viscosity", viscosity);
		p_vascular_network->SetVesselData(data);
		p_vascular_network->SetSegmentData(data);

		PoiseuilleImpedanceCalculator<3> calculator;

		calculator.Calculate(p_vascular_network);

		double expected_impedance = 8*viscosity*5/(M_PI*SmallPow(radius,4u));

		TS_ASSERT_DELTA(p_vessel->GetData<double>("Impedance"),expected_impedance,1e-6);
		TS_ASSERT_DELTA(p_segment->GetData<double>("Impedance"),expected_impedance,1e-6);

		p_segment->SetData("Radius",0.0);

		TS_ASSERT_THROWS_THIS(calculator.Calculate(p_vascular_network),"Radius should be a positive number.");

		p_segment->SetData("Radius",5e-6);
		p_segment->SetData("Viscosity",0.0);

		TS_ASSERT_THROWS_THIS(calculator.Calculate(p_vascular_network),"Viscosity should be a positive number.");

	}

};

#endif /*TestPoiseuilleImpedanceCalculator_HPP_*/
