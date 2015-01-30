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
#include "VascularNode.hpp"
#include "CaVesselSegment.hpp"
#include "CaVessel.hpp"
#include "ChastePoint.hpp"
#include "VasculatureData.hpp"
#include "FakePetscSetup.hpp"

class TestCaVessel : public CxxTest::TestSuite
{
public:

	typedef boost::shared_ptr<VascularNode<2> > NodePtr2;
	typedef boost::shared_ptr<VascularNode<3> > NodePtr3;
	typedef boost::shared_ptr<CaVesselSegment<2> > SegmentPtr2;
	typedef boost::shared_ptr<CaVesselSegment<3> > SegmentPtr3;
	typedef boost::shared_ptr<CaVessel<2> > VesselPtr2;
	typedef boost::shared_ptr<CaVessel<3> > VesselPtr3;

	void TestConstructor() throw(Exception)
	{
    	// Make some nodes
		std::vector<ChastePoint<2> > points;
		points.push_back(ChastePoint<2>(0.0, 0.0));
		points.push_back(ChastePoint<2>(1.0, 2.0));
		points.push_back(ChastePoint<2>(3.0, 4.0));
		points.push_back(ChastePoint<2>(4.0, 5.0));
		points.push_back(ChastePoint<2>(7.0, 8.0));
		points.push_back(ChastePoint<2>(8.0, 9.0));

		std::vector<NodePtr2> nodes;
		for(unsigned i=0; i < points.size(); i++)
		{
			nodes.push_back(NodePtr2 (VascularNode<2>::Create(points[i])));
		}

    	// Make some segments
		SegmentPtr2 pSegment0(CaVesselSegment<2>::Create(nodes[0], nodes[1]));
		SegmentPtr2 pSegment1(CaVesselSegment<2>::Create(nodes[1], nodes[2]));
		SegmentPtr2 pSegment2(CaVesselSegment<2>::Create(nodes[2], nodes[3]));
		SegmentPtr2 pSegment3(CaVesselSegment<2>::Create(nodes[4], nodes[5]));

		// Make a vessel
		VesselPtr2 pVessel1(CaVessel<2>::Create(pSegment1));

		std::vector<SegmentPtr2> good_segments;
		good_segments.push_back(pSegment1);
		good_segments.push_back(pSegment2);
		VesselPtr2 pVessel2(CaVessel<2>::Create(good_segments));

		std::vector<SegmentPtr2> bad_segments = good_segments;
		bad_segments.push_back(pSegment3);
		TS_ASSERT_THROWS_THIS(VesselPtr2 pVessel3(CaVessel<2>::Create(bad_segments));,"Input vessel segments are not attached in the correct order.");

		// Check that locations are correct
		TS_ASSERT(pVessel1->GetStartNode()->IsCoincident(points[1]));
		TS_ASSERT(pVessel2->GetEndNode()->IsCoincident(points[3]));

		// Check that segments are correctly returned
		TS_ASSERT_EQUALS(pVessel2->GetNumberOfSegments(), 2u);
		TS_ASSERT_EQUALS(pVessel2->GetSegments().size(), 2u);
		TS_ASSERT(pVessel2->GetSegments(0)->GetNodes(0)->IsCoincident(points[1]));

		// Test simple Getters and Setters
		pVessel1->SetId(5u);
		std::string label = "Inlet";
		pVessel1->SetLabel(label);
		TS_ASSERT_EQUALS(pVessel1->GetId(), 5u);
		TS_ASSERT_EQUALS(pVessel1->rGetLabel().c_str(), label.c_str());
	}

	void TestAddingAndRemovingSegments() throw(Exception)
	{
    	// Make some nodes
		std::vector<ChastePoint<2> > points;
		points.push_back(ChastePoint<2>(0.0, 0.0));
		points.push_back(ChastePoint<2>(1.0, 2.0));
		points.push_back(ChastePoint<2>(2.0, 3.0));
		points.push_back(ChastePoint<2>(4.0, 5.0));
		points.push_back(ChastePoint<2>(7.0, 8.0));
		points.push_back(ChastePoint<2>(8.0, 9.0));

		std::vector<NodePtr2> nodes;
		for(unsigned i=0; i < points.size(); i++)
		{
			nodes.push_back(NodePtr2 (VascularNode<2>::Create(points[i])));
		}

    	// Make some segments
		SegmentPtr2 pSegment0(CaVesselSegment<2>::Create(nodes[0], nodes[1]));
		SegmentPtr2 pSegment1(CaVesselSegment<2>::Create(nodes[1], nodes[2]));
		SegmentPtr2 pSegment2(CaVesselSegment<2>::Create(nodes[2], nodes[3]));
		SegmentPtr2 pSegment3(CaVesselSegment<2>::Create(nodes[4], nodes[5]));

		// Make a vessel
		VesselPtr2 pVessel1(CaVessel<2>::Create(pSegment1));

		// Try adding a segment to the start and end
		pVessel1->AddSegments(pSegment0);
		pVessel1->AddSegments(pSegment2);
		TS_ASSERT_EQUALS(pVessel1->GetNumberOfSegments(), 3u);

		// Try adding a disconnected segment
		TS_ASSERT_THROWS_THIS(pVessel1->AddSegments(pSegment3),"Input vessel segment does not coincide with any end of the vessel.");

		// Remove the segments from the ends
		pVessel1->RemoveSegments(true, true);
		TS_ASSERT_EQUALS(pVessel1->GetNumberOfSegments(), 1u);

		// Vector version of adding.
		std::vector<SegmentPtr2> good_segments;
		good_segments.push_back(pSegment1);
		good_segments.push_back(pSegment2);

		std::vector<SegmentPtr2> bad_segments = good_segments;
		bad_segments.push_back(pSegment3);
	}

	void TestAccessingData() throw(Exception)
	{
    	// Make some nodes
    	ChastePoint<3> point1(1.0, 2.0, 6.0);
    	ChastePoint<3> point2(3.0, 4.0, 7.0);

    	NodePtr3 pNode1(VascularNode<3>::Create(point1));
    	NodePtr3 pNode2(VascularNode<3>::Create(point2));

    	// Make a segment
    	SegmentPtr3 pSegment1(CaVesselSegment<3>::Create(pNode1, pNode2));

		// Make a vessel
    	VesselPtr3 pVessel(CaVessel<3>::Create(pSegment1));

		// Set some data
		double radius = 5.5;
		pVessel->GetDataContainer()->SetData("radius", radius);
		TS_ASSERT_DELTA(pVessel->GetDataContainer()->GetData<double>("radius"), radius, 1.e-6);

		// Replace the existing data container with a new one
		boost::shared_ptr<VasculatureData> pDataContainer(new VasculatureData());
		double haematocrit = 7.5;
		pDataContainer->SetData("haematocrit", haematocrit);
		pVessel->SetDataContainer(pDataContainer);
		TS_ASSERT_DELTA(pVessel->GetDataContainer()->GetData<double>("haematocrit"), haematocrit, 1.e-6);
	}
};

#endif /*TESTCAVESSEL_HPP_*/
