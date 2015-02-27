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
#include "SmartVasculaturePointers.hpp"
#include "CaVessel.hpp"
#include "ChastePoint.hpp"
#include "VasculatureData.hpp"
#include "FakePetscSetup.hpp"

class TestCaVessel : public CxxTest::TestSuite
{
public:

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

		std::vector<boost::shared_ptr<VascularNode<2> > > nodes;
		for(unsigned i=0; i < points.size(); i++)
		{
			MAKE_VN_PTR_ARGS(VascularNode<2>, pNode, (points[i]));
			nodes.push_back(pNode);
		}

    	// Make some segments
		MAKE_VN_PTR_ARGS(CaVesselSegment<2>, pSegment0, (nodes[0], nodes[1]));
		MAKE_VN_PTR_ARGS(CaVesselSegment<2>, pSegment1, (nodes[1], nodes[2]));
		MAKE_VN_PTR_ARGS(CaVesselSegment<2>, pSegment2, (nodes[2], nodes[3]));
		MAKE_VN_PTR_ARGS(CaVesselSegment<2>, pSegment3, (nodes[4], nodes[5]));

		// Make some vessels
		MAKE_VN_PTR_ARGS(CaVessel<2>, pVessel1, (pSegment1));

		std::vector<boost::shared_ptr<CaVesselSegment<2> > > good_segments;
		good_segments.push_back(pSegment1);
		good_segments.push_back(pSegment2);

		MAKE_VN_PTR_ARGS(CaVessel<2>, pVessel2, (good_segments));
		std::vector<boost::shared_ptr<CaVesselSegment<2> > > bad_segments = good_segments;
		bad_segments.push_back(pSegment3);
		TS_ASSERT_THROWS_THIS(MAKE_VN_PTR_ARGS(CaVessel<2>, pVessel3, (bad_segments)),"Input vessel segments are not attached in the correct order.");
		MAKE_VN_PTR_ARGS(CaVessel<2>, pVessel3, (nodes));

		// Check that locations are correct
		TS_ASSERT(pVessel1->GetStartNode()->IsCoincident(points[1]));
		TS_ASSERT(pVessel2->GetEndNode()->IsCoincident(points[3]));

		// Check that segments are correctly returned
		TS_ASSERT_EQUALS(pVessel2->GetNumberOfSegments(), 2u);
		TS_ASSERT_EQUALS(pVessel3->GetNumberOfSegments(), 5u);
		TS_ASSERT_EQUALS(pVessel2->GetSegments().size(), 2u);
		TS_ASSERT(pVessel2->GetSegment(0)->GetNode(0)->IsCoincident(points[1]));

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

		std::vector<boost::shared_ptr<VascularNode<2> > > nodes;
		for(unsigned i=0; i < points.size(); i++)
		{
			MAKE_VN_PTR_ARGS(VascularNode<2>, pNode, (points[i]));
			nodes.push_back(pNode);
		}

    	// Make some segments
		MAKE_VN_PTR_ARGS(CaVesselSegment<2>, pSegment0, (nodes[0], nodes[1]));
		MAKE_VN_PTR_ARGS(CaVesselSegment<2>, pSegment1, (nodes[1], nodes[2]));
		MAKE_VN_PTR_ARGS(CaVesselSegment<2>, pSegment2, (nodes[2], nodes[3]));
		MAKE_VN_PTR_ARGS(CaVesselSegment<2>, pSegment3, (nodes[4], nodes[5]));

		// Make a vessel
		MAKE_VN_PTR_ARGS(CaVessel<2>, pVessel1, (pSegment1));

		// Try adding a segment to the start and end
		pVessel1->AddSegment(pSegment0);
		pVessel1->AddSegment(pSegment2);
		TS_ASSERT_EQUALS(pVessel1->GetNumberOfSegments(), 3u);

		// Try adding a disconnected segment
		TS_ASSERT_THROWS_THIS(pVessel1->AddSegment(pSegment3),"Input vessel segment does not coincide with any end of the vessel.");

		// Remove the segments from the ends
		pVessel1->RemoveSegments(SegmentLocation::Start);
		pVessel1->RemoveSegments(SegmentLocation::End);
		TS_ASSERT_EQUALS(pVessel1->GetNumberOfSegments(), 1u);

		// Vector version of adding segments
		std::vector<boost::shared_ptr<CaVesselSegment<2> > > good_segments;
		good_segments.push_back(pSegment1);
		good_segments.push_back(pSegment2);

		std::vector<boost::shared_ptr<CaVesselSegment<2> > > bad_segments = good_segments;
		bad_segments.push_back(pSegment3);

		MAKE_VN_PTR_ARGS(CaVessel<2>, pVessel2, (pSegment0));
		pVessel2->AddSegments(good_segments);
		TS_ASSERT_EQUALS(pVessel2->GetNumberOfSegments(), 3u);

		MAKE_VN_PTR_ARGS(CaVessel<2>, pVessel3, (pSegment0));
		TS_ASSERT_THROWS_THIS(pVessel3->AddSegments(bad_segments),"Input vessel segments are not attached in the correct order.");
	}

	void TestAccessingData() throw(Exception)
	{
    	// Make some nodes
    	ChastePoint<3> point1(1.0, 2.0, 6.0);
    	ChastePoint<3> point2(3.0, 4.0, 7.0);

    	MAKE_VN_PTR_ARGS(VascularNode<3>, pNode1, (point1));
    	MAKE_VN_PTR_ARGS(VascularNode<3>, pNode2, (point2));

    	// Make a segment
    	MAKE_VN_PTR_ARGS(CaVesselSegment<3>, pSegment1, (pNode1, pNode2));

		// Make a vessel
    	MAKE_VN_PTR_ARGS(CaVessel<3>, pVessel1, (pSegment1));

        // Set some data
        double radius = 5.5;
        std::string key ="radius";
        pVessel1->SetData(key, radius);

        // Check the key is set
        TS_ASSERT(pVessel1->HasDataKey("radius"));
        TS_ASSERT_EQUALS(pVessel1->GetDataKeys()[0].c_str(), key.c_str());

        bool value_is_castable_to_double = true;
        TS_ASSERT_EQUALS(pVessel1->GetDataKeys(value_is_castable_to_double)[0].c_str(), key.c_str());

        // Check the key value is retrieved
        TS_ASSERT_DELTA(pVessel1->GetData<double>("radius"), radius, 1.e-6);
        TS_ASSERT_DELTA(pVessel1->rGetDataContainer().GetData<double>("radius"), radius, 1.e-6);

        // Replace the existing data container with a new one
        VasculatureData data_container;
        double haematocrit = 7.5;
        data_container.SetData("haematocrit", haematocrit);
        pVessel1->SetDataContainer(data_container);
        TS_ASSERT_DELTA(pVessel1->GetData<double>("haematocrit"), haematocrit, 1.e-6);
	}
};

#endif /*TESTCAVESSEL_HPP_*/
