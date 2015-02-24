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

#ifndef TESTVASCULARSEGMENT_HPP_
#define TESTVASCULARSEGMENT_HPP_

#include <math.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "VascularNode.hpp"
#include "SmartPointers.hpp"
#include "VasculatureData.hpp"
#include "ChastePoint.hpp"
#include "CaVesselSegment.hpp"
#include "CaVessel.hpp"
#include "FakePetscSetup.hpp"

class TestVesselSegment : public AbstractCellBasedTestSuite
{
public:

	typedef boost::shared_ptr<VascularNode<2> > NodePtr2;
	typedef boost::shared_ptr<VascularNode<3> > NodePtr3;
	typedef boost::shared_ptr<CaVesselSegment<2> > SegmentPtr2;
	typedef boost::shared_ptr<CaVesselSegment<3> > SegmentPtr3;

    void TestConstructor() throw(Exception)
    {

    	// Make some nodes
    	ChastePoint<2> point1(1.0, 2.0);
    	ChastePoint<2> point2(3.0, 4.0);
    	NodePtr2 pNode1(VascularNode<2>::Create(point1));
    	NodePtr2 pNode2(VascularNode<2>::Create(point2));

    	// Check for an exception if the segment is defined with the same nodes
    	TS_ASSERT_THROWS_THIS(SegmentPtr2 (CaVesselSegment<2>::Create(pNode1, pNode1)),
    			"Attempted to assign the same node to both ends of a vessel segment.");

    	// Make a segment
    	SegmentPtr2 pSegment(CaVesselSegment<2>::Create(pNode1, pNode2));

    	// Check the locations
    	TS_ASSERT(pSegment->GetNodes().first->IsCoincident(point1));
    	TS_ASSERT(pSegment->GetNode(0)->IsCoincident(point1));
    	TS_ASSERT(pSegment->GetNodes().second->IsCoincident(point2));
    	TS_ASSERT(pSegment->GetNode(1)->IsCoincident(point2));
    	TS_ASSERT_THROWS_THIS(pSegment->GetNode(2),
    			"A node index other than 0 or 1 has been requested for a Vessel Segment.");

    	// Test simple Getters and Setters
    	pSegment->SetId(5u);
		std::string label = "Inlet";
		pSegment->SetLabel(label);
		TS_ASSERT_EQUALS(pSegment->GetId(), 5u);
		TS_ASSERT_EQUALS(pSegment->rGetLabel().c_str(), label.c_str());

		// Test replacing a node
    	ChastePoint<2> point3(6.0, 7.0);
    	ChastePoint<2> point4(8.0, 9.0);
    	NodePtr2 pNode3(VascularNode<2>::Create(point3));
    	NodePtr2 pNode4(VascularNode<2>::Create(point4));

    	pSegment->ReplaceNode(0, pNode3);
    	pSegment->ReplaceNode(1, pNode4);
    	TS_ASSERT(pSegment->GetNodes().first->IsCoincident(point3));
    	TS_ASSERT(pSegment->GetNodes().second->IsCoincident(point4));
    	TS_ASSERT_THROWS_THIS(pSegment->ReplaceNode(2, pNode3), "A node index other than 0 or 1 has been requested for a Vessel Segment.");

    	// Test point distance calculation
    	double distance1 = 0.7071;
    	double distance2 = 2.2360;
    	ChastePoint<2> point5(7.0, 7.0); // projection is inside segment
    	TS_ASSERT_DELTA(pSegment->GetDistance(point5), distance1, 1.e-3);

    	ChastePoint<2> point6(10.0, 10.0); // projection is outside segment
    	TS_ASSERT_DELTA(pSegment->GetDistance(point6), distance2, 1.e-3);

    	ChastePoint<2> point7(5.0, 5.0); // projection is outside segment
    	TS_ASSERT_DELTA(pSegment->GetDistance(point7), distance2, 1.e-3);
    }

	void TestAccessingData() throw(Exception)
	{
    	// Make some nodes
    	ChastePoint<3> point1(1.0, 2.0, 6.0);
    	ChastePoint<3> point2(3.0, 4.0, 7.0);
    	NodePtr3 pNode1(VascularNode<3>::Create(point1));
    	NodePtr3 pNode2(VascularNode<3>::Create(point2));

    	// Make a segment
    	SegmentPtr3 pSegment(CaVesselSegment<3>::Create(pNode1, pNode2));

		// Set some data
		double radius = 5.5;
		pSegment->GetDataContainer().SetData("radius", radius);
		TS_ASSERT_DELTA(pSegment->GetDataContainer().GetData<double>("radius"), radius, 1.e-6);

		// Replace the existing data container with a new one
		VasculatureData dataContainer;
		double haematocrit = 7.5;
		dataContainer.SetData("haematocrit", haematocrit);
		pSegment->SetDataContainer(dataContainer);
		TS_ASSERT_DELTA(pSegment->GetDataContainer().GetData<double>("haematocrit"), haematocrit, 1.e-6);
	}

	void TestGeometricFeatures() throw(Exception)
	{
    	//Check the returned length
    	ChastePoint<2> point1(1.0, 2.0);
    	ChastePoint<2> point2(3.0, 4.0);
    	ChastePoint<3> point3(3.0, 4.0, 5.0);
    	ChastePoint<3> point4(6.0, 7.0, 8.0);

    	NodePtr2 pNode1(VascularNode<2>::Create(point1));
    	NodePtr2 pNode2(VascularNode<2>::Create(point2));
    	NodePtr3 pNode3(VascularNode<3>::Create(point3));
    	NodePtr3 pNode4(VascularNode<3>::Create(point4));

    	SegmentPtr2 pSegment(CaVesselSegment<2>::Create(pNode1, pNode2));
    	SegmentPtr3 pSegment2(CaVesselSegment<3>::Create(pNode3, pNode4));

    	TS_ASSERT_DELTA(pSegment->GetLength(), std::sqrt(8.0), 1.e-6);
    	TS_ASSERT_DELTA(pSegment2->GetLength(), std::sqrt(27.0), 1.e-6);
	}

	void TestAddingAndRemovingVessels() throw(Exception)
	{
		// Make some nodes
		ChastePoint<2> point1(4.0, 3.0);
		ChastePoint<2> point2(4.0, 5.0);
		ChastePoint<2> point3(5.0, 6.0);
		boost::shared_ptr<VascularNode<2> > pNode(new VascularNode<2>(point1));
		boost::shared_ptr<VascularNode<2> > pNode2(new VascularNode<2>(point2));
		boost::shared_ptr<VascularNode<2> > pNode3(new VascularNode<2>(point3));

		// Make some vessel segments
		boost::shared_ptr<CaVesselSegment<2> > pVesselSegment(CaVesselSegment<2>::Create(pNode, pNode2));
		boost::shared_ptr<CaVesselSegment<2> > pVesselSegment2(CaVesselSegment<2>::Create(pNode2, pNode3));

		TS_ASSERT_THROWS_THIS(boost::shared_ptr<CaVessel<2> > vessel = pVesselSegment->GetVessel(),
				"A vessel has been requested but this segment doesn't have one.");

		// Make a vessel and check that it has been suitably added to the segment
		boost::shared_ptr<CaVessel<2> > pVessel(CaVessel<2>::Create(pVesselSegment));
		TS_ASSERT(pVesselSegment->GetNode(0)->IsCoincident(pVesselSegment->GetVessel()->GetSegment(0)->GetNode(0)));

		// Add a different vessel
        boost::shared_ptr<CaVessel<2> > pVessel2(CaVessel<2>::Create(pVesselSegment));
        TS_ASSERT(pVesselSegment->GetNode(0)->IsCoincident(pVesselSegment->GetVessel()->GetSegment(0)->GetNode(0)));

		// Try removing a segment from the vessel
		pVessel->AddSegment(pVesselSegment2);
		pVessel->RemoveSegments(SegmentLocation::Start);
        TS_ASSERT_THROWS_THIS(pVessel->RemoveSegments(SegmentLocation::End),
                "Vessel must have at least one segment.");

		TS_ASSERT_THROWS_THIS(boost::shared_ptr<CaVessel<2> > vessel = pVesselSegment->GetVessel(),
				"A vessel has been requested but this segment doesn't have one.");
	}
};

#endif /*TESTVASCULARSEGMENT_HPP_*/
