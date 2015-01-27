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
#include "CaVascularNetworkNode.hpp"
#include "SmartPointers.hpp"
#include "VascularNetworkData.hpp"
#include "ChastePoint.hpp"
#include "CaVesselSegment.hpp"
#include "FakePetscSetup.hpp"

class TestVesselSegment : public AbstractCellBasedTestSuite
{
public:

	typedef boost::shared_ptr<CaVascularNetworkNode<2> > NodePtr2;
	typedef boost::shared_ptr<CaVascularNetworkNode<3> > NodePtr3;
	typedef boost::shared_ptr<CaVesselSegment<2> > SegmentPtr2;
	typedef boost::shared_ptr<CaVesselSegment<3> > SegmentPtr3;

    void TestConstructor() throw(Exception)
    {
    	typedef boost::shared_ptr<CaVesselSegment<2> > FooPtr;

    	// Make some nodes
    	ChastePoint<2> point1(1.0, 2.0);
    	ChastePoint<2> point2(3.0, 4.0);
    	NodePtr2 pNode1(CaVascularNetworkNode<2>::Create(point1));
    	NodePtr2 pNode2(CaVascularNetworkNode<2>::Create(point2));

    	// Make a segment
    	SegmentPtr2 pSegment(CaVesselSegment<2>::Create(pNode1, pNode2));

    	// Check the locations
    	TS_ASSERT(point1.IsSamePoint(pSegment->GetNodes().first->GetLocation()));
    	TS_ASSERT(point1.IsSamePoint(pSegment->GetNodes(0)->GetLocation()));
    	TS_ASSERT(point2.IsSamePoint(pSegment->GetNodes().second->GetLocation()));
    	TS_ASSERT(point2.IsSamePoint(pSegment->GetNodes(1)->GetLocation()));
    	TS_ASSERT_THROWS_THIS(pSegment->GetNodes(2), "A node index other than 0 or 1 has been requested for a Vessel Segment.");

    	// Test simple Getters and Setters
    	pSegment->SetId(5u);
		std::string label = "Inlet";
		pSegment->SetLabel(label);
		TS_ASSERT_EQUALS(pSegment->GetId(), 5u);
		TS_ASSERT_EQUALS(pSegment->rGetLabel().c_str(), label.c_str());

		// Test replacing a node
    	ChastePoint<2> point3(6.0, 7.0);
    	ChastePoint<2> point4(8.0, 9.0);
    	NodePtr2 pNode3(CaVascularNetworkNode<2>::Create(point3));
    	NodePtr2 pNode4(CaVascularNetworkNode<2>::Create(point4));

    	pSegment->ReplaceNode(0, pNode3);
    	pSegment->ReplaceNode(1, pNode4);
    	TS_ASSERT(point3.IsSamePoint(pSegment->GetNodes().first->GetLocation()));
    	TS_ASSERT(point4.IsSamePoint(pSegment->GetNodes().second->GetLocation()));
    	TS_ASSERT_THROWS_THIS(pSegment->ReplaceNode(2, pNode3), "A node index other than 0 or 1 has been requested for a Vessel Segment.");
    }

	void TestAccessingData() throw(Exception)
	{
    	// Make some nodes
    	ChastePoint<3> point1(1.0, 2.0, 6.0);
    	ChastePoint<3> point2(3.0, 4.0, 7.0);
    	NodePtr3 pNode1(CaVascularNetworkNode<3>::Create(point1));
    	NodePtr3 pNode2(CaVascularNetworkNode<3>::Create(point2));

    	// Make a segment
    	SegmentPtr3 pSegment(CaVesselSegment<3>::Create(pNode1, pNode2));

		// Set some data
		double radius = 5.5;
		pSegment->GetDataContainer()->SetData("radius", radius);
		TS_ASSERT_DELTA(pSegment->GetDataContainer()->GetData<double>("radius"), radius, 1.e-6);

		// Replace the existing data container with a new one
		boost::shared_ptr<VascularNetworkData> pDataContainer(new VascularNetworkData());
		double haematocrit = 7.5;
		pDataContainer->SetData("haematocrit", haematocrit);
		pSegment->SetDataContainer(pDataContainer);
		TS_ASSERT_DELTA(pSegment->GetDataContainer()->GetData<double>("haematocrit"), haematocrit, 1.e-6);
	}

	void TestGeometricFeatures() throw(Exception)
	{
    	//Check the returned length
    	ChastePoint<2> point1(1.0, 2.0);
    	ChastePoint<2> point2(3.0, 4.0);
    	ChastePoint<3> point3(3.0, 4.0, 5.0);
    	ChastePoint<3> point4(6.0, 7.0, 8.0);

    	NodePtr2 pNode1(CaVascularNetworkNode<2>::Create(point1));
    	NodePtr2 pNode2(CaVascularNetworkNode<2>::Create(point2));
    	NodePtr3 pNode3(CaVascularNetworkNode<3>::Create(point3));
    	NodePtr3 pNode4(CaVascularNetworkNode<3>::Create(point4));

    	SegmentPtr2 pSegment(CaVesselSegment<2>::Create(pNode1, pNode2));
    	SegmentPtr3 pSegment2(CaVesselSegment<3>::Create(pNode3, pNode4));

    	TS_ASSERT_DELTA(pSegment->GetLength(), std::sqrt(8.0), 1.e-6);
    	TS_ASSERT_DELTA(pSegment2->GetLength(), std::sqrt(27.0), 1.e-6);
	}
};

#endif /*TESTVASCULARSEGMENT_HPP_*/
