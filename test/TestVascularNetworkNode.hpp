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

#ifndef TESTVASCULARNETWORKNODE_HPP_
#define TESTVASCULARNETWORKNODE_HPP_

#include <cxxtest/TestSuite.h>
#include "SmartPointers.hpp"
#include "CaVascularNetworkNode.hpp"
#include "CaVessel.hpp"
#include "ChastePoint.hpp"
#include "FakePetscSetup.hpp"

class TestVascularNetworkGenerator : public CxxTest::TestSuite
{
public:

	void TestConstructor() throw(Exception)
	{
		// Make a node
		CaVascularNetworkNode<2> node;

		// Test Getters and Setters
		ChastePoint<2> point(1.0, 2.0);
		node.SetLocation(point);

		TS_ASSERT_DELTA(node.GetLocation()[0], point[0], 1.e-6);
		TS_ASSERT_DELTA(node.GetLocation()[1], point[1], 1.e-6);

	}

	void TestAddingAndRemovingVessels() throw(Exception)
	{
		// Make some nodes
		boost::shared_ptr<CaVascularNetworkNode<2> > pNode(new CaVascularNetworkNode<2>);
		boost::shared_ptr<CaVascularNetworkNode<2> > pNode2(new CaVascularNetworkNode<2>);

		// Make some vessels
		boost::shared_ptr<CaVessel<2> > pVessel(new CaVessel<2>);
		boost::shared_ptr<CaVessel<2> > pVessel2(new CaVessel<2>);

		// Add a vessel to the node and remove it again
		TS_ASSERT_EQUALS(pNode->GetNumberOfAdjoiningVessels(), 0u);
		pNode->AddAdjoiningVessel(pVessel);
		TS_ASSERT_EQUALS(pNode->GetNumberOfAdjoiningVessels(), 1u);

        TS_ASSERT_THROWS_THIS(pNode->IsAttachedToVessel(pVessel),
                "The vessel has been added to the node, but the node has not been added to the vessel.");

		pVessel->SetNode1(pNode);
		TS_ASSERT(pNode->IsAttachedToVessel(pVessel));
		TS_ASSERT_EQUALS(pNode->GetAdjoiningVessel(0u), pVessel);

		pNode->RemoveAdjoiningVessel(pVessel);
		TS_ASSERT(!pNode->IsAttachedToVessel(pVessel));
        TS_ASSERT_THROWS_THIS(pNode->GetAdjoiningVessel(1u),
                "Attempted to access a vessel with an out of range index.");

		// Check for a suitable Exception if the vessel is not
		// attached to the node.
        TS_ASSERT_THROWS_THIS(pNode->RemoveAdjoiningVessel(pVessel),
                "Attempted to remove a vessel from a node it is not attached to.");

        // Attempt to attach the same vessel to a node multiple times
        pNode->AddAdjoiningVessel(pVessel);
        TS_ASSERT_THROWS_THIS(pNode->AddAdjoiningVessel(pVessel),
                "Vessels and nodes in inconsistent state.");

        pVessel->SetNode2(pNode);
        pNode->AddAdjoiningVessel(pVessel);
        TS_ASSERT_THROWS_THIS(pNode->AddAdjoiningVessel(pVessel),
                "Vessel is already attached to node twice (at both ends). Cannot attach vessel to same node again.");
	}


	void TestGettingAndSettingProperties()
	{

		boost::shared_ptr<CaVascularNetworkNode<2> > pNode(new CaVascularNetworkNode<2>);

		TS_ASSERT_THROWS_THIS(pNode->GetDoubleDataValue("pressure"),"No double valued property, 'pressure', in property register.");

		TS_ASSERT_THROWS_THIS(pNode->GetBooleanData("isActivelyMigrating"),"No boolean valued property, 'isActivelyMigrating', in property register.");

		pNode->SetDoubleData("pressure", 1e5, "Pa");

		pNode->SetBooleanData("isActivelyMigrating", true);

		TS_ASSERT_EQUALS(pNode->GetDoubleDataValue("pressure"),1e5);
		TS_ASSERT_EQUALS(pNode->GetDoubleDataUnits("pressure"),"Pa");

		TS_ASSERT_EQUALS(pNode->GetBooleanData("isActivelyMigrating"),true);

	}

};

#endif /*TESTVASCULARNETWORKNODE_HPP_*/
