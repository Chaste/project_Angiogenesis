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

#ifndef TESTVESSELNETWORK_HPP_
#define TESTVESSELNETWORK_HPP_

#include <math.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "VascularNode.hpp"
#include "SmartPointers.hpp"
#include "VasculatureData.hpp"
#include "ChastePoint.hpp"
#include "CaVesselSegment.hpp"
#include "CaVascularNetwork.hpp"
#include "OutputFileHandler.hpp"
#include "UblasIncludes.hpp"
#include "FakePetscSetup.hpp"
#include "VasculatureGenerator.hpp"

class TestVesselNetwork : public AbstractCellBasedTestSuite
{
public:

    void TestAddingVessels() throw(Exception)
    {
        // Make some nodes
        std::vector<boost::shared_ptr<VascularNode<3> > > nodes;
        for(unsigned idx=0; idx < 4; idx++)
        {
            nodes.push_back(VascularNode<3>::Create(double(idx), 0.0, 0.0));
        }

        // Make some vessels
        std::vector<boost::shared_ptr<CaVessel<3> > > vessels;
        for(unsigned idx=0; idx < 1; idx++)
        {
            vessels.push_back(CaVessel<3>::Create(nodes[idx], nodes[idx+1]));
        }
        boost::shared_ptr<CaVessel<3> > p_end_vessel = CaVessel<3>::Create(nodes[2], nodes[3]);

        // Make a network
        CaVascularNetwork<3> vessel_network;
        vessel_network.AddVessels(vessels);
        vessel_network.AddVessel(p_end_vessel);

        TS_ASSERT_EQUALS(vessel_network.GetNodes().size(), 4u);
    }

    void TestSettingNetworkData() throw(Exception)
    {
        // Make some nodes
        std::vector<boost::shared_ptr<VascularNode<3> > > nodes;
        for(unsigned idx=0; idx < 4; idx++)
        {
            nodes.push_back(VascularNode<3>::Create(double(idx), 0.0, 0.0));
        }

        // Make some vessels
        std::vector<boost::shared_ptr<CaVessel<3> > > vessels;
        for(unsigned idx=0; idx < 2; idx++)
        {
            vessels.push_back(CaVessel<3>::Create(nodes[idx], nodes[idx+1]));
        }

        // Make a network
        CaVascularNetwork<3> vessel_network;
        vessel_network.AddVessels(vessels);

        // Make some network data
        VasculatureData data;
        double haematocrit = 0.4;
        bool has_flow = true;
        unsigned some_index = 5u;
        data.SetData("Haematocrit", haematocrit);
        data.SetData("Has Flow", has_flow);
        data.SetData("SomeIndex", some_index);
        vessel_network.SetVesselData(data);

        vessel_network.SetNodeRadii(10.0);
        vessel_network.SetSegmentRadii(12.0);
    }

    void TestCopyingAndMovingNetwork() throw(Exception)
    {
        // Make some nodes
        std::vector<boost::shared_ptr<VascularNode<3> > > nodes;
        for(unsigned idx=0; idx < 4; idx++)
        {
            nodes.push_back(VascularNode<3>::Create(double(idx), 0.0, 0.0));
        }

        // Make some vessels
        std::vector<boost::shared_ptr<CaVessel<3> > > vessels;
        for(unsigned idx=0; idx < 3; idx++)
        {
            vessels.push_back(CaVessel<3>::Create(CaVesselSegment<3>::Create(nodes[idx], nodes[idx+1])));
        }

        // Make a network
        CaVascularNetwork<3> vessel_network;
        vessel_network.AddVessels(vessels);
        vessel_network.MergeCoincidentNodes();

        // Move the network
        c_vector<double, 3> translation_vector = zero_vector<double>(3);
        translation_vector[1] = 2.0;
        vessel_network.Translate(translation_vector);
        TS_ASSERT_DELTA(vessel_network.GetVessel(0)->GetSegment(0)->GetNode(0)->GetLocation()[0], 0.0, 1.e-6);
        TS_ASSERT_DELTA(vessel_network.GetVessel(0)->GetSegment(0)->GetNode(0)->GetLocation()[1], 2.0, 1.e-6);
        TS_ASSERT_DELTA(vessel_network.GetVessel(0)->GetSegment(0)->GetNode(1)->GetLocation()[0], 1.0, 1.e-6);
        TS_ASSERT_DELTA(vessel_network.GetVessel(0)->GetSegment(0)->GetNode(1)->GetLocation()[1], 2.0, 1.e-6);
        TS_ASSERT_DELTA(vessel_network.GetVessel(1)->GetSegment(0)->GetNode(1)->GetLocation()[0], 2.0, 1.e-6);
        TS_ASSERT_DELTA(vessel_network.GetVessel(1)->GetSegment(0)->GetNode(1)->GetLocation()[1], 2.0, 1.e-6);

        // Copy the network
        std::vector<boost::shared_ptr<CaVessel<3> > > copied_vessels = vessel_network.CopyVessels();
        TS_ASSERT_EQUALS(vessel_network.GetNumberOfVessels(), 6u);

        // Move the new vessels
        c_vector<double, 3> translation_vector2 = zero_vector<double>(3);
        translation_vector2[2] = 3.0;
        vessel_network.Translate(translation_vector2);
        TS_ASSERT_DELTA(vessel_network.GetVessel(3)->GetSegment(0)->GetNode(1)->GetLocation()[0], 1.0, 1.e-6);
        TS_ASSERT_DELTA(vessel_network.GetVessel(3)->GetSegment(0)->GetNode(1)->GetLocation()[1], 2.0, 1.e-6);
        TS_ASSERT_DELTA(vessel_network.GetVessel(3)->GetSegment(0)->GetNode(1)->GetLocation()[2], 3.0, 1.e-6);
    }

    void TestConnnectedMethods() throw(Exception)
    {
        // Make some nodes
        std::vector<boost::shared_ptr<VascularNode<3> > > nodes;
        nodes.push_back(VascularNode<3>::Create(1.0, 2.0, 6.0));
        nodes.push_back(VascularNode<3>::Create(3.0, 4.0, 7.0));
        nodes.push_back(VascularNode<3>::Create(3.0, 4.0, 7.0));
        nodes.push_back(VascularNode<3>::Create(3.0, 4.0, 8.0));
        nodes.push_back(VascularNode<3>::Create(3.0, 4.0, 9.0));

        // Make some vessels
        boost::shared_ptr<CaVessel<3> > pVessel1(CaVessel<3>::Create(nodes[0], nodes[1]));
        boost::shared_ptr<CaVessel<3> > pVessel2(CaVessel<3>::Create(nodes[2], nodes[3]));
        boost::shared_ptr<CaVessel<3> > pVessel3(CaVessel<3>::Create(nodes[3], nodes[4]));

        std::vector<boost::shared_ptr<CaVessel<3> > > vessels;
        vessels.push_back(pVessel2);
        vessels.push_back(pVessel3);

        // Make a network
        CaVascularNetwork<3> vessel_network;
        vessel_network.AddVessel(pVessel1);
        vessel_network.AddVessels(vessels);

        // Test connectivity
        TS_ASSERT_EQUALS(vessel_network.GetNodes().size(), 5u);
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[0]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[1]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[2]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[3]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[4]));

        // Merge coincident nodes
        vessel_network.MergeCoincidentNodes();
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[0]));
        TS_ASSERT(!vessel_network.NodeIsInNetwork(nodes[1]) != !vessel_network.NodeIsInNetwork(nodes[2]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[3]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[4]));
        TS_ASSERT_EQUALS(vessel_network.GetNodes().size(), 4u);

        boost::shared_ptr<VascularNode<3> > p_node1 = VascularNode<3>::Create(1.0 , 1.0 , 1.0);
        boost::shared_ptr<VascularNode<3> > p_node2 = VascularNode<3>::Create(5.0 , 5.0 , 1.0);
        vessel_network.AddVessel(CaVessel<3>::Create(p_node1, p_node2));

        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[0]));
        // exclusive or (!A != !B)
        TS_ASSERT(!vessel_network.NodeIsInNetwork(nodes[1]) != !vessel_network.NodeIsInNetwork(nodes[2]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[3]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[4]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(p_node1));
        TS_ASSERT(vessel_network.NodeIsInNetwork(p_node2));

        TS_ASSERT(vessel_network.IsConnected(nodes[0], nodes[4]));
        TS_ASSERT(!vessel_network.IsConnected(nodes[0], p_node2));
        TS_ASSERT(vessel_network.IsConnected(p_node1, p_node2));

        std::vector<boost::shared_ptr<VascularNode<3> > > source_nodes;
        source_nodes.push_back(nodes[0]);
        source_nodes.push_back(p_node1);

        std::vector<boost::shared_ptr<VascularNode<3> > > query_nodes;
        query_nodes.push_back(nodes[0]);
        if (vessel_network.NodeIsInNetwork(nodes[1]))
        {
            query_nodes.push_back(nodes[1]);
        }
        if (vessel_network.NodeIsInNetwork(nodes[2]))
        {
            query_nodes.push_back(nodes[2]);
        }
        query_nodes.push_back(nodes[3]);
        query_nodes.push_back(nodes[4]);

        std::vector<bool> connected = vessel_network.IsConnected(source_nodes, query_nodes);
        TS_ASSERT(connected[0]);
        TS_ASSERT(connected[1]);
        TS_ASSERT(connected[2]);
        TS_ASSERT(connected[3]);

        OutputFileHandler output_file_handler("TestVesselNetwork",false);
        std::string output_filename4 = output_file_handler.GetOutputDirectoryFullPath().append("ConnectedTestVesselNetwork.gv");
        vessel_network.WriteConnectivity(output_filename4);
    }

    void TestRemoveVessel() throw(Exception)
    {
        // Make some nodes
        std::vector<boost::shared_ptr<VascularNode<3> > > nodes;
        nodes.push_back(VascularNode<3>::Create(0.0, 0.0, 0.0));
        nodes.push_back(VascularNode<3>::Create(20.0, 0.0, 0.0));
        nodes.push_back(VascularNode<3>::Create(30.0, 0.0, 0.0));
        nodes.push_back(VascularNode<3>::Create(50.0, 0.0, 0.0));

        // Make some vessels
        std::vector<boost::shared_ptr<CaVessel<3> > > vessels;
        for(unsigned idx=0; idx < 3; idx++)
        {
            vessels.push_back(CaVessel<3>::Create(CaVesselSegment<3>::Create(nodes[idx], nodes[idx+1])));
        }

        CaVascularNetwork<3> vessel_network;
        vessel_network.AddVessels(vessels);

        vessel_network.RemoveShortVessels(15.0, false);
        TS_ASSERT_EQUALS(vessel_network.GetNumberOfVessels(), 2u);
    }

    void TestMergeVessel() throw(Exception)
    {
        // Make some nodes
        std::vector<boost::shared_ptr<VascularNode<3> > > nodes;
        nodes.push_back(VascularNode<3>::Create(0.0, 0.0, 0.0));
        nodes.push_back(VascularNode<3>::Create(20.0, 0.0, 0.0));
        nodes.push_back(VascularNode<3>::Create(30.0, 0.0, 0.0));
        nodes.push_back(VascularNode<3>::Create(50.0, 0.0, 0.0));

        // Make some vessels
        std::vector<boost::shared_ptr<CaVessel<3> > > vessels;
        for(unsigned idx=0; idx < 3; idx++)
        {
            vessels.push_back(CaVessel<3>::Create(CaVesselSegment<3>::Create(nodes[idx], nodes[idx+1])));
        }

        CaVascularNetwork<3> vessel_network;
        vessel_network.AddVessels(vessels);

        vessel_network.MergeShortVessels(15.0);
        TS_ASSERT_EQUALS(vessel_network.GetNumberOfVessels(), 2u);
        TS_ASSERT_DELTA(vessels[0]->GetEndNode()->GetLocationVector()[0], 20.0, 1.e-6);
        TS_ASSERT_DELTA(vessels[2]->GetStartNode()->GetLocationVector()[0], 20.0, 1.e-6);
    }

    void TestDivideSingleVessel() throw(Exception)
    {
         // Make some nodes
         std::vector<boost::shared_ptr<VascularNode<3> > > nodes;
         for(unsigned idx=0; idx < 2; idx++)
         {
             nodes.push_back(VascularNode<3>::Create(2.0 * double(idx)));
         }

         // Make a network
         CaVascularNetwork<3> vessel_network;
         vessel_network.AddVessel(CaVessel<3>::Create(nodes[0], nodes[1]));

         TS_ASSERT_EQUALS(vessel_network.GetNumberOfVessels(), 1u);
         TS_ASSERT_EQUALS(vessel_network.GetNumberOfNodes(), 2u);

         // Do the divide
         vessel_network.DivideVessel(vessel_network.GetVessel(0), ChastePoint<3>(0.66, 0.0, 0.0));
         TS_ASSERT_EQUALS(vessel_network.GetNumberOfVessels(), 2u);
         TS_ASSERT_EQUALS(vessel_network.GetNumberOfNodes(), 3u);
         TS_ASSERT_DELTA(vessel_network.GetVessel(0)->GetSegment(0)->GetNode(0)->GetLocation()[0], 0.0, 1.e-6);
         TS_ASSERT_DELTA(vessel_network.GetVessel(0)->GetSegment(0)->GetNode(1)->GetLocation()[0], 0.66, 1.e-6);
         TS_ASSERT_DELTA(vessel_network.GetVessel(1)->GetSegment(0)->GetNode(0)->GetLocation()[0], 0.66, 1.e-6);
         TS_ASSERT_DELTA(vessel_network.GetVessel(1)->GetSegment(0)->GetNode(1)->GetLocation()[0], 2.0, 1.e-6);
    }

    void TestDividingSingleVesselWithMultipleSegments() throw(Exception)
    {

        // Make some nodes
        std::vector<boost::shared_ptr<VascularNode<3> > > nodes;
        for(unsigned idx=1; idx < 6; idx++)
        {
            nodes.push_back(VascularNode<3>::Create(double(idx), 0.0, 0.0));
        }

        // Generate the network
        CaVascularNetwork<3> p_vascular_network;
        p_vascular_network.AddVessel(CaVessel<3>::Create(nodes));

        TS_ASSERT_EQUALS(p_vascular_network.GetNumberOfVessels(), 1u);
        TS_ASSERT_EQUALS(p_vascular_network.GetNumberOfNodes(), 5u);

        // Do the divide
        p_vascular_network.DivideVessel(p_vascular_network.GetVessel(0), ChastePoint<3>(3.0, 0.0, 0.0));
        TS_ASSERT_EQUALS(p_vascular_network.GetNumberOfVessels(), 2u);
        TS_ASSERT_EQUALS(p_vascular_network.GetNumberOfNodes(), 5u);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(0)->GetSegment(0)->GetNode(0)->GetLocation()[0], 1.0, 1.e-6);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(0)->GetSegment(0)->GetNode(1)->GetLocation()[0], 2.0, 1.e-6);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(0)->GetSegment(1)->GetNode(0)->GetLocation()[0], 2.0, 1.e-6);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(0)->GetSegment(1)->GetNode(1)->GetLocation()[0], 3.0, 1.e-6);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(1)->GetSegment(0)->GetNode(0)->GetLocation()[0], 3.0, 1.e-6);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(1)->GetSegment(0)->GetNode(1)->GetLocation()[0], 4.0, 1.e-6);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(1)->GetSegment(1)->GetNode(0)->GetLocation()[0], 4.0, 1.e-6);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(1)->GetSegment(1)->GetNode(1)->GetLocation()[0], 5.0, 1.e-6);
        TS_ASSERT_THROWS_ANYTHING(p_vascular_network.DivideVessel(p_vascular_network.GetVessel(0), ChastePoint<3>(4.5, 0.0, 0.0)));

        // Do the divide
        p_vascular_network.DivideVessel(p_vascular_network.GetVessel(1), ChastePoint<3>(4.5, 0.0, 0.0));
        TS_ASSERT_EQUALS(p_vascular_network.GetNumberOfVessels(), 3u);
        TS_ASSERT_EQUALS(p_vascular_network.GetNumberOfNodes(), 6u);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(0)->GetSegment(0)->GetNode(0)->GetLocation()[0], 1.0, 1.e-6);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(0)->GetSegment(0)->GetNode(1)->GetLocation()[0], 2.0, 1.e-6);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(0)->GetSegment(1)->GetNode(0)->GetLocation()[0], 2.0, 1.e-6);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(0)->GetSegment(1)->GetNode(1)->GetLocation()[0], 3.0, 1.e-6);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(1)->GetSegment(0)->GetNode(0)->GetLocation()[0], 3.0, 1.e-6);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(1)->GetSegment(0)->GetNode(1)->GetLocation()[0], 4.0, 1.e-6);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(1)->GetSegment(1)->GetNode(0)->GetLocation()[0], 4.0, 1.e-6);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(1)->GetSegment(1)->GetNode(1)->GetLocation()[0], 4.5, 1.e-6);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(2)->GetSegment(0)->GetNode(0)->GetLocation()[0], 4.5, 1.e-6);
        TS_ASSERT_DELTA(p_vascular_network.GetVessel(2)->GetSegment(0)->GetNode(1)->GetLocation()[0], 5.0, 1.e-6);
    }

    void TestFormSproutWhichFormsAnastomosisWithExistingNodeImmediately()
    {

        boost::shared_ptr<CaVascularNetwork<2> > p_vessel_network(new CaVascularNetwork<2>());
        std::vector<boost::shared_ptr<VascularNode<2> > > nodes1;
        nodes1.push_back(VascularNode<2>::Create(4,0));
        nodes1.push_back(VascularNode<2>::Create(4,1));
        nodes1.push_back(VascularNode<2>::Create(4,2));
        boost::shared_ptr<CaVessel<2> > p_vessel1 = CaVessel<2>::Create(nodes1);

        std::vector<boost::shared_ptr<VascularNode<2> > > nodes2;
        nodes2.push_back(VascularNode<2>::Create(4,2));
        nodes2.push_back(VascularNode<2>::Create(5,3));
        nodes2.push_back(VascularNode<2>::Create(6,4));
        nodes2.push_back(VascularNode<2>::Create(7,5));
        nodes2.push_back(VascularNode<2>::Create(8,6));
        boost::shared_ptr<CaVessel<2> > p_vessel2 = CaVessel<2>::Create(nodes2);

        std::vector<boost::shared_ptr<VascularNode<2> > > nodes3;
        nodes3.push_back(VascularNode<2>::Create(8,6));
        nodes3.push_back(VascularNode<2>::Create(7,6));
        boost::shared_ptr<CaVessel<2> > p_vessel3 = CaVessel<2>::Create(nodes3);

        std::vector<boost::shared_ptr<VascularNode<2> > > nodes4;
        nodes4.push_back(VascularNode<2>::Create(7,6));
        nodes4.push_back(VascularNode<2>::Create(7,7));
        nodes4.push_back(VascularNode<2>::Create(7,8));
        nodes4.push_back(VascularNode<2>::Create(7,9));
        boost::shared_ptr<CaVessel<2> > p_vessel4 = CaVessel<2>::Create(nodes4);

        std::vector<boost::shared_ptr<VascularNode<2> > > nodes5;
        nodes5.push_back(VascularNode<2>::Create(7,6));
        nodes5.push_back(VascularNode<2>::Create(6,6));
        nodes5.push_back(VascularNode<2>::Create(5,6));
        nodes5.push_back(VascularNode<2>::Create(4,6));
        nodes5.push_back(VascularNode<2>::Create(3,6));
        nodes5.push_back(VascularNode<2>::Create(2,6));
        nodes5.push_back(VascularNode<2>::Create(1,6));
        nodes5.push_back(VascularNode<2>::Create(0,6));
        boost::shared_ptr<CaVessel<2> > p_vessel5 = CaVessel<2>::Create(nodes5);

        std::vector<boost::shared_ptr<VascularNode<2> > > nodes6;
        nodes6.push_back(VascularNode<2>::Create(0,6));
        nodes6.push_back(VascularNode<2>::Create(1,5));
        nodes6.push_back(VascularNode<2>::Create(2,4));
        nodes6.push_back(VascularNode<2>::Create(3,3));
        nodes6.push_back(VascularNode<2>::Create(4,2));
        boost::shared_ptr<CaVessel<2> > p_vessel6 = CaVessel<2>::Create(nodes6);

        p_vessel_network->AddVessel(p_vessel1);
        p_vessel_network->AddVessel(p_vessel2);
        p_vessel_network->AddVessel(p_vessel3);
        p_vessel_network->AddVessel(p_vessel4);
        p_vessel_network->AddVessel(p_vessel5);
        p_vessel_network->AddVessel(p_vessel6);

        p_vessel_network->MergeCoincidentNodes();

        // form sprout
        ChastePoint<2> sproutBaseLocation(7,7);
        ChastePoint<2> sproutTipLocation(8,6);
        boost::shared_ptr<CaVessel<2> > newSprout = p_vessel_network->FormSprout(sproutBaseLocation, sproutTipLocation);

        // check locations of nodes
        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,0)) == 1);
        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,2)) == 1);
        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(8,6)) == 2);
        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(7,6)) == 1);
        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(7,7)) == 1);
        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(7,9)) == 1);
        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(0,6)) == 1);

        // check number of vessels attached to each node
        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,0))->GetNumberOfSegments() == 1)
        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,2))->GetNumberOfSegments() == 3);
        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(7,6))->GetNumberOfSegments() == 3);
        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(7,7))->GetNumberOfSegments() == 3);
        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(7,9))->GetNumberOfSegments() == 1)
        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(0,6))->GetNumberOfSegments() == 2);

        // check active tips
        TS_ASSERT(newSprout->GetSegment(0)->GetNode(1)->IsMigrating());

        // test number of vessels and nodes in network
        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfNodes(),22u);
        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfVessels(),8u);

        // test that vessels are in network
        TS_ASSERT_THROWS_NOTHING(p_vessel_network->GetVesselIndex(p_vessel1));
        TS_ASSERT_THROWS_NOTHING(p_vessel_network->GetVesselIndex(p_vessel2));
        TS_ASSERT_THROWS_NOTHING(p_vessel_network->GetVesselIndex(p_vessel3));
        TS_ASSERT_THROWS_ANYTHING(p_vessel_network->GetVesselIndex(p_vessel4));
        TS_ASSERT_THROWS_NOTHING(p_vessel_network->GetVesselIndex(p_vessel5));
        TS_ASSERT_THROWS_NOTHING(p_vessel_network->GetVesselIndex(p_vessel6));
        TS_ASSERT_THROWS_NOTHING(p_vessel_network->GetVesselIndex(newSprout));

        // test that nodes are in network
        TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel1->GetStartNode()));
        TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel1->GetEndNode()));
        TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel2->GetStartNode()));
        TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel2->GetEndNode()));
        TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel3->GetStartNode()));
        TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel3->GetEndNode()));
        TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel4->GetStartNode()));
        TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel4->GetEndNode()));
        TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel5->GetStartNode()));
        TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel5->GetEndNode()));
        TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel6->GetStartNode()));
        TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel6->GetEndNode()));
        TS_ASSERT(p_vessel_network->NodeIsInNetwork(newSprout->GetStartNode()));
        TS_ASSERT(p_vessel_network->NodeIsInNetwork(newSprout->GetEndNode()));
    }

};

#endif /*TESTVESSELNETWORK_HPP_*/
