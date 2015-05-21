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

class TestVesselNetwork : public AbstractCellBasedTestSuite
{
public:

    typedef boost::shared_ptr<VascularNode<2> > NodePtr2;
    typedef boost::shared_ptr<VascularNode<3> > NodePtr3;
    typedef boost::shared_ptr<CaVesselSegment<2> > SegmentPtr2;
    typedef boost::shared_ptr<CaVesselSegment<3> > SegmentPtr3;
    typedef boost::shared_ptr<CaVessel<2> > VesselPtr2;
    typedef boost::shared_ptr<CaVessel<3> > VesselPtr3;

    void TestAddingVessels() throw(Exception)
    {
        // Make some nodes
        std::vector<NodePtr3> nodes;
        for(unsigned idx=0; idx < 4; idx++)
        {
            nodes.push_back(VascularNode<3>::Create(double(idx), 0.0, 0.0));
        }

        // Make some vessels
        std::vector<VesselPtr3> vessels;
        for(unsigned idx=0; idx < 1; idx++)
        {
            vessels.push_back(CaVessel<3>::Create(CaVesselSegment<3>::Create(nodes[idx], nodes[idx+1])));
        }
        VesselPtr3 p_end_vessel = CaVessel<3>::Create(CaVesselSegment<3>::Create(nodes[2], nodes[3]));

        // Make a network
        CaVascularNetwork<3> vessel_network;
        vessel_network.AddVessels(vessels);
        vessel_network.AddVessel(p_end_vessel);

        TS_ASSERT_EQUALS(vessel_network.GetNodes().size(), 4u);
    }

    void TestSettingNetworkData() throw(Exception)
    {
        // Make some nodes
        std::vector<NodePtr3> nodes;
        for(unsigned idx=0; idx < 4; idx++)
        {
            nodes.push_back(VascularNode<3>::Create(double(idx), 0.0, 0.0));
        }

        // Make some vessels
        std::vector<VesselPtr3> vessels;
        for(unsigned idx=0; idx < 2; idx++)
        {
            vessels.push_back(CaVessel<3>::Create(CaVesselSegment<3>::Create(nodes[idx], nodes[idx+1])));
        }

        // Make a network
        CaVascularNetwork<3> vessel_network;
        vessel_network.AddVessels(vessels);

        // Make some network data
        VasculatureData data;
        double radius = 10.0;
        double haematocrit = 0.4;
        bool has_flow = true;
        unsigned some_index = 5u;
        data.SetData("Radius", radius);
        data.SetData("Haematocrit", haematocrit);
        data.SetData("Has Flow", has_flow);
        data.SetData("SomeIndex", some_index);
        vessel_network.SetVesselData(data);
    }

    void TestCopyingAndMovingNetwork() throw(Exception)
    {
        // Make some nodes
        std::vector<NodePtr3> nodes;
        for(unsigned idx=0; idx < 4; idx++)
        {
            nodes.push_back(VascularNode<3>::Create(double(idx), 0.0, 0.0));
        }

        // Make some vessels
        std::vector<VesselPtr3> vessels;
        for(unsigned idx=0; idx < 3; idx++)
        {
            vessels.push_back(CaVessel<3>::Create(CaVesselSegment<3>::Create(nodes[idx], nodes[idx+1])));
        }

        // Make a network
        CaVascularNetwork<3> vessel_network;
        vessel_network.AddVessels(vessels);

        // Move the network
        c_vector<double, 3> translation_vector;
        translation_vector[0] = 0.0;
        translation_vector[1] = 2.0;
        translation_vector[2] = 0.0;

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
        c_vector<double, 3> translation_vector2;
        translation_vector2[0] = 0.0;
        translation_vector2[1] = 0.0;
        translation_vector2[2] = 3.0;
        vessel_network.Translate(translation_vector2, copied_vessels);
        TS_ASSERT_DELTA(vessel_network.GetVessel(3)->GetSegment(0)->GetNode(1)->GetLocation()[0], 1.0, 1.e-6);
        TS_ASSERT_DELTA(vessel_network.GetVessel(3)->GetSegment(0)->GetNode(1)->GetLocation()[1], 2.0, 1.e-6);
        TS_ASSERT_DELTA(vessel_network.GetVessel(3)->GetSegment(0)->GetNode(1)->GetLocation()[2], 3.0, 1.e-6);

        // Write the network
        OutputFileHandler output_file_handler("TestVesselNetwork", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("CopiedTranslatedNetwork.vtp");
        vessel_network.Write(output_filename);
    }

    void TestConnnectedMethods() throw(Exception)
    {
        // Make some nodes
        std::vector<ChastePoint<3> > points;
        points.push_back(ChastePoint<3>(1.0, 2.0, 6.0));
        points.push_back(ChastePoint<3>(3.0, 4.0, 7.0));
        points.push_back(ChastePoint<3>(3.0, 4.0, 7.0));
        points.push_back(ChastePoint<3>(3.0, 4.0, 8.0));
        points.push_back(ChastePoint<3>(3.0, 4.0, 9.0));

        std::vector<NodePtr3> nodes;
        for(unsigned i=0; i < points.size(); i++)
        {
            nodes.push_back(NodePtr3 (VascularNode<3>::Create(points[i])));
        }
        nodes.push_back(NodePtr3 (VascularNode<3>::Create(points[1])));

        // Make some segments
        SegmentPtr3 pSegment1(CaVesselSegment<3>::Create(nodes[0], nodes[1]));
        SegmentPtr3 pSegment2(CaVesselSegment<3>::Create(nodes[2], nodes[3]));
        SegmentPtr3 pSegment3(CaVesselSegment<3>::Create(nodes[3], nodes[4]));

        // Make some vessels
        VesselPtr3 pVessel1(CaVessel<3>::Create(pSegment1));
        VesselPtr3 pVessel2(CaVessel<3>::Create(pSegment2));
        VesselPtr3 pVessel3(CaVessel<3>::Create(pSegment3));

        std::vector<VesselPtr3> vessels;
        vessels.push_back(pVessel2);
        vessels.push_back(pVessel3);

        // Make a network
        CaVascularNetwork<3> vessel_network;
        vessel_network.AddVessel(pVessel1);
        vessel_network.AddVessels(vessels);

        TS_ASSERT_EQUALS(vessel_network.GetNodes().size(), 5u);

        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[0]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[1]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[2]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[3]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[4]));

        vessel_network.MergeCoincidentNodes();

        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[0]));
        // exclusive or (!A != !B)
        TS_ASSERT(!vessel_network.NodeIsInNetwork(nodes[1]) != !vessel_network.NodeIsInNetwork(nodes[2]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[3]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[4]));

        TS_ASSERT_EQUALS(vessel_network.GetNodes().size(), 4u);

        NodePtr3 node4 = VascularNode<3>::Create(1.0 , 1.0 , 1.0);
        NodePtr3 node5 = VascularNode<3>::Create(5.0 , 5.0 , 1.0);
        SegmentPtr3 pSegment4(CaVesselSegment<3>::Create(node4, node5));

        // Add another vessel
        VesselPtr3 pVessel4(CaVessel<3>::Create(pSegment4));
        vessel_network.AddVessel(pVessel4);

        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[0]));
        // exclusive or (!A != !B)
        TS_ASSERT(!vessel_network.NodeIsInNetwork(nodes[1]) != !vessel_network.NodeIsInNetwork(nodes[2]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[3]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[4]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(node4));
        TS_ASSERT(vessel_network.NodeIsInNetwork(node5));

        TS_ASSERT(vessel_network.IsConnected(nodes[0], nodes[4]));
        TS_ASSERT(!vessel_network.IsConnected(nodes[0], node5));
        TS_ASSERT(vessel_network.IsConnected(node4, node5));

        std::vector<NodePtr3> source_nodes;
        source_nodes.push_back(nodes[0]);
        source_nodes.push_back(node4);

        std::vector<NodePtr3> query_nodes;
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

    void TestDivideSingleVessel() throw(Exception)
    {
         // Make some nodes
         std::vector<NodePtr3> nodes;
         for(unsigned idx=0; idx < 2; idx++)
         {
             nodes.push_back(VascularNode<3>::Create(2.0 * double(idx)));
         }

         // Make a network
         CaVascularNetwork<3> vessel_network;
         vessel_network.AddVessel(CaVessel<3>::Create(CaVesselSegment<3>::Create(nodes[0], nodes[1])));

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
        std::vector<ChastePoint<3> > points;
        points.push_back(ChastePoint<3>(1.0, 0, 0));
        points.push_back(ChastePoint<3>(2.0, 0, 0));
        points.push_back(ChastePoint<3>(3.0, 0, 0));
        points.push_back(ChastePoint<3>(4.0, 0, 0));
        points.push_back(ChastePoint<3>(5.0, 0, 0));

        std::vector<NodePtr3> nodes;
        for(unsigned i=0; i < points.size(); i++)
        {
            nodes.push_back(NodePtr3 (VascularNode<3>::Create(points[i])));
        }

        SegmentPtr3 p_segment1(CaVesselSegment<3>::Create(nodes[0], nodes[1]));
        SegmentPtr3 p_segment2(CaVesselSegment<3>::Create(nodes[1], nodes[2]));
        SegmentPtr3 p_segment3(CaVesselSegment<3>::Create(nodes[2], nodes[3]));
        SegmentPtr3 p_segment4(CaVesselSegment<3>::Create(nodes[3], nodes[4]));

        std::vector<SegmentPtr3>  segments;
        segments.push_back(p_segment1);
        segments.push_back(p_segment2);
        segments.push_back(p_segment3);
        segments.push_back(p_segment4);

        VesselPtr3 p_vessel(CaVessel<3>::Create(segments));

        // Generate the network
        CaVascularNetwork<3> p_vascular_network;

        p_vascular_network.AddVessel(p_vessel);

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

    /**
     * Sprout formation tests
     */

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


//    void testFormSproutOnVesselWhichIsAttachedAtBothEndsToOtherVessels()
//    {
//
//        // set up spatial mesh
//
//        int DomainSize_X = 10;
//        int DomainSize_Y = 10;
//        int DomainSize_Z = 1;
//        double dx = 20*pow(10.0,-6); // metres
//        const int dimensionality = 2;
//
//        /*
//         * Set up SimulationTime object.  This must be set up in order for structural adaptation algorithm
//         * to run.
//         */
//        SimulationTime* p_simulation_time = SimulationTime::Instance();
//        p_simulation_time->SetStartTime(0.0);
//        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(30, 1);
//
//        boost::shared_ptr<OnLatticeMesh<dimensionality,int> > spatialmesh(new OnLatticeMesh<dimensionality,int>(dx,DomainSize_X,DomainSize_Y,DomainSize_Z));
//
//        // set up test vasular network
//
//        boost::shared_ptr<CaVascularNetwork<2> > vesselNetwork(new CaVascularNetwork<2>(spatialmesh));
//
//        boost::shared_ptr<Vessel<2,int> > vessel1(new Vessel<2,int>());
//        vessel1->SetNode1Location(ChastePoint<2>(4,0));
//        vessel1->SetNode2Location(ChastePoint<2>(4,2));
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,0), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,1), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,2), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel2(new Vessel<2,int>());
//        vessel2->SetNode1Location(ChastePoint<2>(4,2));
//        vessel2->SetNode2Location(ChastePoint<2>(8,6));
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,2), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,3), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(6,4), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,5), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(8,6), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel3(new Vessel<2,int>());
//        vessel3->SetNode1Location(ChastePoint<2>(8,6));
//        vessel3->SetNode2Location(ChastePoint<2>(7,6));
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(8,6), spatialmesh->GetSpatialMeshSize());
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,6), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel4(new Vessel<2,int>());
//        vessel4->SetNode1Location(ChastePoint<2>(7,6));
//        vessel4->SetNode2Location(ChastePoint<2>(7,9));
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,6), spatialmesh->GetSpatialMeshSize());
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,7), spatialmesh->GetSpatialMeshSize());
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,8), spatialmesh->GetSpatialMeshSize());
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,9), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel5(new Vessel<2,int>());
//        vessel5->SetNode1Location(ChastePoint<2>(7,6));
//        vessel5->SetNode2Location(ChastePoint<2>(0,6));
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,6), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(6,6), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,6), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,6), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(3,6), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(2,6), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,6), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(0,6), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel6(new Vessel<2,int>());
//        vessel6->SetNode1Location(ChastePoint<2>(0,6));
//        vessel6->SetNode2Location(ChastePoint<2>(4,2));
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(0,6), spatialmesh->GetSpatialMeshSize());
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,5), spatialmesh->GetSpatialMeshSize());
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(2,4), spatialmesh->GetSpatialMeshSize());
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(3,3), spatialmesh->GetSpatialMeshSize());
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,2), spatialmesh->GetSpatialMeshSize());
//
//
//        p_vessel_network->AddVessel(vessel1);
//        p_vessel_network->AddVessel(vessel2);
//        p_vessel_network->AddVessel(vessel3);
//        p_vessel_network->AddVessel(vessel4);
//        p_vessel_network->AddVessel(vessel5);
//        p_vessel_network->AddVessel(vessel6);
//
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel3->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel3->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel4->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel4->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel5->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel5->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel6->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel6->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        p_vessel_network->SetOutputNode(ChastePoint<2>(7,9));
//        p_vessel_network->SetInputNode(ChastePoint<2>(4,0));
//
//        // test initial state of vascular network
//
//        // check locations of nodes
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,2)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(8,6)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(7,6)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(7,9)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(0,6)) == 1);
//
//        // check number of vessels attached to each node
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,0))->GetNumberOfSegments() == 1)
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,2))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(8,6))->GetNumberOfSegments() == 2)
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(7,6))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(7,9))->GetNumberOfSegments() == 1)
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(0,6))->GetNumberOfSegments() == 2);
//
//        // check active tips
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel3->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel3->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel4->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel4->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel5->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel5->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel6->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel6->ActiveTipCellLocatedAtNode2() == false);
//
//        // test number of vessels and nodes in network
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfNodes(),6);
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfVessels(),6);
//
//        // test that vessels are in network
//
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel1));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel2));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel3));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel4));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel5));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel6));
//
//        // test that nodes are in network
//
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel4->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel4->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel5->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel5->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel6->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel6->GetEndNode()));
//
//        p_vessel_network->CheckForConsistency();
//
//        // perform structural adaptation algorithm on network
//
//        boost::shared_ptr<SimpleStructuralAdaptationAlgorithm<2,int> > structuralAdaptationAlgorithm(new SimpleStructuralAdaptationAlgorithm<2,int>());
//        // override default calculations
//        boost::shared_ptr<ConstantHaematocritCalculation<dimensionality,int> > haematocritcalculation(new ConstantHaematocritCalculation<dimensionality,int>());
//        structuralAdaptationAlgorithm->SetHaematocritCalculation(haematocritcalculation);
//
//        p_vessel_network->accept(structuralAdaptationAlgorithm);
//
//        string testoutputdirectory = "VesselNetworkTests/";
//        string testsubdirectory = "FormSproutOnVesselWhichIsAttachedAtBothEndsToOtherVessels/";
//        string directory;
//        directory.append(testoutputdirectory);
//        directory.append(testsubdirectory);
//
//        OutputFileHandler output_file_handler(directory, true);
//
//        string initialStateFilename;
//        initialStateFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        initialStateFilename.append("VesselNetworkDataInitial.vtk");
//
//        string initialDataSummaryFilename;
//        initialDataSummaryFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        initialDataSummaryFilename.append("InitialVesselNetworkDataSummary.txt");
//
//        string finalStateFilename;
//        finalStateFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        finalStateFilename.append("VesselNetworkDataFinal.vtk");
//
//        string finalDataSummaryFilename;
//        finalDataSummaryFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        finalDataSummaryFilename.append("FinalVesselNetworkDataSummary.txt");
//
//        // print out and save some data about the vessel network before forrmation of sprout
//
//        p_vessel_network->SaveVasculatureDataToFile(initialStateFilename);
//
//        PrintVascularNetworkDataToFile(vesselNetwork, initialDataSummaryFilename);
//
//        // form sprout
//
//        ChastePoint<2> sproutBaseLocation(4,6,0);
//        ChastePoint<2> sproutTipLocation(4,7,0);
//
//        boost::shared_ptr<Vessel<2,int> > newSprout(new Vessel<2,int>());
//
//        newSprout = p_vessel_network->FormSprout(p_vessel_network->GetVessel(ChastePoint<2>(4,6,0),0), sproutBaseLocation, sproutTipLocation);
//
//        // perform structural adaptation algorithm on new vessel network
//
//        p_vessel_network->accept(structuralAdaptationAlgorithm);
//
//        // print out and save some data about the vessel network after forrmation of sprout
//
//        p_vessel_network->SaveVasculatureDataToFile(finalStateFilename);
//
//        PrintVascularNetworkDataToFile(vesselNetwork, finalDataSummaryFilename);
//
//        // test final state of vascular network
//
//        // check locations of nodes
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,2)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(8,6)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(7,6)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,7)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,6)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(7,9)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(0,6)) == 1);
//
//        // check number of vessels attached to each node
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,0))->GetNumberOfSegments() == 1)
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,2))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(8,6))->GetNumberOfSegments() == 2)
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(7,6))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,7))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,6))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(7,9))->GetNumberOfSegments() == 1)
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(0,6))->GetNumberOfSegments() == 2);
//
//        // check active tips
//        for (int i = 0; i < p_vessel_network->GetNumberOfVessels(); i++)
//        {
//            if (newSprout == p_vessel_network->GetVessel(i))
//            {
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode1() == false);
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode2() == true);
//            }
//            else
//            {
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode1() == false);
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode2() == false);
//            }
//        }
//
//        // test number of vessels and nodes in network
//
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfNodes(),8);
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfVessels(),8);
//
//        // test that vessels are in network
//
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel1));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel2));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel3));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel4));
//        TS_ASSERT(!p_vessel_network->VesselIsInNetwork(vessel5));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel6));
//
//        // test that nodes are in network
//
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel4->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel4->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel5->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel5->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel6->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel6->GetEndNode()));
//
//        // check for consistency between node and vessel circular references
//
//        for (int i = 0; i < p_vessel_network->GetNumberOfVessels(); i++)
//        {
//            TS_ASSERT(p_vessel_network->GetVessel(i)->GetStartNode()->IsAttachedToVessel(p_vessel_network->GetVessel(i)));
//            TS_ASSERT(p_vessel_network->GetVessel(i)->GetEndNode()->IsAttachedToVessel(p_vessel_network->GetVessel(i)));
//            TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel_network->GetVessel(i)->GetStartNode()));
//            TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel_network->GetVessel(i)->GetEndNode()));
//        }
//
//        for (int i = 0; i < p_vessel_network->GetNumberOfNodes(); i++)
//        {
//            for (int j = 0; j < p_vessel_network->GetNode(i)->GetNumberOfSegments(); j++)
//            {
//                TS_ASSERT(p_vessel_network->GetNode(i)->GetAdjoiningVessel(j)->IsAttachedToNode(p_vessel_network->GetNode(i)));
//                TS_ASSERT(p_vessel_network->VesselIsInNetwork(p_vessel_network->GetNode(i)->GetAdjoiningVessel(j)));
//            }
//        }
//
//        p_vessel_network->CheckForConsistency();
//
//        SimulationTime::Destroy();
//    }
//
//    void testFormSproutOnVesselWhichIsAttachedAtOneEndToOtherVessels()
//    {
//
//        // set up spatial mesh
//
//        int DomainSize_X = 10;
//        int DomainSize_Y = 10;
//        int DomainSize_Z = 1;
//        double dx = 20*pow(10.0,-6); // metres
//        const int dimensionality = 2;
//
//        /*
//         * Set up SimulationTime object.  This must be set up in order for structural adaptation algorithm
//         * to run.
//         */
//        SimulationTime* p_simulation_time = SimulationTime::Instance();
//        p_simulation_time->SetStartTime(0.0);
//        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(30, 1);
//
//        boost::shared_ptr<OnLatticeMesh<dimensionality,int> > spatialmesh(new OnLatticeMesh<dimensionality,int>(dx,DomainSize_X,DomainSize_Y,DomainSize_Z));
//
//        // set up test vasular network
//
//        boost::shared_ptr<CaVascularNetwork<2> > vesselNetwork(new CaVascularNetwork<2>(spatialmesh));
//
//        boost::shared_ptr<Vessel<2,int> > vessel1(new Vessel<2,int>());
//        vessel1->SetNode1Location(ChastePoint<2>(4,0));
//        vessel1->SetNode2Location(ChastePoint<2>(4,2));
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,0), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,1), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,2), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel2(new Vessel<2,int>());
//        vessel2->SetNode1Location(ChastePoint<2>(4,2));
//        vessel2->SetNode2Location(ChastePoint<2>(8,6));
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,2), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,3), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(6,4), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,5), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(8,6), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel3(new Vessel<2,int>());
//        vessel3->SetNode1Location(ChastePoint<2>(8,6));
//        vessel3->SetNode2Location(ChastePoint<2>(7,6));
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(8,6), spatialmesh->GetSpatialMeshSize());
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,6), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel4(new Vessel<2,int>());
//        vessel4->SetNode1Location(ChastePoint<2>(7,6));
//        vessel4->SetNode2Location(ChastePoint<2>(7,9));
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,6), spatialmesh->GetSpatialMeshSize());
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,7), spatialmesh->GetSpatialMeshSize());
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,8), spatialmesh->GetSpatialMeshSize());
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,9), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel5(new Vessel<2,int>());
//        vessel5->SetNode1Location(ChastePoint<2>(7,6));
//        vessel5->SetNode2Location(ChastePoint<2>(0,6));
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,6), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(6,6), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,6), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,6), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(3,6), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(2,6), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,6), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(0,6), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel6(new Vessel<2,int>());
//        vessel6->SetNode1Location(ChastePoint<2>(0,6));
//        vessel6->SetNode2Location(ChastePoint<2>(4,2));
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(0,6), spatialmesh->GetSpatialMeshSize());
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,5), spatialmesh->GetSpatialMeshSize());
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(2,4), spatialmesh->GetSpatialMeshSize());
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(3,3), spatialmesh->GetSpatialMeshSize());
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,2), spatialmesh->GetSpatialMeshSize());
//
//
//        p_vessel_network->AddVessel(vessel1);
//        p_vessel_network->AddVessel(vessel2);
//        p_vessel_network->AddVessel(vessel3);
//        p_vessel_network->AddVessel(vessel4);
//        p_vessel_network->AddVessel(vessel5);
//        p_vessel_network->AddVessel(vessel6);
//
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel3->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel3->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel4->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel4->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel5->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel5->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel6->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel6->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        p_vessel_network->SetOutputNode(ChastePoint<2>(7,9));
//        p_vessel_network->SetInputNode(ChastePoint<2>(4,0));
//
//        // test initial state of vascular network
//
//        // check locations of nodes
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,2)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(8,6)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(7,6)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(7,9)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(0,6)) == 1);
//
//        // check number of vessels attached to each node
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,0))->GetNumberOfSegments() == 1)
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,2))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(8,6))->GetNumberOfSegments() == 2)
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(7,6))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(7,9))->GetNumberOfSegments() == 1)
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(0,6))->GetNumberOfSegments() == 2);
//
//        // check active tips
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel3->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel3->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel4->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel4->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel5->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel5->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel6->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel6->ActiveTipCellLocatedAtNode2() == false);
//
//        // test number of vessels and nodes in network
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfNodes(),6);
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfVessels(),6);
//
//        // test that vessels are in network
//
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel1));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel2));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel3));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel4));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel5));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel6));
//
//        // test that nodes are in network
//
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel4->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel4->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel5->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel5->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel6->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel6->GetEndNode()));
//
//        p_vessel_network->CheckForConsistency();
//
//        // perform structural adaptation algorithm on network
//
//        boost::shared_ptr<SimpleStructuralAdaptationAlgorithm<2,int> > structuralAdaptationAlgorithm(new SimpleStructuralAdaptationAlgorithm<2,int>());
//        // override default calculations
//        boost::shared_ptr<ConstantHaematocritCalculation<dimensionality,int> > haematocritcalculation(new ConstantHaematocritCalculation<dimensionality,int>());
//        structuralAdaptationAlgorithm->SetHaematocritCalculation(haematocritcalculation);
//
//        p_vessel_network->accept(structuralAdaptationAlgorithm);
//
//        string testoutputdirectory = "VesselNetworkTests/";
//        string testsubdirectory = "FormSproutOnVesselWhichIsAttachedAtOneEndToOtherVessels/";
//        string directory;
//        directory.append(testoutputdirectory);
//        directory.append(testsubdirectory);
//
//        OutputFileHandler output_file_handler(directory, true);
//
//        string initialStateFilename;
//        initialStateFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        initialStateFilename.append("VesselNetworkDataInitial.vtk");
//
//        string initialDataSummaryFilename;
//        initialDataSummaryFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        initialDataSummaryFilename.append("InitialVesselNetworkDataSummary.txt");
//
//        string finalStateFilename;
//        finalStateFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        finalStateFilename.append("VesselNetworkDataFinal.vtk");
//
//        string finalDataSummaryFilename;
//        finalDataSummaryFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        finalDataSummaryFilename.append("FinalVesselNetworkDataSummary.txt");
//
//        // print out and save some data about the vessel network before forrmation of sprout
//
//        p_vessel_network->SaveVasculatureDataToFile(initialStateFilename);
//
//        PrintVascularNetworkDataToFile(vesselNetwork, initialDataSummaryFilename);
//
//        // form sprout
//
//        ChastePoint<2> sproutBaseLocation(7,7,0);
//        ChastePoint<2> sproutTipLocation(8,7,0);
//
//        boost::shared_ptr<Vessel<2,int> > newSprout(new Vessel<2,int>());
//
//        newSprout = p_vessel_network->FormSprout(p_vessel_network->GetVessel(ChastePoint<2>(7,7,0),0), sproutBaseLocation, sproutTipLocation);
//
//        // perform structural adaptation algorithm on new vessel network
//
//        p_vessel_network->accept(structuralAdaptationAlgorithm);
//
//        // print out and save some data about the vessel network after forrmation of sprout
//
//        p_vessel_network->SaveVasculatureDataToFile(finalStateFilename);
//
//        PrintVascularNetworkDataToFile(vesselNetwork, finalDataSummaryFilename);
//
//        // test final state of vascular network
//
//        // check locations of nodes
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,2)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(8,6)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(7,6)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(7,7)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(8,7)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(7,9)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(0,6)) == 1);
//
//        // check number of vessels attached to each node
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,0))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,2))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(8,6))->GetNumberOfSegments() == 2);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(7,6))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(7,7))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(8,7))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(7,9))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(0,6))->GetNumberOfSegments() == 2);
//
//        // check active tips
//        for (int i = 0; i < p_vessel_network->GetNumberOfVessels(); i++)
//        {
//            if (newSprout == p_vessel_network->GetVessel(i))
//            {
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode1() == false);
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode2() == true);
//            }
//            else
//            {
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode1() == false);
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode2() == false);
//            }
//        }
//
//        // test number of vessels and nodes in network
//
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfNodes(),8);
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfVessels(),8);
//
//        // test that vessels are in network
//
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel1));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel2));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel3));
//        TS_ASSERT(!p_vessel_network->VesselIsInNetwork(vessel4));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel5));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel6));
//
//        // test that nodes are in network
//
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel4->GetStartNode()));
//        TS_ASSERT(!p_vessel_network->NodeIsInNetwork(vessel4->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel5->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel5->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel6->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel6->GetEndNode()));
//
//        // check for consistency between node and vessel circular references
//
//        for (int i = 0; i < p_vessel_network->GetNumberOfVessels(); i++)
//        {
//            TS_ASSERT(p_vessel_network->GetVessel(i)->GetStartNode()->IsAttachedToVessel(p_vessel_network->GetVessel(i)));
//            TS_ASSERT(p_vessel_network->GetVessel(i)->GetEndNode()->IsAttachedToVessel(p_vessel_network->GetVessel(i)));
//            TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel_network->GetVessel(i)->GetStartNode()));
//            TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel_network->GetVessel(i)->GetEndNode()));
//        }
//
//        for (int i = 0; i < p_vessel_network->GetNumberOfNodes(); i++)
//        {
//            for (int j = 0; j < p_vessel_network->GetNode(i)->GetNumberOfSegments(); j++)
//            {
//                TS_ASSERT(p_vessel_network->GetNode(i)->GetAdjoiningVessel(j)->IsAttachedToNode(p_vessel_network->GetNode(i)));
//                TS_ASSERT(p_vessel_network->VesselIsInNetwork(p_vessel_network->GetNode(i)->GetAdjoiningVessel(j)));
//            }
//        }
//
//        p_vessel_network->CheckForConsistency();
//
//        SimulationTime::Destroy();
//
//    }
//
//    void testFormSproutOnVesselWhichFormsASelfLoop()
//    {
//
//        // set up spatial mesh
//
//        int DomainSize_X = 10;
//        int DomainSize_Y = 10;
//        int DomainSize_Z = 1;
//        double dx = 20*pow(10.0,-6); // metres
//        const int dimensionality = 2;
//
//        /*
//         * Set up SimulationTime object.  This must be set up in order for structural adaptation algorithm
//         * to run.
//         */
//        SimulationTime* p_simulation_time = SimulationTime::Instance();
//        p_simulation_time->SetStartTime(0.0);
//        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(30, 1);
//
//        boost::shared_ptr<OnLatticeMesh<dimensionality,int> > spatialmesh(new OnLatticeMesh<dimensionality,int>(dx,DomainSize_X,DomainSize_Y,DomainSize_Z));
//
//        // set up test vasular network
//
//        boost::shared_ptr<CaVascularNetwork<2> > vesselNetwork(new CaVascularNetwork<2>(spatialmesh));
//
//        boost::shared_ptr<Vessel<2,int> > vessel1(new Vessel<2,int>());
//        vessel1->SetNode1Location(ChastePoint<2>(1,0));
//        vessel1->SetNode2Location(ChastePoint<2>(1,4));
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,0), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,1), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,2), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,3), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,4), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel2(new Vessel<2,int>());
//        vessel2->SetNode1Location(ChastePoint<2>(1,4));
//        vessel2->SetNode2Location(ChastePoint<2>(1,7));
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,4), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,5), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,6), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,7), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel3(new Vessel<2,int>());
//        vessel3->SetNode1Location(ChastePoint<2>(1,7));
//        vessel3->SetNode2Location(ChastePoint<2>(1,9));
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,7), spatialmesh->GetSpatialMeshSize());
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,8), spatialmesh->GetSpatialMeshSize());
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,9), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel4(new Vessel<2,int>());
//        vessel4->SetNode1Location(ChastePoint<2>(1,4));
//        vessel4->SetNode2Location(ChastePoint<2>(3,4));
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,4), spatialmesh->GetSpatialMeshSize());
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(2,4), spatialmesh->GetSpatialMeshSize());
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(3,4), spatialmesh->GetSpatialMeshSize());
//
//        boost::shared_ptr<Vessel<2,int> > vessel5(new Vessel<2,int>());
//        vessel5->SetNode1Location(ChastePoint<2>(3,4));
//        vessel5->SetNode2Location(ChastePoint<2>(1,7));
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(3,4), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(3,5), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(3,6), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(3,7), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(2,7), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,7), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel6(new Vessel<2,int>());
//        vessel6->SetNode1Location(ChastePoint<2>(3,4));
//        vessel6->SetNode2Location(ChastePoint<2>(6,4));
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(3,4), spatialmesh->GetSpatialMeshSize());
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,4), spatialmesh->GetSpatialMeshSize());
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,4), spatialmesh->GetSpatialMeshSize());
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(6,4), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel7(new Vessel<2,int>());
//        vessel7->SetNode1Location(ChastePoint<2>(6,4));
//        vessel7->SetNode2Location(ChastePoint<2>(6,4));
//        vessel7->SetNextVesselSegmentCoordinate(ChastePoint<2>(6,4), spatialmesh->GetSpatialMeshSize());
//        vessel7->SetNextVesselSegmentCoordinate(ChastePoint<2>(6,5), spatialmesh->GetSpatialMeshSize());
//        vessel7->SetNextVesselSegmentCoordinate(ChastePoint<2>(6,6), spatialmesh->GetSpatialMeshSize());
//        vessel7->SetNextVesselSegmentCoordinate(ChastePoint<2>(6,7), spatialmesh->GetSpatialMeshSize());
//        vessel7->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,7), spatialmesh->GetSpatialMeshSize());
//        vessel7->SetNextVesselSegmentCoordinate(ChastePoint<2>(8,7), spatialmesh->GetSpatialMeshSize());
//        vessel7->SetNextVesselSegmentCoordinate(ChastePoint<2>(8,6), spatialmesh->GetSpatialMeshSize());
//        vessel7->SetNextVesselSegmentCoordinate(ChastePoint<2>(8,5), spatialmesh->GetSpatialMeshSize());
//        vessel7->SetNextVesselSegmentCoordinate(ChastePoint<2>(8,4), spatialmesh->GetSpatialMeshSize());
//        vessel7->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,4), spatialmesh->GetSpatialMeshSize());
//        vessel7->SetNextVesselSegmentCoordinate(ChastePoint<2>(6,4), spatialmesh->GetSpatialMeshSize());
//
//
//        p_vessel_network->AddVessel(vessel1);
//        p_vessel_network->AddVessel(vessel2);
//        p_vessel_network->AddVessel(vessel3);
//        p_vessel_network->AddVessel(vessel4);
//        p_vessel_network->AddVessel(vessel5);
//        p_vessel_network->AddVessel(vessel6);
//        p_vessel_network->AddVessel(vessel7);
//
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel3->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel3->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel4->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel4->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel5->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel5->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel6->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel6->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel7->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel7->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        p_vessel_network->SetOutputNode(ChastePoint<2>(1,0));
//        p_vessel_network->SetInputNode(ChastePoint<2>(1,9));
//
//
//        // test initial state of vascular network
//
//        // check locations of nodes
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(1,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(1,4)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(1,7)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(1,9)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(3,4)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(6,4)) == 1);
//
//        // check number of vessels attached to each node
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(1,0))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(1,4))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(1,7))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(1,9))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(3,4))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(6,4))->GetNumberOfSegments() == 3); // two of these are the same
//
//        // check active tips
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel3->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel3->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel4->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel4->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel5->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel5->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel6->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel6->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel7->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel7->ActiveTipCellLocatedAtNode2() == false);
//
//        // test number of vessels and nodes in network
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfNodes(),6);
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfVessels(),7);
//
//        // test that vessels are in network
//
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel1));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel2));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel3));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel4));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel5));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel6));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel7));
//
//        // test that nodes are in network
//
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel4->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel4->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel5->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel5->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel6->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel6->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel7->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel7->GetEndNode()));
//
//        p_vessel_network->CheckForConsistency();
//
//        // perform structural adaptation algorithm on network
//
//        boost::shared_ptr<SimpleStructuralAdaptationAlgorithm<2,int> > structuralAdaptationAlgorithm(new SimpleStructuralAdaptationAlgorithm<2,int>());
//        // override default calculations
//        boost::shared_ptr<ConstantHaematocritCalculation<dimensionality,int> > haematocritcalculation(new ConstantHaematocritCalculation<dimensionality,int>());
//        structuralAdaptationAlgorithm->SetHaematocritCalculation(haematocritcalculation);
//
//        p_vessel_network->accept(structuralAdaptationAlgorithm);
//
//        string testoutputdirectory = "VesselNetworkTests/";
//        string testsubdirectory = "FormSproutOnVesselWhichFormsASelfLoop/";
//        string directory;
//        directory.append(testoutputdirectory);
//        directory.append(testsubdirectory);
//
//        OutputFileHandler output_file_handler(directory, true);
//
//        string initialStateFilename;
//        initialStateFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        initialStateFilename.append("VesselNetworkDataInitial.vtk");
//
//        string initialDataSummaryFilename;
//        initialDataSummaryFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        initialDataSummaryFilename.append("InitialVesselNetworkDataSummary.txt");
//
//        string finalStateFilename;
//        finalStateFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        finalStateFilename.append("VesselNetworkDataFinal.vtk");
//
//        string finalDataSummaryFilename;
//        finalDataSummaryFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        finalDataSummaryFilename.append("FinalVesselNetworkDataSummary.txt");
//
//        // print out and save some data about the vessel network before forrmation of sprout
//
//        p_vessel_network->SaveVasculatureDataToFile(initialStateFilename);
//
//        PrintVascularNetworkDataToFile(vesselNetwork, initialDataSummaryFilename);
//
//        // form sprout
//
//        ChastePoint<2> sproutBaseLocation(8,7,0);
//        ChastePoint<2> sproutTipLocation(9,8,0);
//
//        boost::shared_ptr<Vessel<2,int> > newSprout(new Vessel<2,int>());
//
//        newSprout = p_vessel_network->FormSprout(p_vessel_network->GetVessel(ChastePoint<2>(8,7,0),0), sproutBaseLocation, sproutTipLocation);
//
//        // perform structural adaptation algorithm on new vessel network
//
//        p_vessel_network->accept(structuralAdaptationAlgorithm);
//
//        // print out and save some data about the vessel network after forrmation of sprout
//
//        p_vessel_network->SaveVasculatureDataToFile(finalStateFilename);
//
//        PrintVascularNetworkDataToFile(vesselNetwork, finalDataSummaryFilename);
//
//        // test final state of vascular network
//
//        // check locations of nodes
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(1,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(1,4)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(1,7)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(1,9)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(3,4)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(6,4)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(8,7)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(9,8)) == 1);
//
//        // check number of vessels attached to each node
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(1,0))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(1,4))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(1,7))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(1,9))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(3,4))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(6,4))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(8,7))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(9,8))->GetNumberOfSegments() == 1);
//
//        // check active tips
//        for (int i = 0; i < p_vessel_network->GetNumberOfVessels(); i++)
//        {
//            if (newSprout == p_vessel_network->GetVessel(i))
//            {
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode1() == false);
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode2() == true);
//            }
//            else
//            {
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode1() == false);
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode2() == false);
//            }
//        }
//
//        // test number of vessels and nodes in network
//
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfNodes(),8);
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfVessels(),9);
//
//        // test that vessels are in network
//
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel1));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel2));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel3));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel4));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel5));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel6));
//        TS_ASSERT(!p_vessel_network->VesselIsInNetwork(vessel7));
//
//        // test that nodes are in network
//
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel4->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel4->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel5->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel5->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel6->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel6->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel7->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel7->GetEndNode()));
//
//        // check for consistency between node and vessel circular references
//
//        for (int i = 0; i < p_vessel_network->GetNumberOfVessels(); i++)
//        {
//            TS_ASSERT(p_vessel_network->GetVessel(i)->GetStartNode()->IsAttachedToVessel(p_vessel_network->GetVessel(i)));
//            TS_ASSERT(p_vessel_network->GetVessel(i)->GetEndNode()->IsAttachedToVessel(p_vessel_network->GetVessel(i)));
//            TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel_network->GetVessel(i)->GetStartNode()));
//            TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel_network->GetVessel(i)->GetEndNode()));
//        }
//
//        for (int i = 0; i < p_vessel_network->GetNumberOfNodes(); i++)
//        {
//            for (int j = 0; j < p_vessel_network->GetNode(i)->GetNumberOfSegments(); j++)
//            {
//                TS_ASSERT(p_vessel_network->GetNode(i)->GetAdjoiningVessel(j)->IsAttachedToNode(p_vessel_network->GetNode(i)));
//                TS_ASSERT(p_vessel_network->VesselIsInNetwork(p_vessel_network->GetNode(i)->GetAdjoiningVessel(j)));
//            }
//        }
//
//        p_vessel_network->CheckForConsistency();
//
//        SimulationTime::Destroy();
//
//    }
//
//    void testFormSproutOnNodeWhereTwoVesselsIntersect()
//    {
//
//        // set up spatial mesh
//
//        int DomainSize_X = 10;
//        int DomainSize_Y = 10;
//        int DomainSize_Z = 1;
//        double dx = 20*pow(10.0,-6); // metres
//        const int dimensionality = 2;
//
//        /*
//         * Set up SimulationTime object.  This must be set up in order for structural adaptation algorithm
//         * to run.
//         */
//        SimulationTime* p_simulation_time = SimulationTime::Instance();
//        p_simulation_time->SetStartTime(0.0);
//        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(30, 1);
//
//        boost::shared_ptr<OnLatticeMesh<dimensionality,int> > spatialmesh(new OnLatticeMesh<dimensionality,int>(dx,DomainSize_X,DomainSize_Y,DomainSize_Z));
//        spatialmesh->SetVolumeFractionOccupiedByVessel(0.3);
//        // set up test vasular network
//
//        boost::shared_ptr<CaVascularNetwork<2> > vesselNetwork(new CaVascularNetwork<2>(spatialmesh));
//
//        boost::shared_ptr<Vessel<2,int> > vessel1(new Vessel<2,int>());
//        vessel1->SetNode1Location(ChastePoint<2>(4,0));
//        vessel1->SetNode2Location(ChastePoint<2>(4,5));
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,0), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,1), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,2), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,3), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,4), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,5), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel2(new Vessel<2,int>());
//        vessel2->SetNode1Location(ChastePoint<2>(4,5));
//        vessel2->SetNode2Location(ChastePoint<2>(4,9));
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,5), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,6), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,7), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,8), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,9), spatialmesh->GetSpatialMeshSize());
//
//
//
//
//        p_vessel_network->AddVessel(vessel1);
//        p_vessel_network->AddVessel(vessel2);
//
//
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//
//        p_vessel_network->SetOutputNode(ChastePoint<2>(4,9));
//        p_vessel_network->SetInputNode(ChastePoint<2>(4,0));
//
//        // test initial state of vascular network
//
//        // check locations of nodes
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,5)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,9)) == 1);
//
//
//        // check number of vessels attached to each node
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,0))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,5))->GetNumberOfSegments() == 2);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,9))->GetNumberOfSegments() == 1);
//
//
//        // check active tips
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode2() == false);
//
//
//        // test number of vessels and nodes in network
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfNodes(),3);
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfVessels(),2);
//
//        // test that vessels are in network
//
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel1));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel2));
//
//
//        // test that nodes are in network
//
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetEndNode()));
//
//
//        p_vessel_network->CheckForConsistency();
//
//        // perform structural adaptation algorithm on network
//
//        boost::shared_ptr<SimpleStructuralAdaptationAlgorithm<2,int> > structuralAdaptationAlgorithm(new SimpleStructuralAdaptationAlgorithm<2,int>());
//        // override default calculations
//        boost::shared_ptr<ConstantHaematocritCalculation<dimensionality,int> > haematocritcalculation(new ConstantHaematocritCalculation<dimensionality,int>());
//        structuralAdaptationAlgorithm->SetHaematocritCalculation(haematocritcalculation);
//
//        p_vessel_network->accept(structuralAdaptationAlgorithm);
//
//        string testoutputdirectory = "VesselNetworkTests/";
//        string testsubdirectory = "FormSproutOnNodeWhereTwoVesselsIntersect/";
//        string directory;
//        directory.append(testoutputdirectory);
//        directory.append(testsubdirectory);
//
//        OutputFileHandler output_file_handler(directory, true);
//
//        string initialStateFilename;
//        initialStateFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        initialStateFilename.append("VesselNetworkDataInitial.vtk");
//
//        string initialDataSummaryFilename;
//        initialDataSummaryFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        initialDataSummaryFilename.append("InitialVesselNetworkDataSummary.txt");
//
//        string finalStateFilename;
//        finalStateFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        finalStateFilename.append("VesselNetworkDataFinal.vtk");
//
//        string finalDataSummaryFilename;
//        finalDataSummaryFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        finalDataSummaryFilename.append("FinalVesselNetworkDataSummary.txt");
//
//        // print out and save some data about the vessel network before forrmation of sprout
//
//        p_vessel_network->SaveVasculatureDataToFile(initialStateFilename);
//
//        PrintVascularNetworkDataToFile(vesselNetwork, initialDataSummaryFilename);
//
//        // form sprout
//
//        ChastePoint<2> sproutBaseLocation(4,5,0);
//        ChastePoint<2> sproutTipLocation(5,5,0);
//
//        boost::shared_ptr<Vessel<2,int> > newSprout(new Vessel<2,int>());
//
//        newSprout = p_vessel_network->FormSprout(p_vessel_network->GetVessel(ChastePoint<2>(4,5,0),0), sproutBaseLocation, sproutTipLocation);
//
//        // perform structural adaptation algorithm on new vessel network
//
//        p_vessel_network->accept(structuralAdaptationAlgorithm);
//
//        // print out and save some data about the vessel network after forrmation of sprout
//
//        p_vessel_network->SaveVasculatureDataToFile(finalStateFilename);
//
//        PrintVascularNetworkDataToFile(vesselNetwork, finalDataSummaryFilename);
//
//        // test final state of vascular network
//
//        // check locations of nodes
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,5)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(5,5)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,9)) == 1);
//
//
//        // check number of vessels attached to each node
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,0))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,5))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(5,5))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,9))->GetNumberOfSegments() == 1);
//
//        // check active tips
//        for (int i = 0; i < p_vessel_network->GetNumberOfVessels(); i++)
//        {
//            if (newSprout == p_vessel_network->GetVessel(i))
//            {
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode1() == false);
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode2() == true);
//            }
//            else
//            {
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode1() == false);
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode2() == false);
//            }
//        }
//
//        // test number of vessels and nodes in network
//
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfNodes(),4);
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfVessels(),3);
//
//        // test that vessels are in network
//
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel1));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel2));
//
//
//        // test that nodes are in network
//
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetEndNode()));
//
//
//        // check for consistency between node and vessel circular references
//
//        for (int i = 0; i < p_vessel_network->GetNumberOfVessels(); i++)
//        {
//            TS_ASSERT(p_vessel_network->GetVessel(i)->GetStartNode()->IsAttachedToVessel(p_vessel_network->GetVessel(i)));
//            TS_ASSERT(p_vessel_network->GetVessel(i)->GetEndNode()->IsAttachedToVessel(p_vessel_network->GetVessel(i)));
//            TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel_network->GetVessel(i)->GetStartNode()));
//            TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel_network->GetVessel(i)->GetEndNode()));
//        }
//
//        for (int i = 0; i < p_vessel_network->GetNumberOfNodes(); i++)
//        {
//            for (int j = 0; j < p_vessel_network->GetNode(i)->GetNumberOfSegments(); j++)
//            {
//                TS_ASSERT(p_vessel_network->GetNode(i)->GetAdjoiningVessel(j)->IsAttachedToNode(p_vessel_network->GetNode(i)));
//                TS_ASSERT(p_vessel_network->VesselIsInNetwork(p_vessel_network->GetNode(i)->GetAdjoiningVessel(j)));
//            }
//        }
//
//        p_vessel_network->CheckForConsistency();
//
//        SimulationTime::Destroy();
//
//    }
//
//    void testFormSproutOnTipCell()
//    {
//
//        // set up spatial mesh
//
//        int DomainSize_X = 10;
//        int DomainSize_Y = 10;
//        int DomainSize_Z = 1;
//        double dx = 20*pow(10.0,-6); // metres
//        const int dimensionality = 2;
//
//        /*
//         * Set up SimulationTime object.  This must be set up in order for structural adaptation algorithm
//         * to run.
//         */
//        SimulationTime* p_simulation_time = SimulationTime::Instance();
//        p_simulation_time->SetStartTime(0.0);
//        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(30, 1);
//
//        boost::shared_ptr<OnLatticeMesh<dimensionality,int> > spatialmesh(new OnLatticeMesh<dimensionality,int>(dx,DomainSize_X,DomainSize_Y,DomainSize_Z));
//        spatialmesh->SetVolumeFractionOccupiedByVessel(0.3);
//        // set up test vasular network
//
//        boost::shared_ptr<CaVascularNetwork<2> > vesselNetwork(new CaVascularNetwork<2>(spatialmesh));
//
//        boost::shared_ptr<Vessel<2,int> > vessel1(new Vessel<2,int>());
//        vessel1->SetNode1Location(ChastePoint<2>(4,0));
//        vessel1->SetNode2Location(ChastePoint<2>(4,5));
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,0), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,1), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,2), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,3), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,4), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,5), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel2(new Vessel<2,int>());
//        vessel2->SetNode1Location(ChastePoint<2>(4,5));
//        vessel2->SetNode2Location(ChastePoint<2>(4,9));
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,5), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,6), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,7), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,8), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,9), spatialmesh->GetSpatialMeshSize());
//
//
//        p_vessel_network->AddVessel(vessel1);
//        p_vessel_network->AddVessel(vessel2);
//
//
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//
//        p_vessel_network->SetOutputNode(ChastePoint<2>(4,9));
//        p_vessel_network->SetInputNode(ChastePoint<2>(4,0));
//
//        // test initial state of vascular network
//
//        // check locations of nodes
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,5)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,9)) == 1);
//
//
//        // check number of vessels attached to each node
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,0))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,5))->GetNumberOfSegments() == 2);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,9))->GetNumberOfSegments() == 1);
//
//
//        // check active tips
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode2() == false);
//
//
//        // test number of vessels and nodes in network
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfNodes(),3);
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfVessels(),2);
//
//        // test that vessels are in network
//
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel1));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel2));
//
//
//        // test that nodes are in network
//
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetEndNode()));
//
//
//        p_vessel_network->CheckForConsistency();
//
//        // perform structural adaptation algorithm on network
//
//        boost::shared_ptr<SimpleStructuralAdaptationAlgorithm<2,int> > structuralAdaptationAlgorithm(new SimpleStructuralAdaptationAlgorithm<2,int>());
//        // override default calculations
//        boost::shared_ptr<ConstantHaematocritCalculation<dimensionality,int> > haematocritcalculation(new ConstantHaematocritCalculation<dimensionality,int>());
//        structuralAdaptationAlgorithm->SetHaematocritCalculation(haematocritcalculation);
//
//        p_vessel_network->accept(structuralAdaptationAlgorithm);
//
//        // form sprout
//
//        ChastePoint<2> sproutBaseLocation(4,5,0);
//        ChastePoint<2> sproutTipLocation(5,5,0);
//
//        boost::shared_ptr<Vessel<2,int> > newSprout(new Vessel<2,int>());
//
//        newSprout = p_vessel_network->FormSprout(p_vessel_network->GetVessel(ChastePoint<2>(4,5,0),0), sproutBaseLocation, sproutTipLocation);
//
//        string testoutputdirectory = "VesselNetworkTests/";
//        string testsubdirectory = "FormSproutOnTipCell/";
//        string directory;
//        directory.append(testoutputdirectory);
//        directory.append(testsubdirectory);
//
//        OutputFileHandler output_file_handler(directory, true);
//
//        string initialStateFilename;
//        initialStateFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        initialStateFilename.append("VesselNetworkDataInitial.vtk");
//
//        string initialDataSummaryFilename;
//        initialDataSummaryFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        initialDataSummaryFilename.append("InitialVesselNetworkDataSummary.txt");
//        // print out and save some data about the vessel network before forrmation of sprout
//
//        p_vessel_network->SaveVasculatureDataToFile(initialStateFilename);
//
//        PrintVascularNetworkDataToFile(vesselNetwork, initialDataSummaryFilename);
//
//
//        // perform structural adaptation algorithm on new vessel network
//
//        p_vessel_network->accept(structuralAdaptationAlgorithm);
//
//        p_vessel_network->CheckForConsistency();
//
//        ChastePoint<2> sproutBaseLocation2(5,5,0);
//        ChastePoint<2> sproutTipLocation2(5,6,0);
//
//        boost::shared_ptr<Vessel<2,int> > newSprout2(new Vessel<2,int>());
//
//        TS_ASSERT_THROWS(p_vessel_network->FormSprout(newSprout, sproutBaseLocation2, sproutTipLocation2),Exception);
//
//
//        SimulationTime::Destroy();
//
//    }
//
//    void testFormSproutOnOutputNode()
//    {
//
//        // set up spatial mesh
//
//        int DomainSize_X = 10;
//        int DomainSize_Y = 10;
//        int DomainSize_Z = 1;
//        double dx = 20*pow(10.0,-6); // metres
//        const int dimensionality = 2;
//
//        /*
//         * Set up SimulationTime object.  This must be set up in order for structural adaptation algorithm
//         * to run.
//         */
//        SimulationTime* p_simulation_time = SimulationTime::Instance();
//        p_simulation_time->SetStartTime(0.0);
//        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(30, 1);
//
//        boost::shared_ptr<OnLatticeMesh<dimensionality,int> > spatialmesh(new OnLatticeMesh<dimensionality,int>(dx,DomainSize_X,DomainSize_Y,DomainSize_Z));
//        spatialmesh->SetVolumeFractionOccupiedByVessel(0.3);
//        // set up test vasular network
//
//        boost::shared_ptr<CaVascularNetwork<2> > vesselNetwork(new CaVascularNetwork<2>(spatialmesh));
//
//        boost::shared_ptr<Vessel<2,int> > vessel1(new Vessel<2,int>());
//        vessel1->SetNode1Location(ChastePoint<2>(4,0));
//        vessel1->SetNode2Location(ChastePoint<2>(4,5));
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,0), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,1), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,2), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,3), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,4), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,5), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel2(new Vessel<2,int>());
//        vessel2->SetNode1Location(ChastePoint<2>(4,5));
//        vessel2->SetNode2Location(ChastePoint<2>(4,9));
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,5), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,6), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,7), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,8), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,9), spatialmesh->GetSpatialMeshSize());
//
//
//        p_vessel_network->AddVessel(vessel1);
//        p_vessel_network->AddVessel(vessel2);
//
//
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//
//        p_vessel_network->SetOutputNode(ChastePoint<2>(4,9));
//        p_vessel_network->SetInputNode(ChastePoint<2>(4,0));
//
//        // test initial state of vascular network
//
//        // check locations of nodes
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,5)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,9)) == 1);
//
//
//        // check number of vessels attached to each node
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,0))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,5))->GetNumberOfSegments() == 2);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,9))->GetNumberOfSegments() == 1);
//
//
//        // check active tips
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode2() == false);
//
//
//        // test number of vessels and nodes in network
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfNodes(),3);
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfVessels(),2);
//
//        // test that vessels are in network
//
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel1));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel2));
//
//
//        // test that nodes are in network
//
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetEndNode()));
//
//
//        p_vessel_network->CheckForConsistency();
//
//        // perform structural adaptation algorithm on network
//
//        boost::shared_ptr<SimpleStructuralAdaptationAlgorithm<2,int> > structuralAdaptationAlgorithm(new SimpleStructuralAdaptationAlgorithm<2,int>());
//        // override default calculations
//        boost::shared_ptr<ConstantHaematocritCalculation<dimensionality,int> > haematocritcalculation(new ConstantHaematocritCalculation<dimensionality,int>());
//        structuralAdaptationAlgorithm->SetHaematocritCalculation(haematocritcalculation);
//
//        p_vessel_network->accept(structuralAdaptationAlgorithm);
//
//        // form sprout
//
//        ChastePoint<2> sproutBaseLocation(4,9,0);
//        ChastePoint<2> sproutTipLocation(5,9,0);
//
//        TS_ASSERT_THROWS(p_vessel_network->FormSprout(p_vessel_network->GetVessel(sproutBaseLocation,0), sproutBaseLocation, sproutTipLocation),Exception);
//
//
//        SimulationTime::Destroy();
//
//    }
//
//    void testFormSproutOnInputNode()
//    {
//
//        // set up spatial mesh
//
//        int DomainSize_X = 10;
//        int DomainSize_Y = 10;
//        int DomainSize_Z = 1;
//        double dx = 20*pow(10.0,-6); // metres
//        const int dimensionality = 2;
//
//        /*
//         * Set up SimulationTime object.  This must be set up in order for structural adaptation algorithm
//         * to run.
//         */
//        SimulationTime* p_simulation_time = SimulationTime::Instance();
//        p_simulation_time->SetStartTime(0.0);
//        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(30, 1);
//
//        boost::shared_ptr<OnLatticeMesh<dimensionality,int> > spatialmesh(new OnLatticeMesh<dimensionality,int>(dx,DomainSize_X,DomainSize_Y,DomainSize_Z));
//        spatialmesh->SetVolumeFractionOccupiedByVessel(0.3);
//        // set up test vasular network
//
//        boost::shared_ptr<CaVascularNetwork<2> > vesselNetwork(new CaVascularNetwork<2>(spatialmesh));
//
//        boost::shared_ptr<Vessel<2,int> > vessel1(new Vessel<2,int>());
//        vessel1->SetNode1Location(ChastePoint<2>(4,0));
//        vessel1->SetNode2Location(ChastePoint<2>(4,5));
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,0), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,1), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,2), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,3), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,4), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,5), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel2(new Vessel<2,int>());
//        vessel2->SetNode1Location(ChastePoint<2>(4,5));
//        vessel2->SetNode2Location(ChastePoint<2>(4,9));
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,5), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,6), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,7), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,8), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,9), spatialmesh->GetSpatialMeshSize());
//
//
//        p_vessel_network->AddVessel(vessel1);
//        p_vessel_network->AddVessel(vessel2);
//
//
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//
//        p_vessel_network->SetOutputNode(ChastePoint<2>(4,9));
//        p_vessel_network->SetInputNode(ChastePoint<2>(4,0));
//
//        // test initial state of vascular network
//
//        // check locations of nodes
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,5)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,9)) == 1);
//
//
//        // check number of vessels attached to each node
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,0))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,5))->GetNumberOfSegments() == 2);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,9))->GetNumberOfSegments() == 1);
//
//
//        // check active tips
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode2() == false);
//
//
//        // test number of vessels and nodes in network
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfNodes(),3);
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfVessels(),2);
//
//        // test that vessels are in network
//
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel1));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel2));
//
//
//        // test that nodes are in network
//
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetEndNode()));
//
//
//        p_vessel_network->CheckForConsistency();
//
//        // perform structural adaptation algorithm on network
//
//        boost::shared_ptr<SimpleStructuralAdaptationAlgorithm<2,int> > structuralAdaptationAlgorithm(new SimpleStructuralAdaptationAlgorithm<2,int>());
//        // override default calculations
//        boost::shared_ptr<ConstantHaematocritCalculation<dimensionality,int> > haematocritcalculation(new ConstantHaematocritCalculation<dimensionality,int>());
//        structuralAdaptationAlgorithm->SetHaematocritCalculation(haematocritcalculation);
//
//        p_vessel_network->accept(structuralAdaptationAlgorithm);
//
//        // form sprout
//
//        ChastePoint<2> sproutBaseLocation(4,0,0);
//        ChastePoint<2> sproutTipLocation(5,0,0);
//
//        TS_ASSERT_THROWS(p_vessel_network->FormSprout(p_vessel_network->GetVessel(sproutBaseLocation,0), sproutBaseLocation, sproutTipLocation),Exception);
//
//
//        SimulationTime::Destroy();
//
//    }
//
//    void testFormSproutOnNodeWhereThreeVesselsIntersect()
//    {
//
//        // set up spatial mesh
//
//        int DomainSize_X = 10;
//        int DomainSize_Y = 10;
//        int DomainSize_Z = 1;
//        double dx = 20*pow(10.0,-6); // metres
//        const int dimensionality = 2;
//
//        /*
//         * Set up SimulationTime object.  This must be set up in order for structural adaptation algorithm
//         * to run.
//         */
//        SimulationTime* p_simulation_time = SimulationTime::Instance();
//        p_simulation_time->SetStartTime(0.0);
//        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(30, 1);
//
//        boost::shared_ptr<OnLatticeMesh<dimensionality,int> > spatialmesh(new OnLatticeMesh<dimensionality,int>(dx,DomainSize_X,DomainSize_Y,DomainSize_Z));
//        spatialmesh->SetVolumeFractionOccupiedByVessel(0.25);
//        // set up test vasular network
//
//        boost::shared_ptr<CaVascularNetwork<2> > vesselNetwork(new CaVascularNetwork<2>(spatialmesh));
//
//        boost::shared_ptr<Vessel<2,int> > vessel1(new Vessel<2,int>());
//        vessel1->SetNode1Location(ChastePoint<2>(4,0));
//        vessel1->SetNode2Location(ChastePoint<2>(4,5));
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,0), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,1), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,2), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,3), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,4), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,5), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel2(new Vessel<2,int>());
//        vessel2->SetNode1Location(ChastePoint<2>(4,5));
//        vessel2->SetNode2Location(ChastePoint<2>(4,9));
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,5), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,6), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,7), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,8), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,9), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel3(new Vessel<2,int>());
//        vessel3->SetNode1Location(ChastePoint<2>(4,5));
//        vessel3->SetNode2Location(ChastePoint<2>(0,5));
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,5), spatialmesh->GetSpatialMeshSize());
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(3,5), spatialmesh->GetSpatialMeshSize());
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(2,5), spatialmesh->GetSpatialMeshSize());
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,5), spatialmesh->GetSpatialMeshSize());
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(0,5), spatialmesh->GetSpatialMeshSize());
//
//        p_vessel_network->AddVessel(vessel1);
//        p_vessel_network->AddVessel(vessel2);
//        p_vessel_network->AddVessel(vessel3);
//
//
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel3->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel3->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        p_vessel_network->SetOutputNode(ChastePoint<2>(4,9));
//        p_vessel_network->SetInputNode(ChastePoint<2>(4,0));
//        p_vessel_network->SetInputNode(ChastePoint<2>(0,5));
//
//        // test initial state of vascular network
//
//        // check locations of nodes
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,5)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,9)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(0,5)) == 1);
//
//        // check number of vessels attached to each node
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,0))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,5))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,9))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(0,5))->GetNumberOfSegments() == 1);
//
//
//        // check active tips
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel3->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel3->ActiveTipCellLocatedAtNode2() == false);
//
//
//        // test number of vessels and nodes in network
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfNodes(),4);
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfVessels(),3);
//
//        // test that vessels are in network
//
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel1));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel2));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel3));
//
//        // test that nodes are in network
//
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetEndNode()));
//
//
//        p_vessel_network->CheckForConsistency();
//
//        // perform structural adaptation algorithm on network
//
//        boost::shared_ptr<SimpleStructuralAdaptationAlgorithm<2,int> > structuralAdaptationAlgorithm(new SimpleStructuralAdaptationAlgorithm<2,int>());
//        // override default calculations
//        boost::shared_ptr<ConstantHaematocritCalculation<dimensionality,int> > haematocritcalculation(new ConstantHaematocritCalculation<dimensionality,int>());
//        structuralAdaptationAlgorithm->SetHaematocritCalculation(haematocritcalculation);
//
//        p_vessel_network->accept(structuralAdaptationAlgorithm);
//
//        string testoutputdirectory = "VesselNetworkTests/";
//        string testsubdirectory = "FormSproutOnNodeWhereThreeVesselsIntersect/";
//        string directory;
//        directory.append(testoutputdirectory);
//        directory.append(testsubdirectory);
//
//        OutputFileHandler output_file_handler(directory, true);
//
//        string initialStateFilename;
//        initialStateFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        initialStateFilename.append("VesselNetworkDataInitial.vtk");
//
//        string initialDataSummaryFilename;
//        initialDataSummaryFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        initialDataSummaryFilename.append("InitialVesselNetworkDataSummary.txt");
//
//        string finalStateFilename;
//        finalStateFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        finalStateFilename.append("VesselNetworkDataFinal.vtk");
//
//        string finalDataSummaryFilename;
//        finalDataSummaryFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        finalDataSummaryFilename.append("FinalVesselNetworkDataSummary.txt");
//
//        // print out and save some data about the vessel network before forrmation of sprout
//
//        p_vessel_network->SaveVasculatureDataToFile(initialStateFilename);
//
//        PrintVascularNetworkDataToFile(vesselNetwork, initialDataSummaryFilename);
//
//        // form sprout
//
//        ChastePoint<2> sproutBaseLocation(4,5,0);
//        ChastePoint<2> sproutTipLocation(5,5,0);
//
//        boost::shared_ptr<Vessel<2,int> > newSprout(new Vessel<2,int>());
//
//        newSprout = p_vessel_network->FormSprout(p_vessel_network->GetVessel(ChastePoint<2>(4,5,0),0), sproutBaseLocation, sproutTipLocation);
//
//        // perform structural adaptation algorithm on new vessel network
//
//        p_vessel_network->accept(structuralAdaptationAlgorithm);
//
//        // print out and save some data about the vessel network after formation of sprout
//
//        p_vessel_network->SaveVasculatureDataToFile(finalStateFilename);
//
//        PrintVascularNetworkDataToFile(vesselNetwork, finalDataSummaryFilename);
//
//        // test final state of vascular network
//
//        // check locations of nodes
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,5)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(5,5)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,9)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(0,5)) == 1);
//
//        // check number of vessels attached to each node
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,0))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,5))->GetNumberOfSegments() == 4);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(5,5))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,9))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(0,5))->GetNumberOfSegments() == 1);
//
//
//        // check active tips
//        for (int i = 0; i < p_vessel_network->GetNumberOfVessels(); i++)
//        {
//            if (newSprout == p_vessel_network->GetVessel(i))
//            {
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode1() == false);
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode2() == true);
//            }
//            else
//            {
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode1() == false);
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode2() == false);
//            }
//        }
//
//        // test number of vessels and nodes in network
//
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfNodes(),5);
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfVessels(),4);
//
//        // test that vessels are in network
//
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel1));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel2));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel3));
//
//        // test that nodes are in network
//
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetEndNode()));
//
//        // check for consistency between node and vessel circular references
//
//        for (int i = 0; i < p_vessel_network->GetNumberOfVessels(); i++)
//        {
//            TS_ASSERT(p_vessel_network->GetVessel(i)->GetStartNode()->IsAttachedToVessel(p_vessel_network->GetVessel(i)));
//            TS_ASSERT(p_vessel_network->GetVessel(i)->GetEndNode()->IsAttachedToVessel(p_vessel_network->GetVessel(i)));
//            TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel_network->GetVessel(i)->GetStartNode()));
//            TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel_network->GetVessel(i)->GetEndNode()));
//        }
//
//        for (int i = 0; i < p_vessel_network->GetNumberOfNodes(); i++)
//        {
//            for (int j = 0; j < p_vessel_network->GetNode(i)->GetNumberOfSegments(); j++)
//            {
//                TS_ASSERT(p_vessel_network->GetNode(i)->GetAdjoiningVessel(j)->IsAttachedToNode(p_vessel_network->GetNode(i)));
//                TS_ASSERT(p_vessel_network->VesselIsInNetwork(p_vessel_network->GetNode(i)->GetAdjoiningVessel(j)));
//            }
//        }
//
//        p_vessel_network->CheckForConsistency();
//
//        SimulationTime::Destroy();
//
//
//    }
//
//    void testFormSproutOnNodeWhereThreeVesselsIntersect2()
//    {
//
//        // set up spatial mesh
//
//        int DomainSize_X = 10;
//        int DomainSize_Y = 10;
//        int DomainSize_Z = 1;
//        double dx = 20*pow(10.0,-6); // metres
//        const int dimensionality = 2;
//
//        /*
//         * Set up SimulationTime object.  This must be set up in order for structural adaptation algorithm
//         * to run.
//         */
//        SimulationTime* p_simulation_time = SimulationTime::Instance();
//        p_simulation_time->SetStartTime(0.0);
//        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(30, 1);
//
//        boost::shared_ptr<OnLatticeMesh<dimensionality,int> > spatialmesh(new OnLatticeMesh<dimensionality,int>(dx,DomainSize_X,DomainSize_Y,DomainSize_Z));
//        spatialmesh->SetVolumeFractionOccupiedByVessel(0.25);
//        // set up test vasular network
//
//        boost::shared_ptr<CaVascularNetwork<2> > vesselNetwork(new CaVascularNetwork<2>(spatialmesh));
//
//        boost::shared_ptr<Vessel<2,int> > vessel1(new Vessel<2,int>());
//        vessel1->SetNode1Location(ChastePoint<2>(4,0));
//        vessel1->SetNode2Location(ChastePoint<2>(4,5));
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,0), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,1), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,2), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,3), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,4), spatialmesh->GetSpatialMeshSize());
//        vessel1->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,5), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel2(new Vessel<2,int>());
//        vessel2->SetNode1Location(ChastePoint<2>(4,5));
//        vessel2->SetNode2Location(ChastePoint<2>(4,9));
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,5), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,6), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,7), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,8), spatialmesh->GetSpatialMeshSize());
//        vessel2->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,9), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel3(new Vessel<2,int>());
//        vessel3->SetNode1Location(ChastePoint<2>(4,5));
//        vessel3->SetNode2Location(ChastePoint<2>(0,5));
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(4,5), spatialmesh->GetSpatialMeshSize());
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(3,5), spatialmesh->GetSpatialMeshSize());
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(2,5), spatialmesh->GetSpatialMeshSize());
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(1,5), spatialmesh->GetSpatialMeshSize());
//        vessel3->SetNextVesselSegmentCoordinate(ChastePoint<2>(0,5), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel4(new Vessel<2,int>());
//        vessel4->SetNode1Location(ChastePoint<2>(5,0));
//        vessel4->SetNode2Location(ChastePoint<2>(5,5));
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,0), spatialmesh->GetSpatialMeshSize());
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,1), spatialmesh->GetSpatialMeshSize());
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,2), spatialmesh->GetSpatialMeshSize());
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,3), spatialmesh->GetSpatialMeshSize());
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,4), spatialmesh->GetSpatialMeshSize());
//        vessel4->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,5), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel5(new Vessel<2,int>());
//        vessel5->SetNode1Location(ChastePoint<2>(5,5));
//        vessel5->SetNode2Location(ChastePoint<2>(5,9));
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,5), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,6), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,7), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,8), spatialmesh->GetSpatialMeshSize());
//        vessel5->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,9), spatialmesh->GetSpatialMeshSize());
//
//
//        boost::shared_ptr<Vessel<2,int> > vessel6(new Vessel<2,int>());
//        vessel6->SetNode1Location(ChastePoint<2>(5,5));
//        vessel6->SetNode2Location(ChastePoint<2>(9,5));
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(5,5), spatialmesh->GetSpatialMeshSize());
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(6,5), spatialmesh->GetSpatialMeshSize());
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(7,5), spatialmesh->GetSpatialMeshSize());
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(8,5), spatialmesh->GetSpatialMeshSize());
//        vessel6->SetNextVesselSegmentCoordinate(ChastePoint<2>(9,5), spatialmesh->GetSpatialMeshSize());
//
//
//        p_vessel_network->AddVessel(vessel1);
//        p_vessel_network->AddVessel(vessel2);
//        p_vessel_network->AddVessel(vessel3);
//        p_vessel_network->AddVessel(vessel4);
//        p_vessel_network->AddVessel(vessel5);
//        p_vessel_network->AddVessel(vessel6);
//
//
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel1->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel2->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel3->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel3->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel4->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel4->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel5->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel5->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        vessel6->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(20.0,"unitless"),3.0*pow(10.0,(-6)));
//        vessel6->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("VEGF",Concentration(50.0,"micromolar"),4.0*pow(10.0,(-6)));
//
//        p_vessel_network->SetOutputNode(ChastePoint<2>(4,9));
//        p_vessel_network->SetInputNode(ChastePoint<2>(4,0));
//        p_vessel_network->SetInputNode(ChastePoint<2>(0,5));
//        p_vessel_network->SetInputNode(ChastePoint<2>(9,5));
//
//        // test initial state of vascular network
//
//        // check locations of nodes
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,5)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,9)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(0,5)) == 1);
//
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(5,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(5,5)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(5,9)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(9,5)) == 1);
//
//        // check number of vessels attached to each node
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,0))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,5))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,9))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(0,5))->GetNumberOfSegments() == 1);
//
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(5,0))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(5,5))->GetNumberOfSegments() == 3);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(5,9))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(9,5))->GetNumberOfSegments() == 1);
//
//        // check active tips
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel1->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel2->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel3->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel3->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel4->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel4->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel5->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel5->ActiveTipCellLocatedAtNode2() == false);
//        TS_ASSERT(vessel6->ActiveTipCellLocatedAtNode1() == false);
//        TS_ASSERT(vessel6->ActiveTipCellLocatedAtNode2() == false);
//
//        // test number of vessels and nodes in network
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfNodes(),8);
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfVessels(),6);
//
//        // test that vessels are in network
//
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel1));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel2));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel3));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel4));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel5));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel6));
//
//        // test that nodes are in network
//
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel4->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel4->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel5->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel5->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel6->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel6->GetEndNode()));
//
//
//        p_vessel_network->CheckForConsistency();
//
//        // perform structural adaptation algorithm on network
//
//        boost::shared_ptr<SimpleStructuralAdaptationAlgorithm<2,int> > structuralAdaptationAlgorithm(new SimpleStructuralAdaptationAlgorithm<2,int>());
//        // override default calculations
//        boost::shared_ptr<ConstantHaematocritCalculation<dimensionality,int> > haematocritcalculation(new ConstantHaematocritCalculation<dimensionality,int>());
//        structuralAdaptationAlgorithm->SetHaematocritCalculation(haematocritcalculation);
//
//        p_vessel_network->accept(structuralAdaptationAlgorithm);
//
//        string testoutputdirectory = "VesselNetworkTests/";
//        string testsubdirectory = "FormSproutOnNodeWhereThreeVesselsIntersect2/";
//        string directory;
//        directory.append(testoutputdirectory);
//        directory.append(testsubdirectory);
//
//        OutputFileHandler output_file_handler(directory, true);
//
//        string initialStateFilename;
//        initialStateFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        initialStateFilename.append("VesselNetworkDataInitial.vtk");
//
//        string initialDataSummaryFilename;
//        initialDataSummaryFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        initialDataSummaryFilename.append("InitialVesselNetworkDataSummary.txt");
//
//        string finalStateFilename;
//        finalStateFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        finalStateFilename.append("VesselNetworkDataFinal.vtk");
//
//        string finalDataSummaryFilename;
//        finalDataSummaryFilename.append(output_file_handler.GetOutputDirectoryFullPath());
//        finalDataSummaryFilename.append("FinalVesselNetworkDataSummary.txt");
//
//        // print out and save some data about the vessel network before forrmation of sprout
//
//        p_vessel_network->SaveVasculatureDataToFile(initialStateFilename);
//
//        PrintVascularNetworkDataToFile(vesselNetwork, initialDataSummaryFilename);
//
//        // form sprout
//
//        ChastePoint<2> sproutBaseLocation(4,5,0);
//        ChastePoint<2> sproutTipLocation(5,5,0);
//
//        boost::shared_ptr<Vessel<2,int> > newSprout(new Vessel<2,int>());
//
//        newSprout = p_vessel_network->FormSprout(p_vessel_network->GetVessel(ChastePoint<2>(4,5,0),0), sproutBaseLocation, sproutTipLocation);
//
//        // perform structural adaptation algorithm on new vessel network
//
//        p_vessel_network->accept(structuralAdaptationAlgorithm);
//
//        // print out and save some data about the vessel network after formation of sprout
//
//        p_vessel_network->SaveVasculatureDataToFile(finalStateFilename);
//
//        PrintVascularNetworkDataToFile(vesselNetwork, finalDataSummaryFilename);
//
//        // test final state of vascular network
//
//        // check locations of nodes
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,5)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(5,5)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(4,9)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(0,5)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(5,0)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(5,9)) == 1);
//        TS_ASSERT(p_vessel_network->NumberOfNodesNearLocation(ChastePoint<2>(9,5)) == 1);
//
//        // check number of vessels attached to each node
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,0))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,5))->GetNumberOfSegments() == 4);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(5,5))->GetNumberOfSegments() == 4);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(4,9))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(0,5))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(9,5))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(5,0))->GetNumberOfSegments() == 1);
//        TS_ASSERT(p_vessel_network->GetNearestNode(ChastePoint<2>(5,9))->GetNumberOfSegments() == 1);
//
//
//        // check active tips
//        for (int i = 0; i < p_vessel_network->GetNumberOfVessels(); i++)
//        {
//            if (newSprout == p_vessel_network->GetVessel(i))
//            {
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode1() == false);
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode2() == false);
//            }
//            else
//            {
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode1() == false);
//                TS_ASSERT(p_vessel_network->GetVessel(i)->ActiveTipCellLocatedAtNode2() == false);
//            }
//        }
//
//        // test number of vessels and nodes in network
//
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfNodes(),8);
//        TS_ASSERT_EQUALS(p_vessel_network->GetNumberOfVessels(),7);
//
//        // test that vessels are in network
//
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel1));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel2));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel3));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel4));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel5));
//        TS_ASSERT(p_vessel_network->VesselIsInNetwork(vessel6));
//
//        // test that nodes are in network
//
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel1->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel2->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel3->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel4->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel4->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel5->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel5->GetEndNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel6->GetStartNode()));
//        TS_ASSERT(p_vessel_network->NodeIsInNetwork(vessel6->GetEndNode()));
//
//        // check for consistency between node and vessel circular references
//
//        for (int i = 0; i < p_vessel_network->GetNumberOfVessels(); i++)
//        {
//            TS_ASSERT(p_vessel_network->GetVessel(i)->GetStartNode()->IsAttachedToVessel(p_vessel_network->GetVessel(i)));
//            TS_ASSERT(p_vessel_network->GetVessel(i)->GetEndNode()->IsAttachedToVessel(p_vessel_network->GetVessel(i)));
//            TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel_network->GetVessel(i)->GetStartNode()));
//            TS_ASSERT(p_vessel_network->NodeIsInNetwork(p_vessel_network->GetVessel(i)->GetEndNode()));
//        }
//
//        for (int i = 0; i < p_vessel_network->GetNumberOfNodes(); i++)
//        {
//            for (int j = 0; j < p_vessel_network->GetNode(i)->GetNumberOfSegments(); j++)
//            {
//                TS_ASSERT(p_vessel_network->GetNode(i)->GetAdjoiningVessel(j)->IsAttachedToNode(p_vessel_network->GetNode(i)));
//                TS_ASSERT(p_vessel_network->VesselIsInNetwork(p_vessel_network->GetNode(i)->GetAdjoiningVessel(j)));
//            }
//        }
//
//        p_vessel_network->CheckForConsistency();
//
//        SimulationTime::Destroy();
//
//
//    }


};

#endif /*TESTVESSELNETWORK_HPP_*/
