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

        // Do the divide
        vessel_network.DivideVessel(vessel_network.GetVessel(0), ChastePoint<3>(0.66, 0.0, 0.0));
        TS_ASSERT_EQUALS(vessel_network.GetNumberOfVessels(), 2u);
        TS_ASSERT_DELTA(vessel_network.GetVessel(0)->GetSegment(0)->GetNode(0)->GetLocation()[0], 0.0, 1.e-6);
        TS_ASSERT_DELTA(vessel_network.GetVessel(0)->GetSegment(0)->GetNode(1)->GetLocation()[0], 0.66, 1.e-6);
        TS_ASSERT_DELTA(vessel_network.GetVessel(1)->GetSegment(0)->GetNode(0)->GetLocation()[0], 0.66, 1.e-6);
        TS_ASSERT_DELTA(vessel_network.GetVessel(1)->GetSegment(0)->GetNode(1)->GetLocation()[0], 2.0, 1.e-6);
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
        TS_ASSERT(!vessel_network.NodeIsInNetwork(nodes[1]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[2]));
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
        TS_ASSERT(!vessel_network.NodeIsInNetwork(nodes[1]));
        TS_ASSERT(vessel_network.NodeIsInNetwork(nodes[2]));
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
        query_nodes.push_back(nodes[2]);
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
};

#endif /*TESTVESSELNETWORK_HPP_*/
