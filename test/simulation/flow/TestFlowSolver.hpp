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

#ifndef TESTFlowSolver_HPP_
#define TESTFlowSolver_HPP_

#include <cxxtest/TestSuite.h>
#include "VesselImpedanceCalculator.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VasculatureGenerator.hpp"
#include "FlowSolver.hpp"
#include "FakePetscSetup.hpp"
#include "UnitCollection.hpp"

class TestFlowSolver : public CxxTest::TestSuite
{

    typedef boost::shared_ptr<VesselNode<2> > NodePtr2;
    typedef boost::shared_ptr<VesselNode<3> > NodePtr3;
    typedef boost::shared_ptr<VesselSegment<2> > SegmentPtr2;
    typedef boost::shared_ptr<VesselSegment<3> > SegmentPtr3;
    typedef boost::shared_ptr<Vessel<2> > VesselPtr2;
    typedef boost::shared_ptr<Vessel<3> > VesselPtr3;

public:

    void TestFlowThroughSingleSegment() throw (Exception)
    {

        // Make some nodes
        std::vector<ChastePoint<3> > points;
        points.push_back(ChastePoint<3>(0, 0, 0));
        points.push_back(ChastePoint<3>(5, 0, 0));

        boost::shared_ptr<VesselNode<3> > pn1 = VesselNode<3>::Create(0, 0, 0);
        boost::shared_ptr<VesselNode<3> > pn2 = VesselNode<3>::Create(5, 0, 0);

        SegmentPtr3 p_segment(VesselSegment<3>::Create(pn1, pn2));
        VesselPtr3 p_vessel(Vessel<3>::Create(p_segment));

        // Generate the network
        boost::shared_ptr<VesselNetwork<3> > p_vascular_network = VesselNetwork<3>::Create();
        p_vascular_network->AddVessel(p_vessel);

        units::quantity<unit::flow_impedance> impedance = 1.e12 * unit::unit_flow_impedance;
        p_segment->GetFlowProperties()->SetDimensionalImpedance(impedance);
        p_vascular_network->SetSegmentProperties(p_segment);

        p_vessel->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
        p_vessel->GetStartNode()->GetFlowProperties()->SetPressure(3393);

        p_vessel->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
        p_vessel->GetEndNode()->GetFlowProperties()->SetPressure(1000.5);

        FlowSolver<3> solver;
        solver.SetVesselNetwork(p_vascular_network);
        solver.SetUp();
        solver.Solve();

        TS_ASSERT_DELTA(p_vessel->GetStartNode()->GetFlowProperties()->GetPressure(), 3393, 1e-6);
        TS_ASSERT_DELTA(p_vessel->GetEndNode()->GetFlowProperties()->GetPressure(), 1000.5, 1e-6);

        TS_ASSERT_DELTA(p_vessel->GetFlowProperties()->GetDimensionalFlowRate(p_vessel->GetSegments())/unit::unit_flow_rate, (3393 - 1000.5) *unit::unit_flow_impedance/ impedance, 1e-6);
        TS_ASSERT_DELTA(p_segment->GetFlowProperties()->GetDimensionalFlowRate()/unit::unit_flow_rate, (3393 - 1000.5)*unit::unit_flow_impedance / impedance, 1e-6);

        p_segment->GetFlowProperties()->SetDimensionalImpedance(-1.0*unit::unit_flow_impedance);
        TS_ASSERT_THROWS_THIS(solver.Update(), "Impedance should be a positive number.");
    }

    void TestFlowThroughSingleVesselWithMultipleSegments() throw (Exception)
    {
        // Make some nodes
        std::vector<NodePtr3> nodes;
        nodes.push_back(NodePtr3(VesselNode<3>::Create(1.0, 0, 0)));
        nodes.push_back(NodePtr3(VesselNode<3>::Create(2.0, 0, 0)));
        nodes.push_back(NodePtr3(VesselNode<3>::Create(3.0, 0, 0)));
        nodes.push_back(NodePtr3(VesselNode<3>::Create(4.0, 0, 0)));
        nodes.push_back(NodePtr3(VesselNode<3>::Create(5.0, 0, 0)));
        SegmentPtr3 p_segment1(VesselSegment<3>::Create(nodes[0], nodes[1]));
        SegmentPtr3 p_segment2(VesselSegment<3>::Create(nodes[1], nodes[2]));
        SegmentPtr3 p_segment3(VesselSegment<3>::Create(nodes[2], nodes[3]));
        SegmentPtr3 p_segment4(VesselSegment<3>::Create(nodes[3], nodes[4]));
        std::vector<SegmentPtr3> segments;
        segments.push_back(p_segment1);
        segments.push_back(p_segment2);
        segments.push_back(p_segment3);
        segments.push_back(p_segment4);

        VesselPtr3 p_vessel(Vessel<3>::Create(segments));

        // Generate the network
        boost::shared_ptr<VesselNetwork<3> > p_vascular_network(new VesselNetwork<3>());
        p_vascular_network->AddVessel(p_vessel);
        double impedance = 1.e14;
        p_segment1->GetFlowProperties()->SetImpedance(1.e14);
        p_vascular_network->SetSegmentProperties(p_segment1);

        p_vessel->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
        p_vessel->GetStartNode()->GetFlowProperties()->SetPressure(3393);

        p_vessel->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
        p_vessel->GetEndNode()->GetFlowProperties()->SetPressure(1000.5);

        FlowSolver<3> solver;
        solver.SetVesselNetwork(p_vascular_network);
        solver.SetUp();
        solver.Solve();

        for (unsigned i = 0; i < nodes.size(); i++)
        {
            TS_ASSERT_DELTA(nodes[i]->GetFlowProperties()->GetPressure(),
                            3393 - (3393 - 1000.5) * i / (nodes.size() - 1), 1e-6);
        }

        TS_ASSERT_DELTA(p_vessel->GetStartNode()->GetFlowProperties()->GetPressure(), 3393, 1e-6);
        TS_ASSERT_DELTA(p_vessel->GetEndNode()->GetFlowProperties()->GetPressure(), 1000.5, 1e-6);

        TS_ASSERT_DELTA(p_vessel->GetFlowProperties()->GetFlowRate(p_vessel->GetSegments()), (3393 - 1000.5) / (segments.size() * impedance), 1e-6);

        for (unsigned i = 0; i < segments.size(); i++)
        {
            TS_ASSERT_DELTA(segments[i]->GetFlowProperties()->GetFlowRate(),
                            (3393 - 1000.5) / (segments.size() * impedance), 1e-6);
        }

    }

    void TestFlowThroughMultipleVessels() throw (Exception)
    {
        // Make some nodes
        std::vector<NodePtr3> nodes;
        nodes.push_back(NodePtr3(VesselNode<3>::Create(1.0, 0, 0)));
        nodes.push_back(NodePtr3(VesselNode<3>::Create(2.0, 0, 0)));
        nodes.push_back(NodePtr3(VesselNode<3>::Create(3.0, 0, 0)));

        SegmentPtr3 p_segment1(VesselSegment<3>::Create(nodes[0], nodes[1]));
        SegmentPtr3 p_segment2(VesselSegment<3>::Create(nodes[1], nodes[2]));

        VesselPtr3 p_vessel1(Vessel<3>::Create(p_segment1));
        VesselPtr3 p_vessel2(Vessel<3>::Create(p_segment2));

        // Generate the network
        boost::shared_ptr<VesselNetwork<3> > p_vascular_network(new VesselNetwork<3>());

        p_vascular_network->AddVessel(p_vessel1);
        p_vascular_network->AddVessel(p_vessel2);

        double impedance = 1.e14;
        p_segment1->GetFlowProperties()->SetImpedance(impedance);
        p_vascular_network->SetSegmentProperties(p_segment1);

        nodes[0]->GetFlowProperties()->SetIsInputNode(true);
        nodes[0]->GetFlowProperties()->SetPressure(3393);
        nodes[2]->GetFlowProperties()->SetIsOutputNode(true);
        nodes[2]->GetFlowProperties()->SetPressure(1000.5);

        FlowSolver<3> solver;
        solver.SetVesselNetwork(p_vascular_network);
        solver.SetUp();
        solver.Solve();

        TS_ASSERT_DELTA(p_vessel1->GetStartNode()->GetFlowProperties()->GetPressure(), 3393, 1e-6);
        TS_ASSERT_DELTA(p_vessel2->GetEndNode()->GetFlowProperties()->GetPressure(), 1000.5, 1e-6);
        TS_ASSERT_DELTA(nodes[1]->GetFlowProperties()->GetPressure(), (3393 + 1000.5) / 2.0, 1e-6);
        TS_ASSERT_DELTA(p_vessel1->GetFlowProperties()->GetFlowRate(p_vessel1->GetSegments()), (3393 - 1000.5) / (2.0 * impedance), 1e-6);

    }

    void TestFlowThroughBifurcation() throw (Exception)
    {

        // Make some nodes
        std::vector<NodePtr3> nodes;
        nodes.push_back(NodePtr3(VesselNode<3>::Create(0.0, 0, 0)));
        nodes.push_back(NodePtr3(VesselNode<3>::Create(0.0, 1, 0)));
        nodes.push_back(NodePtr3(VesselNode<3>::Create(1.0, 0.5, 0)));
        nodes.push_back(NodePtr3(VesselNode<3>::Create(1.0, 1, 0)));
        nodes.push_back(NodePtr3(VesselNode<3>::Create(2.0, 1, 0)));
        nodes.push_back(NodePtr3(VesselNode<3>::Create(3.0, 1, 0)));
        nodes.push_back(NodePtr3(VesselNode<3>::Create(4.0, 1, 0)));

        SegmentPtr3 p_segment1(VesselSegment<3>::Create(nodes[0], nodes[2]));
        SegmentPtr3 p_segment2(VesselSegment<3>::Create(nodes[1], nodes[2]));
        SegmentPtr3 p_segment3(VesselSegment<3>::Create(nodes[3], nodes[2]));
        SegmentPtr3 p_segment4(VesselSegment<3>::Create(nodes[4], nodes[3]));
        SegmentPtr3 p_segment5(VesselSegment<3>::Create(nodes[5], nodes[4]));
        SegmentPtr3 p_segment6(VesselSegment<3>::Create(nodes[6], nodes[5]));

        VesselPtr3 p_vessel1(Vessel<3>::Create(p_segment1));
        VesselPtr3 p_vessel2(Vessel<3>::Create(p_segment2));
        VesselPtr3 p_vessel3(Vessel<3>::Create(p_segment3));
        VesselPtr3 p_vessel4(Vessel<3>::Create(p_segment4));
        VesselPtr3 p_vessel5(Vessel<3>::Create(p_segment5));
        VesselPtr3 p_vessel6(Vessel<3>::Create(p_segment6));

        std::vector<VesselPtr3> vessels;
        vessels.push_back(p_vessel1); // lower input vessel
        vessels.push_back(p_vessel2); // upper input vessel
        vessels.push_back(p_vessel3);
        vessels.push_back(p_vessel4);
        vessels.push_back(p_vessel5);
        vessels.push_back(p_vessel6);

        // Generate the network
        boost::shared_ptr<VesselNetwork<3> > p_vascular_network(new VesselNetwork<3>());

        p_vascular_network->AddVessels(vessels);

        double impedance = 1.e14;
        p_segment1->GetFlowProperties()->SetImpedance(impedance);
        p_vascular_network->SetSegmentProperties(p_segment1);

        nodes[0]->GetFlowProperties()->SetIsInputNode(true);
        nodes[0]->GetFlowProperties()->SetPressure(3393);

        nodes[1]->GetFlowProperties()->SetIsInputNode(true);
        nodes[1]->GetFlowProperties()->SetPressure(3393);

        nodes[6]->GetFlowProperties()->SetIsOutputNode(true);
        nodes[6]->GetFlowProperties()->SetPressure(1000.5);

        FlowSolver<3> solver;
        solver.SetVesselNetwork(p_vascular_network);
        solver.SetUp();
        solver.Solve();

        TS_ASSERT_DELTA(nodes[0]->GetFlowProperties()->GetPressure(), 3393, 1e-6);
        TS_ASSERT_DELTA(nodes[1]->GetFlowProperties()->GetPressure(), 3393, 1e-6);
        TS_ASSERT_DELTA(nodes[2]->GetFlowProperties()->GetPressure(), (2.0 * 3393.0 / 10.0 + 1000.5 / 40.0) / (1.0 / 40.0 + 2.0 / 10.0), 1e-6);
        TS_ASSERT_DELTA(nodes[6]->GetFlowProperties()->GetPressure(), 1000.5, 1e-6);
        TS_ASSERT_DELTA(vessels[0]->GetFlowProperties()->GetFlowRate(vessels[0]->GetSegments()), (3393 - nodes[2]->GetFlowProperties()->GetPressure())/ impedance, 1e-6);
        TS_ASSERT_DELTA(vessels[1]->GetFlowProperties()->GetFlowRate(vessels[1]->GetSegments()), (3393 - nodes[2]->GetFlowProperties()->GetPressure()) / impedance, 1e-6);
        TS_ASSERT_DELTA(vessels[5]->GetFlowProperties()->GetFlowRate(vessels[5]->GetSegments()), -(nodes[2]->GetFlowProperties()->GetPressure() - 1000.5) / (4.0 * impedance), 1e-6);

        TS_ASSERT_DELTA(p_segment1->GetFlowProperties()->GetFlowRate(), (3393 - nodes[2]->GetFlowProperties()->GetPressure()) / impedance, 1e-6);
        TS_ASSERT_DELTA(p_segment2->GetFlowProperties()->GetFlowRate(),(3393 - nodes[2]->GetFlowProperties()->GetPressure()) / impedance, 1e-6);
        TS_ASSERT_DELTA(p_segment6->GetFlowProperties()->GetFlowRate(),-(nodes[2]->GetFlowProperties()->GetPressure() - 1000.5)/ (4.0 * impedance), 1e-6);

        double kirchoff_residual = vessels[0]->GetFlowProperties()->GetFlowRate(vessels[0]->GetSegments()) +
                vessels[1]->GetFlowProperties()->GetFlowRate(vessels[1]->GetSegments()) +
                vessels[5]->GetFlowProperties()->GetFlowRate(vessels[5]->GetSegments());

        TS_ASSERT_DELTA(kirchoff_residual, 0, 1e-6);

    }

    void TestFlowThroughBifurcationHavingSwappedNodeLabels() throw (Exception)
    {
        std::vector<NodePtr3> nodes;
        nodes.push_back(NodePtr3(VesselNode<3>::Create(0.0, 0, 0)));
        nodes.push_back(NodePtr3(VesselNode<3>::Create(0.0, 1, 0)));
        nodes.push_back(NodePtr3(VesselNode<3>::Create(1.0, 0.5, 0)));
        nodes.push_back(NodePtr3(VesselNode<3>::Create(1.0, 1, 0)));

        SegmentPtr3 p_segment1(VesselSegment<3>::Create(nodes[0], nodes[2]));
        SegmentPtr3 p_segment2(VesselSegment<3>::Create(nodes[1], nodes[2]));
        SegmentPtr3 p_segment3(VesselSegment<3>::Create(nodes[2], nodes[3]));

        VesselPtr3 p_vessel1(Vessel<3>::Create(p_segment1));
        VesselPtr3 p_vessel2(Vessel<3>::Create(p_segment2));
        VesselPtr3 p_vessel3(Vessel<3>::Create(p_segment3));

        std::vector<VesselPtr3> vessels;
        vessels.push_back(p_vessel1); // lower input vessel
        vessels.push_back(p_vessel2); // upper input vessel
        vessels.push_back(p_vessel3); // output vessel

        // Generate the network
        boost::shared_ptr<VesselNetwork<3> > p_vascular_network(new VesselNetwork<3>());

        p_vascular_network->AddVessels(vessels);

        double impedance = 1.e14;
        p_segment1->GetFlowProperties()->SetImpedance(impedance);
        p_vascular_network->SetSegmentProperties(p_segment1);

        nodes[0]->GetFlowProperties()->SetIsInputNode(true);
        nodes[0]->GetFlowProperties()->SetPressure(3393);

        nodes[1]->GetFlowProperties()->SetIsInputNode(true);
        nodes[1]->GetFlowProperties()->SetPressure(3393);

        nodes[3]->GetFlowProperties()->SetIsOutputNode(true);
        nodes[3]->GetFlowProperties()->SetPressure(1000.5);

        FlowSolver<3> solver;
        solver.SetVesselNetwork(p_vascular_network);
        solver.SetUp();
        solver.Solve();
        TS_ASSERT_DELTA(nodes[0]->GetFlowProperties()->GetPressure(), 3393, 1e-6);
        TS_ASSERT_DELTA(nodes[1]->GetFlowProperties()->GetPressure(), 3393, 1e-6);
        TS_ASSERT_DELTA(nodes[2]->GetFlowProperties()->GetPressure(), (2 * 3393 + 1000.5) / 3, 1e-6);
        TS_ASSERT_DELTA(nodes[3]->GetFlowProperties()->GetPressure(), 1000.5, 1e-6);

        TS_ASSERT_DELTA(vessels[0]->GetFlowProperties()->GetFlowRate(vessels[0]->GetSegments()), (3393 - nodes[2]->GetFlowProperties()->GetPressure())/ impedance,
                        1e-6);
        TS_ASSERT_DELTA(vessels[1]->GetFlowProperties()->GetFlowRate(vessels[1]->GetSegments()), (3393 - nodes[2]->GetFlowProperties()->GetPressure()) / impedance,
                        1e-6);
        TS_ASSERT_DELTA(vessels[2]->GetFlowProperties()->GetFlowRate(vessels[2]->GetSegments()), (nodes[2]->GetFlowProperties()->GetPressure() - 1000.5) / impedance,
                        1e-6);

        TS_ASSERT_DELTA(p_segment1->GetFlowProperties()->GetFlowRate(),
                        (3393 - nodes[2]->GetFlowProperties()->GetPressure())/ impedance, 1e-6);
        TS_ASSERT_DELTA(p_segment2->GetFlowProperties()->GetFlowRate(),
                        (3393 - nodes[2]->GetFlowProperties()->GetPressure())/ impedance, 1e-6);
        TS_ASSERT_DELTA(p_segment3->GetFlowProperties()->GetFlowRate(),
                        (nodes[2]->GetFlowProperties()->GetPressure() - 1000.5) / impedance, 1e-6);

        double kirchoff_residual = vessels[0]->GetFlowProperties()->GetFlowRate(vessels[0]->GetSegments()) +
                vessels[1]->GetFlowProperties()->GetFlowRate(vessels[1]->GetSegments()) -
                vessels[2]->GetFlowProperties()->GetFlowRate(vessels[2]->GetSegments());

        TS_ASSERT_DELTA(kirchoff_residual, 0, 1e-6);

    }

    void TestFlowThroughHexagonalNetwork() throw (Exception)
    {
        // Specify the network dimensions
        double vessel_length = 80.0;

        // Generate the network
        VasculatureGenerator<2> vascular_network_generator;
        boost::shared_ptr<VesselNetwork<2> > vascular_network = vascular_network_generator.GenerateHexagonalNetwork(
                1000, 1000, vessel_length);

        // Make some nodes
        std::vector<NodePtr2> nodes;
        nodes.push_back(NodePtr2(VesselNode<2>::Create(0.0, 0)));
        nodes.push_back(NodePtr2(VesselNode<2>::Create(0.0, 1)));
        SegmentPtr2 p_segment1(VesselSegment<2>::Create(nodes[0], nodes[1]));

        double impedance = 0.001;
        p_segment1->GetFlowProperties()->SetImpedance(impedance);
        vascular_network->SetSegmentProperties(p_segment1);

        std::vector<std::pair<double, double> > extents = vascular_network->GetExtents();
        double y_middle = (extents[1].first + extents[1].second) / 2.0;

        std::vector<boost::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
        std::vector<boost::shared_ptr<Vessel<2> > > vessels = vascular_network->GetVessels();
        for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
        {
            if ((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
            {
                if ((*vessel_iterator)->GetStartNode()->rGetLocation()[1] > y_middle)
                {
                    (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
                    (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3393);
                }
            }
            if ((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
            {
                if ((*vessel_iterator)->GetEndNode()->rGetLocation()[1] > y_middle)
                {
                    (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
                    (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3393);
                }
            }
            if ((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
            {
                if ((*vessel_iterator)->GetStartNode()->rGetLocation()[1] <= y_middle)
                {
                    (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
                    (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(1993);
                }
            }
            if ((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
            {
                if ((*vessel_iterator)->GetEndNode()->rGetLocation()[1] <= y_middle)
                {
                    (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
                    (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(1993);
                }
            }

        }

        FlowSolver<2> solver;
        solver.SetVesselNetwork(vascular_network);
        solver.SetUp();
        solver.Solve();

        // Write the network to file
        OutputFileHandler output_file_handler("TestFlowSolver", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexagonalVesselNetwork.vtp");
        vascular_network->Write(output_filename);
    }
//
//    void TestLoop() throw(Exception)
//    {
//        // Make a network
//        std::vector<boost::shared_ptr<VesselNode<3> > > bottom_nodes;
//        for(unsigned idx=0; idx<5; idx++)
//        {
//            bottom_nodes.push_back(VesselNode<3>::Create(double(idx)*10, 10.0, 0.0));
//        }
//        bottom_nodes[0]->GetFlowProperties()->SetIsInputNode(true);
//        bottom_nodes[0]->GetFlowProperties()->SetPressure(3000);
//        bottom_nodes[4]->GetFlowProperties()->SetIsOutputNode(true);
//        bottom_nodes[4]->GetFlowProperties()->SetPressure(1000);
//
//        std::vector<boost::shared_ptr<VesselNode<3> > > top_nodes;
//        for(unsigned idx=1; idx<3; idx+=1)
//        {
//            top_nodes.push_back(VesselNode<3>::Create(double(idx)*10, 20.0, 0.0));
//        }
//
//        boost::shared_ptr<Vessel<3> > p_vessel1 = Vessel<3>::Create(bottom_nodes[0], bottom_nodes[1]);
//        boost::shared_ptr<Vessel<3> > p_vessel2 = Vessel<3>::Create(bottom_nodes[1], bottom_nodes[2]);
//        boost::shared_ptr<Vessel<3> > p_vessel3 = Vessel<3>::Create(bottom_nodes[2], bottom_nodes[3]);
//        boost::shared_ptr<Vessel<3> > p_vessel7 = Vessel<3>::Create(bottom_nodes[3], bottom_nodes[4]);
//        boost::shared_ptr<Vessel<3> > p_vessel4 = Vessel<3>::Create(bottom_nodes[1], top_nodes[0]);
//        boost::shared_ptr<Vessel<3> > p_vessel5 = Vessel<3>::Create(bottom_nodes[2], top_nodes[1]);
//        boost::shared_ptr<Vessel<3> > p_vessel6 = Vessel<3>::Create(top_nodes[0], top_nodes[1]);
//
//        boost::shared_ptr<VesselNetwork<3> > p_network = VesselNetwork<3>::Create();
//        p_network->AddVessel(p_vessel1);
//        p_network->AddVessel(p_vessel2);
//        p_network->AddVessel(p_vessel3);
//        p_network->AddVessel(p_vessel4);
//        p_network->AddVessel(p_vessel5);
//        p_network->AddVessel(p_vessel6);
//        p_network->AddVessel(p_vessel7);
//        p_network->SetSegmentRadii(10.0);
//        std::vector<boost::shared_ptr<VesselSegment<3> > > segments = p_network->GetVesselSegments();
//        for(unsigned idx=0; idx<segments.size(); idx++)
//        {
//            segments[idx]->GetFlowProperties()->SetViscosity(1.e-3);
//        }
//
//        // Grow the vessel
//        PoiseuilleImpedanceCalculator<3> impedance_calculator;
//        impedance_calculator.Calculate(p_network);
//        FlowSolver<3> solver;
//        solver.SetVesselNetwork(p_network);
//        solver.SetUp();
//        solver.Solve();
//
//        // Write the network to file
//        OutputFileHandler output_file_handler("TestFlowSolver", false);
//        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("LoopFlow.vtp");
//        p_network->Write(output_filename);
//    }
//
//    void TestSproutingWithFlow() throw(Exception)
//    {
//        // Make a network
//        std::vector<boost::shared_ptr<VesselNode<3> > > bottom_nodes;
//        for(unsigned idx=0; idx<6; idx++)
//        {
//            bottom_nodes.push_back(VesselNode<3>::Create(double(idx)*10, 10.0, 0.0));
//        }
//        bottom_nodes[0]->GetFlowProperties()->SetIsInputNode(true);
//        bottom_nodes[0]->GetFlowProperties()->SetPressure(3000);
//        bottom_nodes[5]->GetFlowProperties()->SetIsOutputNode(true);
//        bottom_nodes[5]->GetFlowProperties()->SetPressure(1000);
//
//        boost::shared_ptr<Vessel<3> > p_vessel1 = Vessel<3>::Create(bottom_nodes);
//        boost::shared_ptr<VesselNetwork<3> > p_network = VesselNetwork<3>::Create();
//        p_network->AddVessel(p_vessel1);
//        p_network->SetSegmentRadii(10.0);
//
//        for(unsigned idx=1; idx<4; idx+=1)
//        {
//            p_network->FormSprout(ChastePoint<3>(double(idx)*10, 10.0, 0.0), ChastePoint<3>(double(idx)*10, 20.0, 0.0));
//        }
//
//        p_network->FormSprout(ChastePoint<3>(10, 20.0, 0.0), ChastePoint<3>(20, 20.0, 0.0));
//        p_network->MergeCoincidentNodes();
//        p_network->UpdateSegments();
//        p_network->UpdateNodes();
//        p_network->UpdateVesselNodes();
//        std::vector<boost::shared_ptr<VesselSegment<3> > > segments = p_network->GetVesselSegments();
//        for(unsigned idx=0; idx<segments.size(); idx++)
//        {
//            segments[idx]->GetFlowProperties()->SetViscosity(1.e-3);
//        }
//
//        // Grow the vessel
//        PoiseuilleImpedanceCalculator<3> impedance_calculator;
//        impedance_calculator.Calculate(p_network);
//        FlowSolver<3> solver;
//        solver.SetVesselNetwork(p_network);
//        solver.SetUp();
//        solver.Solve();
//
//        // Write the network to file
//        OutputFileHandler output_file_handler("TestFlowSolver", false);
//        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("SproutingFlow.vtp");
//        p_network->Write(output_filename);
//    }
//
//    void TestFlowThroughHexagonalNetworkWithSprouting() throw (Exception)
//    {
//        // Specify the network dimensions
//        double vessel_length = 80.0;
//
//        // Generate the network
//        VasculatureGenerator<2> vascular_network_generator;
//        boost::shared_ptr<VesselNetwork<2> > vascular_network = vascular_network_generator.GenerateHexagonalNetwork(
//                1000, 1000, vessel_length);
//
//        // Make some nodes
//        std::vector<ChastePoint<2> > points;
//        points.push_back(ChastePoint<2>(0, 0, 0)); // input
//        points.push_back(ChastePoint<2>(0, 1.0, 0)); // input
//
//        std::vector<NodePtr2> nodes;
//        for (unsigned i = 0; i < points.size(); i++)
//        {
//            nodes.push_back(NodePtr2(VesselNode<2>::Create(points[i])));
//        }
//
//        SegmentPtr2 p_segment1(VesselSegment<2>::Create(nodes[0], nodes[1]));
//
//        double impedance = 1.e14;
//        p_segment1->GetFlowProperties()->SetImpedance(impedance);
//        vascular_network->SetSegmentProperties(p_segment1);
//
//        std::vector<std::pair<double, double> > extents = vascular_network->GetExtents();
//        double y_middle = (extents[1].first + extents[1].second) / 2.0;
//
//        std::vector<boost::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
//
//        std::vector<boost::shared_ptr<Vessel<2> > > vessels = vascular_network->GetVessels();
//
//        for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
//        {
//            if ((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
//            {
//                if ((*vessel_iterator)->GetStartNode()->GetLocation()[1] > y_middle)
//                {
//                    (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
//                    (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3393);
//                }
//            }
//            if ((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
//            {
//                if ((*vessel_iterator)->GetEndNode()->GetLocation()[1] > y_middle)
//                {
//                    (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
//                    (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3393);
//                }
//            }
//            if ((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
//            {
//                if ((*vessel_iterator)->GetStartNode()->GetLocation()[1] <= y_middle)
//                {
//                    (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
//                    (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(1993);
//                }
//            }
//            if ((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
//            {
//                if ((*vessel_iterator)->GetEndNode()->GetLocation()[1] <= y_middle)
//                {
//                    (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
//                    (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(1993);
//                }
//            }
//        }
//
//        FlowSolver<2> solver;
//        solver.SetVesselNetwork(vascular_network);
//        solver.SetUp();
//        solver.Solve();
//
//        // Write the network to file
//        OutputFileHandler output_file_handler("TestFlowSolver", false);
//        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexagonalVesselNetwork.vtp");
//        vascular_network->Write(output_filename);
//    }
};

#endif /*TESTFlowSolver_HPP_*/
