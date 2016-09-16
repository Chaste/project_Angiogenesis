/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTALARCONHAEMATOCRITSOLVER_HPP
#define TESTALARCONHAEMATOCRITSOLVER_HPP

#include <cxxtest/TestSuite.h>
#include "VesselImpedanceCalculator.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VesselNetworkGenerator.hpp"
#include "FlowSolver.hpp"
#include "SimulationTime.hpp"
#include "AlarconHaematocritSolver.hpp"
#include "UnitCollection.hpp"

#include "FakePetscSetup.hpp"

class TestAlarconHaematocritSolver : public CxxTest::TestSuite
{

public:

    void TestTwoVesselNetwork() throw(Exception)
    {
        boost::shared_ptr<VesselNode<2> > p_node1 = VesselNode<2>::Create(0.0, 0.0);
        boost::shared_ptr<VesselNode<2> > p_node2 = VesselNode<2>::Create(80, 0.0);
        boost::shared_ptr<VesselNode<2> > p_node3 = VesselNode<2>::Create(160, 0.0);
        p_node1->GetFlowProperties()->SetIsInputNode(true);

        boost::shared_ptr<VesselSegment<2> > p_segment1(VesselSegment<2>::Create(p_node1, p_node2));
        boost::shared_ptr<VesselSegment<2> > p_segment2(VesselSegment<2>::Create(p_node2, p_node3));

        boost::shared_ptr<Vessel<2> > p_vessel1(Vessel<2>::Create(p_segment1));
        boost::shared_ptr<Vessel<2> > p_vessel2(Vessel<2>::Create(p_segment2));
        p_segment1->GetFlowProperties()->SetFlowRate(1.0*unit::metre_cubed_per_second);
        p_segment2->GetFlowProperties()->SetFlowRate(2.0*unit::metre_cubed_per_second);

        boost::shared_ptr<VesselNetwork<2> > p_network = boost::shared_ptr<VesselNetwork<2> >(new VesselNetwork<2>);
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);

        boost::shared_ptr<AlarconHaematocritSolver<2> > p_haematocrit_calculator(new AlarconHaematocritSolver<2>());
        p_haematocrit_calculator->SetVesselNetwork(p_network);
        p_haematocrit_calculator->Calculate();

        TS_ASSERT_DELTA(p_vessel1->GetFlowProperties()->GetHaematocrit(p_vessel1->GetSegments()),0.45, 1e-6);
        TS_ASSERT_DELTA(p_vessel2->GetFlowProperties()->GetHaematocrit(p_vessel2->GetSegments()),0.45, 1e-6);
    }

    void TestBifurcationInflowNetwork() throw(Exception)
    {
        boost::shared_ptr<VesselNode<2> > p_node1 = VesselNode<2>::Create(0.0, 0.0);
        boost::shared_ptr<VesselNode<2> > p_node2 = VesselNode<2>::Create(80, 0.0);
        boost::shared_ptr<VesselNode<2> > p_node3 = VesselNode<2>::Create(160, 0.0);
        boost::shared_ptr<VesselNode<2> > p_node4 = VesselNode<2>::Create(200, 0.0);
        p_node1->GetFlowProperties()->SetIsInputNode(true);
        p_node2->GetFlowProperties()->SetIsInputNode(true);

        boost::shared_ptr<Vessel<2> > p_vessel1(Vessel<2>::Create(p_node1, p_node3));
        boost::shared_ptr<Vessel<2> > p_vessel2(Vessel<2>::Create(p_node2, p_node3));
        boost::shared_ptr<Vessel<2> > p_vessel3(Vessel<2>::Create(p_node3, p_node4));
        p_vessel1->GetFlowProperties()->SetFlowRate(1.0*unit::metre_cubed_per_second, p_vessel1->GetSegments());
        p_vessel2->GetFlowProperties()->SetFlowRate(1.0*unit::metre_cubed_per_second, p_vessel2->GetSegments());
        p_vessel3->GetFlowProperties()->SetFlowRate(1.0*unit::metre_cubed_per_second, p_vessel3->GetSegments() );

        boost::shared_ptr<VesselNetwork<2> > p_network = boost::shared_ptr<VesselNetwork<2> >(new VesselNetwork<2>);
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);
        p_network->AddVessel(p_vessel3);

        boost::shared_ptr<AlarconHaematocritSolver<2> > p_haematocrit_calculator(new AlarconHaematocritSolver<2>());
        p_haematocrit_calculator->SetVesselNetwork(p_network);
        p_haematocrit_calculator->Calculate();

        TS_ASSERT_DELTA(p_vessel1->GetFlowProperties()->GetHaematocrit(p_vessel1->GetSegments()),0.45, 1e-6);
        TS_ASSERT_DELTA(p_vessel2->GetFlowProperties()->GetHaematocrit(p_vessel2->GetSegments()),0.45, 1e-6);
        TS_ASSERT_DELTA(p_vessel3->GetFlowProperties()->GetHaematocrit(p_vessel3->GetSegments()),0.9, 1e-6);
    }

    void TestBifurcationOutflowNetwork() throw(Exception)
    {
        boost::shared_ptr<VesselNode<2> > p_node1 = VesselNode<2>::Create(0.0, 0.0);
        boost::shared_ptr<VesselNode<2> > p_node2 = VesselNode<2>::Create(80.0e-6, 0.0);
        boost::shared_ptr<VesselNode<2> > p_node3 = VesselNode<2>::Create(160.0e-6, 0.0);
        boost::shared_ptr<VesselNode<2> > p_node4 = VesselNode<2>::Create(200.0e-6, 0.0);
        p_node4->GetFlowProperties()->SetIsInputNode(true);

        boost::shared_ptr<VesselSegment<2> > p_segment1(VesselSegment<2>::Create(p_node1, p_node3));
        boost::shared_ptr<VesselSegment<2> > p_segment2(VesselSegment<2>::Create(p_node2, p_node3));
        boost::shared_ptr<VesselSegment<2> > p_segment3(VesselSegment<2>::Create(p_node3, p_node4));

        boost::shared_ptr<Vessel<2> > p_vessel1(Vessel<2>::Create(p_segment1));
        boost::shared_ptr<Vessel<2> > p_vessel2(Vessel<2>::Create(p_segment2));
        boost::shared_ptr<Vessel<2> > p_vessel3(Vessel<2>::Create(p_segment3));
        p_vessel1->GetFlowProperties()->SetFlowRate(-1.0*unit::metre_cubed_per_second, p_vessel1->GetSegments());
        p_vessel2->GetFlowProperties()->SetFlowRate(-1.0*unit::metre_cubed_per_second, p_vessel2->GetSegments());
        p_vessel3->GetFlowProperties()->SetFlowRate(-1.0*unit::metre_cubed_per_second, p_vessel3->GetSegments());

        boost::shared_ptr<VesselNetwork<2> > p_network = boost::shared_ptr<VesselNetwork<2> >(new VesselNetwork<2>);
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);
        p_network->AddVessel(p_vessel3);

        boost::shared_ptr<AlarconHaematocritSolver<2> > p_haematocrit_calculator(new AlarconHaematocritSolver<2>());
        p_haematocrit_calculator->SetVesselNetwork(p_network);
        p_haematocrit_calculator->Calculate();

        TS_ASSERT_DELTA(p_vessel1->GetFlowProperties()->GetHaematocrit(p_vessel1->GetSegments()),0.15, 1e-6);
        TS_ASSERT_DELTA(p_vessel2->GetFlowProperties()->GetHaematocrit(p_vessel2->GetSegments()),0.3, 1e-6);
        TS_ASSERT_DELTA(p_vessel3->GetFlowProperties()->GetHaematocrit(p_vessel3->GetSegments()),0.45, 1e-6);
    }

    void TestBifurcationOutflowNetworkBiasedFlow() throw(Exception)
    {
        boost::shared_ptr<VesselNode<2> > p_node1 = VesselNode<2>::Create(0.0, 0.0);
        boost::shared_ptr<VesselNode<2> > p_node2 = VesselNode<2>::Create(80.0e-6, 0.0);
        boost::shared_ptr<VesselNode<2> > p_node3 = VesselNode<2>::Create(160.0e-6, 0.0);
        boost::shared_ptr<VesselNode<2> > p_node4 = VesselNode<2>::Create(200.0e-6, 0.0);
        p_node4->GetFlowProperties()->SetIsInputNode(true);

        boost::shared_ptr<VesselSegment<2> > p_segment1(VesselSegment<2>::Create(p_node1, p_node3));
        boost::shared_ptr<VesselSegment<2> > p_segment2(VesselSegment<2>::Create(p_node2, p_node3));
        boost::shared_ptr<VesselSegment<2> > p_segment3(VesselSegment<2>::Create(p_node3, p_node4));

        boost::shared_ptr<Vessel<2> > p_vessel1(Vessel<2>::Create(p_segment1));
        boost::shared_ptr<Vessel<2> > p_vessel2(Vessel<2>::Create(p_segment2));
        boost::shared_ptr<Vessel<2> > p_vessel3(Vessel<2>::Create(p_segment3));
        p_vessel1->GetFlowProperties()->SetFlowRate(-1.0*unit::metre_cubed_per_second, p_vessel1->GetSegments());
        p_vessel2->GetFlowProperties()->SetFlowRate(-3.0*unit::metre_cubed_per_second, p_vessel2->GetSegments());
        p_vessel3->GetFlowProperties()->SetFlowRate(-1.0*unit::metre_cubed_per_second, p_vessel3->GetSegments());

        boost::shared_ptr<VesselNetwork<2> > p_network = boost::shared_ptr<VesselNetwork<2> >(new VesselNetwork<2>);
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);
        p_network->AddVessel(p_vessel3);

        boost::shared_ptr<AlarconHaematocritSolver<2> > p_haematocrit_calculator(new AlarconHaematocritSolver<2>());
        p_haematocrit_calculator->SetVesselNetwork(p_network);
        p_haematocrit_calculator->Calculate();

        TS_ASSERT_DELTA(p_vessel1->GetFlowProperties()->GetHaematocrit(p_vessel1->GetSegments()),0.0, 1e-6);
        TS_ASSERT_DELTA(p_vessel2->GetFlowProperties()->GetHaematocrit(p_vessel2->GetSegments()),0.45, 1e-6);
        TS_ASSERT_DELTA(p_vessel3->GetFlowProperties()->GetHaematocrit(p_vessel3->GetSegments()),0.45, 1e-6);
    }

    void TestHexagonalNetwork() throw(Exception)
    {
        // Specify the network dimensions
        units::quantity<unit::length> vessel_length = 80.0 * 1.e-6 * unit::metres;

        // Generate the network
        VesselNetworkGenerator<2> vascular_network_generator;
        boost::shared_ptr<VesselNetwork<2> > vascular_network = vascular_network_generator.GenerateHexagonalNetwork(800.0* 1.e-6 * unit::metres,
                                                                                                                        1000.0* 1.e-6 * unit::metres,
                                                                                                                        vessel_length);
        std::vector<boost::shared_ptr<VesselNode<2> > > nodes;
        nodes.push_back(boost::shared_ptr<VesselNode<2> > (VesselNode<2>::Create(0.0, 0.0)));
        nodes.push_back(boost::shared_ptr<VesselNode<2> > (VesselNode<2>::Create(5.0, 0.0)));
        boost::shared_ptr<VesselSegment<2> > p_segment(VesselSegment<2>::Create(nodes[0], nodes[1]));

        double radius = 10.0;
        p_segment->SetRadius(radius*1.e-6*unit::metres);
        double haematocrit = 0.45;
        p_segment->GetFlowProperties()->SetHaematocrit(haematocrit);
        vascular_network->SetSegmentProperties(p_segment);

        std::pair<DimensionalChastePoint<2>, DimensionalChastePoint<2> > network_extents = vascular_network->GetExtents();
        double y_middle = (network_extents.first[1] + network_extents.second[1]) / 2.0;
        double x_middle = (network_extents.first[0] + network_extents.second[0]) / 2.0;

        std::vector<boost::shared_ptr<Vessel<2> > >::iterator vessel_iterator;

        std::vector<boost::shared_ptr<Vessel<2> > > vessels = vascular_network->GetVessels();

        for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
        {
            if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetStartNode()->rGetLocation()[1] >  y_middle)
                {
                    if((*vessel_iterator)->GetStartNode()->rGetLocation()[0] >  x_middle)
                    {
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3320*unit::pascals);
                    }
                }
            }
            if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetEndNode()->rGetLocation()[1] >  y_middle)
                {
                    if((*vessel_iterator)->GetStartNode()->rGetLocation()[0] >  x_middle)
                    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3320*unit::pascals);
                    }
                }
            }
            if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetStartNode()->rGetLocation()[1] <=  y_middle)
                {
                    if((*vessel_iterator)->GetStartNode()->rGetLocation()[0] <  x_middle)
                    {
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2090*unit::pascals);
                    }
                }
            }
            if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetEndNode()->rGetLocation()[1] <=  y_middle)
                {
                    if((*vessel_iterator)->GetStartNode()->rGetLocation()[0] <  x_middle)
                    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2090*unit::pascals);
                    }
                }
            }
        }

        std::vector<boost::shared_ptr<VesselSegment<2> > > segments = vascular_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(1.e-3*unit::poiseuille);
        }

        VesselImpedanceCalculator<2> impedance_calculator;
        impedance_calculator.SetVesselNetwork(vascular_network);
        impedance_calculator.Calculate();
        FlowSolver<2> solver;
        solver.SetVesselNetwork(vascular_network);
        solver.SetUp();
        solver.Solve();

        OutputFileHandler output_file_handler("TestHaematocritSolver", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexNet.vtp");
        vascular_network->Write(output_filename);

        boost::shared_ptr<AlarconHaematocritSolver<2> > p_haematocrit_calculator(new AlarconHaematocritSolver<2>());
        p_haematocrit_calculator->SetVesselNetwork(vascular_network);
        p_haematocrit_calculator->Calculate();

        std::string output_filename2 = output_file_handler.GetOutputDirectoryFullPath().append("HexNetHemo.vtp");
        vascular_network->Write(output_filename2);
    }
};

#endif // TESTALARCONHAEMATOCRITSOLVER_HPP
