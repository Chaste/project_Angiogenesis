//
//  TestSimpleStructuralAdaptationSolver.hpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 03/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#ifndef TestSimpleStructuralAdaptationSolver_hpp
#define TestSimpleStructuralAdaptationSolver_hpp

#include <cxxtest/TestSuite.h>
#include "StructuralAdaptationSolver.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VasculatureGenerator.hpp"
#include "FlowSolver.hpp"
#include "VasculatureData.hpp"
#include "SimulationTime.hpp"
#include "FakePetscSetup.hpp"
#include "Alarcon03HaematocritSolver.hpp"
#include "PoiseuilleImpedanceCalculator.hpp"

class TestSimpleStructuralAdaptationSolver : public CxxTest::TestSuite
{

public:

    void TestTwoVesselNetwork() throw(Exception)
    {
        boost::shared_ptr<VascularNode<2> > p_node1 = VascularNode<2>::Create(0.0, 0.0);
        boost::shared_ptr<VascularNode<2> > p_node2 = VascularNode<2>::Create(80, 0.0);
        boost::shared_ptr<VascularNode<2> > p_node3 = VascularNode<2>::Create(160, 0.0);
        p_node1->GetFlowProperties()->SetIsInputNode(true);

        boost::shared_ptr<CaVesselSegment<2> > p_segment1(CaVesselSegment<2>::Create(p_node1, p_node2));
        boost::shared_ptr<CaVesselSegment<2> > p_segment2(CaVesselSegment<2>::Create(p_node2, p_node3));

        boost::shared_ptr<CaVessel<2> > p_vessel1(CaVessel<2>::Create(p_segment1));
        boost::shared_ptr<CaVessel<2> > p_vessel2(CaVessel<2>::Create(p_segment2));
        p_vessel1->SetFlowRate(1.0);
        p_vessel2->SetFlowRate(2.0);

        boost::shared_ptr<CaVascularNetwork<2> > p_network = boost::shared_ptr<CaVascularNetwork<2> >(new CaVascularNetwork<2>);
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);

        boost::shared_ptr<Alarcon03HaematocritSolver<2> > p_haematocrit_calculator(new Alarcon03HaematocritSolver<2>());
        p_haematocrit_calculator->Calculate(p_network);

        TS_ASSERT_DELTA(p_vessel1->GetHaematocrit(),0.45, 1e-6);
        TS_ASSERT_DELTA(p_vessel2->GetHaematocrit(),0.45, 1e-6);
    }

    void TestBifurcationInflowNetwork() throw(Exception)
    {
        boost::shared_ptr<VascularNode<2> > p_node1 = VascularNode<2>::Create(0.0, 0.0);
        boost::shared_ptr<VascularNode<2> > p_node2 = VascularNode<2>::Create(80, 0.0);
        boost::shared_ptr<VascularNode<2> > p_node3 = VascularNode<2>::Create(160, 0.0);
        boost::shared_ptr<VascularNode<2> > p_node4 = VascularNode<2>::Create(200, 0.0);
        p_node1->GetFlowProperties()->SetIsInputNode(true);
        p_node2->GetFlowProperties()->SetIsInputNode(true);

        boost::shared_ptr<CaVessel<2> > p_vessel1(CaVessel<2>::Create(p_node1, p_node3));
        boost::shared_ptr<CaVessel<2> > p_vessel2(CaVessel<2>::Create(p_node2, p_node3));
        boost::shared_ptr<CaVessel<2> > p_vessel3(CaVessel<2>::Create(p_node3, p_node4));
        p_vessel1->SetFlowRate(1.0);
        p_vessel2->SetFlowRate(1.0);
        p_vessel3->SetFlowRate(1.0);

        boost::shared_ptr<CaVascularNetwork<2> > p_network = boost::shared_ptr<CaVascularNetwork<2> >(new CaVascularNetwork<2>);
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);
        p_network->AddVessel(p_vessel3);

        boost::shared_ptr<Alarcon03HaematocritSolver<2> > p_haematocrit_calculator(new Alarcon03HaematocritSolver<2>());
        p_haematocrit_calculator->Calculate(p_network);

        TS_ASSERT_DELTA(p_vessel1->GetHaematocrit(),0.45, 1e-6);
        TS_ASSERT_DELTA(p_vessel2->GetHaematocrit(),0.45, 1e-6);
        TS_ASSERT_DELTA(p_vessel3->GetHaematocrit(),0.9, 1e-6);
    }

    void TestBifurcationOutflowNetwork() throw(Exception)
    {
        boost::shared_ptr<VascularNode<2> > p_node1 = VascularNode<2>::Create(0.0, 0.0);
        boost::shared_ptr<VascularNode<2> > p_node2 = VascularNode<2>::Create(80.0e-6, 0.0);
        boost::shared_ptr<VascularNode<2> > p_node3 = VascularNode<2>::Create(160.0e-6, 0.0);
        boost::shared_ptr<VascularNode<2> > p_node4 = VascularNode<2>::Create(200.0e-6, 0.0);
        p_node4->GetFlowProperties()->SetIsInputNode(true);

        boost::shared_ptr<CaVesselSegment<2> > p_segment1(CaVesselSegment<2>::Create(p_node1, p_node3));
        boost::shared_ptr<CaVesselSegment<2> > p_segment2(CaVesselSegment<2>::Create(p_node2, p_node3));
        boost::shared_ptr<CaVesselSegment<2> > p_segment3(CaVesselSegment<2>::Create(p_node3, p_node4));

        boost::shared_ptr<CaVessel<2> > p_vessel1(CaVessel<2>::Create(p_segment1));
        boost::shared_ptr<CaVessel<2> > p_vessel2(CaVessel<2>::Create(p_segment2));
        boost::shared_ptr<CaVessel<2> > p_vessel3(CaVessel<2>::Create(p_segment3));
        p_vessel1->SetFlowRate(-1.0);
        p_vessel2->SetFlowRate(-1.0);
        p_vessel3->SetFlowRate(-1.0);

        boost::shared_ptr<CaVascularNetwork<2> > p_network = boost::shared_ptr<CaVascularNetwork<2> >(new CaVascularNetwork<2>);
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);
        p_network->AddVessel(p_vessel3);

        boost::shared_ptr<Alarcon03HaematocritSolver<2> > p_haematocrit_calculator(new Alarcon03HaematocritSolver<2>());
        p_haematocrit_calculator->Calculate(p_network);

        TS_ASSERT_DELTA(p_vessel1->GetHaematocrit(),0.15, 1e-6);
        TS_ASSERT_DELTA(p_vessel2->GetHaematocrit(),0.3, 1e-6);
        TS_ASSERT_DELTA(p_vessel3->GetHaematocrit(),0.45, 1e-6);
    }

    void TestBifurcationOutflowNetworkBiasedFlow() throw(Exception)
    {
        boost::shared_ptr<VascularNode<2> > p_node1 = VascularNode<2>::Create(0.0, 0.0);
        boost::shared_ptr<VascularNode<2> > p_node2 = VascularNode<2>::Create(80.0e-6, 0.0);
        boost::shared_ptr<VascularNode<2> > p_node3 = VascularNode<2>::Create(160.0e-6, 0.0);
        boost::shared_ptr<VascularNode<2> > p_node4 = VascularNode<2>::Create(200.0e-6, 0.0);
        p_node4->GetFlowProperties()->SetIsInputNode(true);

        boost::shared_ptr<CaVesselSegment<2> > p_segment1(CaVesselSegment<2>::Create(p_node1, p_node3));
        boost::shared_ptr<CaVesselSegment<2> > p_segment2(CaVesselSegment<2>::Create(p_node2, p_node3));
        boost::shared_ptr<CaVesselSegment<2> > p_segment3(CaVesselSegment<2>::Create(p_node3, p_node4));

        boost::shared_ptr<CaVessel<2> > p_vessel1(CaVessel<2>::Create(p_segment1));
        boost::shared_ptr<CaVessel<2> > p_vessel2(CaVessel<2>::Create(p_segment2));
        boost::shared_ptr<CaVessel<2> > p_vessel3(CaVessel<2>::Create(p_segment3));
        p_vessel1->SetFlowRate(-1.0);
        p_vessel2->SetFlowRate(-3.0);
        p_vessel3->SetFlowRate(-1.0);

        boost::shared_ptr<CaVascularNetwork<2> > p_network = boost::shared_ptr<CaVascularNetwork<2> >(new CaVascularNetwork<2>);
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);
        p_network->AddVessel(p_vessel3);

        boost::shared_ptr<Alarcon03HaematocritSolver<2> > p_haematocrit_calculator(new Alarcon03HaematocritSolver<2>());
        p_haematocrit_calculator->Calculate(p_network);

        TS_ASSERT_DELTA(p_vessel1->GetHaematocrit(),0.0, 1e-6);
        TS_ASSERT_DELTA(p_vessel2->GetHaematocrit(),0.45, 1e-6);
        TS_ASSERT_DELTA(p_vessel3->GetHaematocrit(),0.45, 1e-6);
    }
    void TestHexagonalNetworkAlarcon03Haematocrit() throw(Exception)
    {
        // Specify the network dimensions
        double vessel_length = 80.0;

        // Generate the network
        VasculatureGenerator<2> vascular_network_generator;
        boost::shared_ptr<CaVascularNetwork<2> > vascular_network = vascular_network_generator.GenerateHexagonalNetwork(800.0,
                                                                                                                        1000.0,
                                                                                                                        vessel_length);

        std::vector<ChastePoint<2> > points;
        points.push_back(ChastePoint<2>(0, 0));
        points.push_back(ChastePoint<2>(5, 0));
        std::vector<boost::shared_ptr<VascularNode<2> > > nodes;
        for(unsigned i=0; i < points.size(); i++)
        {
            nodes.push_back(boost::shared_ptr<VascularNode<2> > (VascularNode<2>::Create(points[i])));
        }

        boost::shared_ptr<CaVesselSegment<2> > p_segment(CaVesselSegment<2>::Create(nodes[0], nodes[1]));

        double radius = 10.0;
        p_segment->SetRadius(radius);
        double haematocrit = 0.45;
        p_segment->GetFlowProperties()->SetHaematocrit(haematocrit);
        vascular_network->SetSegmentProperties(p_segment);

        std::vector<std::pair<double, double> > extents = vascular_network->GetExtents();
        double y_middle = (extents[1].first + extents[1].second) /2.0;
        double x_middle = (extents[0].first + extents[0].second) /2.0;

        typename std::vector<boost::shared_ptr<CaVessel<2> > >::iterator vessel_iterator;

        std::vector<boost::shared_ptr<CaVessel<2> > > vessels = vascular_network->GetVessels();

        for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
        {
            if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetStartNode()->GetLocation()[1] >  y_middle)
                {
                    if((*vessel_iterator)->GetStartNode()->GetLocation()[0] >  x_middle)
                    {
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3320);
                    }
                }
            }
            if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetEndNode()->GetLocation()[1] >  y_middle)
                {
                    if((*vessel_iterator)->GetStartNode()->GetLocation()[0] >  x_middle)
                    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3320);
                    }
                }
            }
            if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetStartNode()->GetLocation()[1] <=  y_middle)
                {
                    if((*vessel_iterator)->GetStartNode()->GetLocation()[0] <  x_middle)
                    {
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2090);
                    }
                }
            }
            if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetEndNode()->GetLocation()[1] <=  y_middle)
                {
                    if((*vessel_iterator)->GetStartNode()->GetLocation()[0] <  x_middle)
                    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2090);
                    }
                }
            }
        }

        std::vector<boost::shared_ptr<CaVesselSegment<2> > > segments = vascular_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(1.e-3);
        }


        PoiseuilleImpedanceCalculator<2> impedance_calculator;
        impedance_calculator.Calculate(vascular_network);
        FlowSolver<2> solver;
        solver.SetUp();
        solver.Solve();

        OutputFileHandler output_file_handler("TestHaematocritSolver", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexNet.vtp");
        vascular_network->Write(output_filename);

        boost::shared_ptr<Alarcon03HaematocritSolver<2> > p_haematocrit_calculator(new Alarcon03HaematocritSolver<2>());
        p_haematocrit_calculator->Calculate(vascular_network);

        std::string output_filename2 = output_file_handler.GetOutputDirectoryFullPath().append("HexNetHemo.vtp");
        vascular_network->Write(output_filename2);
    }
};

#endif