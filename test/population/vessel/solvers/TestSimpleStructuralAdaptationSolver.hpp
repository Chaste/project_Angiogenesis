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
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VasculatureGenerator.hpp"
#include "SimpleFlowSolver.hpp"
#include "VasculatureData.hpp"
#include "SimpleStructuralAdaptationSolver.hpp"
#include "SimulationTime.hpp"
#include "FakePetscSetup.hpp"
#include "Alarcon03HaematocritSolver.hpp"

#include "Debug.hpp"

class TestSimpleStructuralAdaptationSolver : public CxxTest::TestSuite
{

    typedef boost::shared_ptr<VascularNode<2> > NodePtr2;
    typedef boost::shared_ptr<VascularNode<3> > NodePtr3;
    typedef boost::shared_ptr<CaVesselSegment<2> > SegmentPtr2;
    typedef boost::shared_ptr<CaVesselSegment<3> > SegmentPtr3;
    typedef boost::shared_ptr<CaVessel<2> > VesselPtr2;
    typedef boost::shared_ptr<CaVessel<3> > VesselPtr3;

public:

    void TestStructuralAdaptationOfTwoVesselNetwork() throw(Exception)
    {
        NodePtr2 p_node1 = VascularNode<2>::Create(0.0, 0.0);
        NodePtr2 p_node2 = VascularNode<2>::Create(80.0e-6, 0.0);
        NodePtr2 p_node3 = VascularNode<2>::Create(160.0e-6, 0.0);

        SegmentPtr2 p_segment1(CaVesselSegment<2>::Create(p_node1, p_node2));
        SegmentPtr2 p_segment2(CaVesselSegment<2>::Create(p_node2, p_node3));

        p_node1->GetFlowProperties()->SetIsInputNode(true);
        p_node1->GetFlowProperties()->SetPressure(3322);
        p_node3->GetFlowProperties()->SetIsOutputNode(true);
        p_node3->GetFlowProperties()->SetPressure(1993);

        VesselPtr2 p_vessel1(CaVessel<2>::Create(p_segment1));
        VesselPtr2 p_vessel2(CaVessel<2>::Create(p_segment2));

        boost::shared_ptr<CaVascularNetwork<2> > p_network = boost::shared_ptr<CaVascularNetwork<2> >(new CaVascularNetwork<2>);
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);

        double radius = 10.0e-6;
        p_segment1->SetRadius(radius);
        double haematocrit = 0.45;
        p_segment1->GetFlowProperties()->SetHaematocrit(haematocrit);
        p_network->SetSegmentProperties(p_segment1);

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(30, 1);

        // Write the network to file
        OutputFileHandler output_file_handler("TestSimpleStructuralAdaptationSolver", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("TwoVesselNetwork.vtp");
        std::string progress_output_filename = output_file_handler.GetOutputDirectoryFullPath().append("TwoVesselNetwork_SAAProgress.dat");

        SimpleStructuralAdaptationSolver<2> solver;
        solver.SetWriteOutput(true);
        solver.SetOutputFileName(progress_output_filename);
        solver.SetTolerance(0.0001);
        solver.SetTimeIncrement(0.0001);
        //		solver.SetMaxIterations(1);
        solver.Implement(p_network);

        // Write the network to file
        p_network->Write(output_filename);

        TS_ASSERT_DELTA(p_node2->GetFlowProperties()->GetPressure(),(3322.0 + 1993.0) / 2.0, 1e-6);
        TS_ASSERT_DELTA(abs(p_segment1->GetFlowProperties()->GetFlowRate()),abs(p_segment2->GetFlowProperties()->GetFlowRate()),1e-6);
        p_simulation_time->Destroy();
    }

    void TestStructuralAdaptationOfOneVesselNetwork() throw(Exception)
    {
        NodePtr2 p_node1 = VascularNode<2>::Create(0.0, 0.0);
        NodePtr2 p_node2 = VascularNode<2>::Create(80.0e-6, 0.0);

        SegmentPtr2 p_segment1(CaVesselSegment<2>::Create(p_node1, p_node2));

        p_node1->GetFlowProperties()->SetIsInputNode(true);
        p_node1->GetFlowProperties()->SetPressure(3322);
        p_node2->GetFlowProperties()->SetIsOutputNode(true);
        p_node2->GetFlowProperties()->SetPressure(1993);

        VesselPtr2 p_vessel1(CaVessel<2>::Create(p_segment1));

        boost::shared_ptr<CaVascularNetwork<2> > p_network = boost::shared_ptr<CaVascularNetwork<2> >(new CaVascularNetwork<2>);
        p_network->AddVessel(p_vessel1);

        double radius = 10.0e-6;
        p_segment1->SetRadius(radius);
        double haematocrit = 0.45;
        p_segment1->GetFlowProperties()->SetHaematocrit(haematocrit);
        p_network->SetSegmentProperties(p_segment1);

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(30, 1);

        // Write the network to file
        OutputFileHandler output_file_handler("TestSimpleStructuralAdaptationSolver", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("OneVesselNetwork.vtp");
        std::string progress_output_filename = output_file_handler.GetOutputDirectoryFullPath().append("OneVesselNetwork_SAAProgress.dat");

        SimpleStructuralAdaptationSolver<2> solver;
        solver.SetWriteOutput(true);
        solver.SetOutputFileName(progress_output_filename);
        solver.SetTolerance(0.0001);
        solver.SetTimeIncrement(0.0001);
        //        solver.SetMaxIterations(1);
        solver.Implement(p_network);

        // Write the network to file
        p_network->Write(output_filename);

        p_simulation_time->Destroy();
    }

    void TestStructuralAdaptationOfHexagonalNetworkConstantHaematocrit() throw(Exception)
	{
        // Specify the network dimensions
        double vessel_length = 82.0e-6;

        // Generate the network
        VasculatureGenerator<2> vascular_network_generator;
        boost::shared_ptr<CaVascularNetwork<2> > vascular_network = vascular_network_generator.GenerateHexagonalNetwork(800.0e-6,
                                                                                                                        1000.0e-6, vessel_length);

        std::vector<ChastePoint<2> > points;
        points.push_back(ChastePoint<2>(0, 0));
        points.push_back(ChastePoint<2>(5, 0));

        std::vector<NodePtr2> nodes;
        for(unsigned i=0; i < points.size(); i++)
        {
            nodes.push_back(NodePtr2 (VascularNode<2>::Create(points[i])));
        }

        SegmentPtr2 p_segment(CaVesselSegment<2>::Create(nodes[0], nodes[1]));

        double radius = 10.0e-6;
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

        /*
         * Set up SimulationTime object.  This must be set up in order for structural adaptation algorithm
         * to run.
         */
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(30, 1);

        // Write the network to file
        OutputFileHandler output_file_handler("TestSimpleStructuralAdaptationSolver", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexagonalVesselNetwork_ConstantHaematocrit.vtp");
        std::string progress_output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexagonalVesselNetwork_ConstantHaematocrit_SAAProgress.dat");

        SimpleStructuralAdaptationSolver<2> solver;
        solver.SetWriteOutput(true);
        solver.SetOutputFileName(progress_output_filename);
        solver.SetTolerance(0.001);
        solver.SetTimeIncrement(0.0001);
        solver.SetMaxIterations(1000);
        solver.Implement(vascular_network);

        // Write the network to file
        vascular_network->Write(output_filename);

        p_simulation_time->Destroy();
	}
};

#endif
