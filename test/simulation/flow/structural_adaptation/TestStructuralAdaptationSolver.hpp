//
//  TestStructuralAdaptationSolver.hpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 03/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#ifndef TestStructuralAdaptationSolver_hpp
#define TestStructuralAdaptationSolver_hpp

#include <cxxtest/TestSuite.h>
#include "StructuralAdaptationSolver.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VasculatureGenerator.hpp"
#include "FlowSolver.hpp"
#include "VasculatureData.hpp"
#include "SimulationTime.hpp"
#include "Alarcon03HaematocritSolver.hpp"
#include "FakePetscSetup.hpp"

class TestStructuralAdaptationSolver : public CxxTest::TestSuite
{

public:

    void TestMultiVesselNetwork() throw(Exception)
    {
        std::vector<boost::shared_ptr<VascularNode<2> > > nodes;
        for(unsigned idx=0; idx<8; idx++)
        {
            nodes.push_back(VascularNode<2>::Create(double(idx*10.0), 0.0));
        }

        nodes[0]->GetFlowProperties()->SetIsInputNode(true);
        nodes[0]->GetFlowProperties()->SetPressure(3322);
        nodes[7]->GetFlowProperties()->SetIsOutputNode(true);
        nodes[7]->GetFlowProperties()->SetPressure(1993);

        double haematocrit = 0.45;
        std::vector<boost::shared_ptr<VesselSegment<2> > > segments;
        for(unsigned idx=0; idx<7; idx++)
        {
            segments.push_back(VesselSegment<2>::Create(nodes[idx], nodes[idx+1]));
            segments[idx]->GetFlowProperties()->SetHaematocrit(haematocrit);
            segments[idx]->SetRadius(10.0);
        }

        std::vector<boost::shared_ptr<Vessel<2> > > vessels;
        for(unsigned idx=0; idx<7; idx++)
        {
            vessels.push_back(Vessel<2>::Create(segments[idx]));
        }

        boost::shared_ptr<VascularNetwork<2> > p_network = VascularNetwork<2>::Create();
        p_network->AddVessels(vessels);

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(30, 1);

        // Write the network to file
        OutputFileHandler output_file_handler("TestStructuralAdaptationSolver", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("MultiVesselNetwork.vtp");
        std::string progress_output_filename = output_file_handler.GetOutputDirectoryFullPath().append("MultiVesselNetwork_SAAProgress.dat");

        StructuralAdaptationSolver<2> solver;
        solver.SetVesselNetwork(p_network);
        solver.SetWriteOutput(true);
        solver.SetOutputFileName(progress_output_filename);
        solver.SetTolerance(0.0001);
        solver.SetTimeIncrement(0.0001);
        solver.Solve();

        // Write the network to file
        p_network->Write(output_filename);

        TS_ASSERT_DELTA((nodes[3]->GetFlowProperties()->GetPressure() + nodes[4]->GetFlowProperties()->GetPressure())/2.0,(3322.0 + 1993.0) / 2.0, 1e-6);
        TS_ASSERT_DELTA(fabs(segments[0]->GetFlowProperties()->GetFlowRate()),fabs(segments[1]->GetFlowProperties()->GetFlowRate()),1e-6);
        p_simulation_time->Destroy();
    }

    void TestOneVesselNetwork() throw(Exception)
    {
        boost::shared_ptr<VascularNode<2> > p_node1 = VascularNode<2>::Create(0.0, 0.0);
        boost::shared_ptr<VascularNode<2> > p_node2 = VascularNode<2>::Create(80.0e-6, 0.0);
        boost::shared_ptr<VesselSegment<2> > p_segment1(VesselSegment<2>::Create(p_node1, p_node2));

        p_node1->GetFlowProperties()->SetIsInputNode(true);
        p_node1->GetFlowProperties()->SetPressure(3322);
        p_node2->GetFlowProperties()->SetIsOutputNode(true);
        p_node2->GetFlowProperties()->SetPressure(1993);

        boost::shared_ptr<Vessel<2> > p_vessel1(Vessel<2>::Create(p_segment1));

        boost::shared_ptr<VascularNetwork<2> > p_network = boost::shared_ptr<VascularNetwork<2> >(new VascularNetwork<2>);
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
        OutputFileHandler output_file_handler("TestStructuralAdaptationSolver", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("OneVesselNetwork.vtp");
        std::string progress_output_filename = output_file_handler.GetOutputDirectoryFullPath().append("OneVesselNetwork_SAAProgress.dat");

        StructuralAdaptationSolver<2> solver;
        solver.SetVesselNetwork(p_network);
        solver.SetWriteOutput(true);
        solver.SetOutputFileName(progress_output_filename);
        solver.SetTolerance(0.0001);
        solver.SetTimeIncrement(0.0001);
        solver.Solve();

        // Write the network to file
        p_network->Write(output_filename);
        p_simulation_time->Destroy();
    }

    void TestHexagonalNetwork() throw(Exception)
	{
        // Specify the network dimensions
        double vessel_length = 82.0;

        // Generate the network
        VasculatureGenerator<2> vascular_network_generator;
        boost::shared_ptr<VascularNetwork<2> > vascular_network = vascular_network_generator.GenerateHexagonalNetwork(800.0,
                                                                                                                        1000.0, vessel_length);

        std::vector<ChastePoint<2> > points;
        points.push_back(ChastePoint<2>(0, 0));
        points.push_back(ChastePoint<2>(5, 0));

        std::vector<boost::shared_ptr<VascularNode<2> > > nodes;
        for(unsigned i=0; i < points.size(); i++)
        {
            nodes.push_back(boost::shared_ptr<VascularNode<2> > (VascularNode<2>::Create(points[i])));
        }

        boost::shared_ptr<VesselSegment<2> > p_segment(VesselSegment<2>::Create(nodes[0], nodes[1]));

        double radius = 10.0;
        p_segment->SetRadius(radius);
        double haematocrit = 0.45;
        p_segment->GetFlowProperties()->SetHaematocrit(haematocrit);
        vascular_network->SetSegmentProperties(p_segment);

        std::vector<std::pair<double, double> > extents = vascular_network->GetExtents();
        double y_middle = (extents[1].first + extents[1].second) /2.0;
        double x_middle = (extents[0].first + extents[0].second) /2.0;

        std::vector<boost::shared_ptr<Vessel<2> > >::iterator vessel_iterator;

        std::vector<boost::shared_ptr<Vessel<2> > > vessels = vascular_network->GetVessels();

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

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(30, 1);

        // Write the network to file
        OutputFileHandler output_file_handler("TestStructuralAdaptationSolver", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexagonalVesselNetwork.vtp");
        std::string progress_output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexagonalVesselNetwork_SAAProgress.dat");

        StructuralAdaptationSolver<2> solver;
        solver.SetVesselNetwork(vascular_network);
        solver.SetWriteOutput(true);
        solver.SetOutputFileName(progress_output_filename);
        solver.SetTolerance(0.0001);
        solver.SetTimeIncrement(0.001);
        solver.SetMaxIterations(10000);
        solver.Solve();

        // Write the network to file
        vascular_network->Write(output_filename);
        p_simulation_time->Destroy();
	}
};

#endif
