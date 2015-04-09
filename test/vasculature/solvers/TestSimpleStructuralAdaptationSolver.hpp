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
#include "FakePetscSetup.hpp"
#include "SimpleFlowSolver.hpp"
#include "VasculatureData.hpp"
#include "SimpleStructuralAdaptationSolver.hpp"
#include "SimulationTime.hpp"

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

	void TestStructuralAdaptationOfHexagonalNetwork() throw(Exception)
	{
		// Specify the network dimensions
		double vessel_length = 80e-6;

		// Generate the network
		VasculatureGenerator<2> vascular_network_generator;
		boost::shared_ptr<CaVascularNetwork<2> > vascular_network = vascular_network_generator.GenerateHexagonalNetwork(1e-3,1e-3,vessel_length);

		VasculatureData data;
		double radius = 1e-6;
		data.SetData("Radius", radius);
		double haematocrit = 0.45;
		data.SetData("Haematocrit Level", haematocrit);
		data.SetData("Upstream Conducted Stimulus", 0.0); // Take this out when the calculator has been implemented
		data.SetData("Downstream Conducted Stimulus", 0.0); // Take this out when the calculator has been implemented
		vascular_network->SetSegmentData(data);

		VasculatureData node_data;
		bool is_input = false;
		node_data.SetData("Is Input", is_input);
		bool is_output = false;
		node_data.SetData("Is Output", is_output);
		vascular_network->SetNodeData(node_data);

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
						(*vessel_iterator)->GetStartNode()->SetData<bool>("Is Input", true);
						(*vessel_iterator)->GetStartNode()->SetData<double>("Pressure", 3322);
					}
				}
			}
			if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
			{
				if((*vessel_iterator)->GetEndNode()->GetLocation()[1] >  y_middle)
				{
					if((*vessel_iterator)->GetStartNode()->GetLocation()[0] >  x_middle)
					{
						(*vessel_iterator)->GetEndNode()->SetData<bool>("Is Input", true);
						(*vessel_iterator)->GetEndNode()->SetData<double>("Pressure", 3322);
					}
				}
			}
			if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
			{
				if((*vessel_iterator)->GetStartNode()->GetLocation()[1] <=  y_middle)
				{
					if((*vessel_iterator)->GetStartNode()->GetLocation()[0] <  x_middle)
					{
						(*vessel_iterator)->GetStartNode()->SetData<bool>("Is Output", true);
						(*vessel_iterator)->GetStartNode()->SetData<double>("Pressure", 1993);
					}
				}
			}
			if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
			{
				if((*vessel_iterator)->GetEndNode()->GetLocation()[1] <=  y_middle)
				{
					if((*vessel_iterator)->GetStartNode()->GetLocation()[0] <  x_middle)
					{
						(*vessel_iterator)->GetEndNode()->SetData<bool>("Is Output", true);
						(*vessel_iterator)->GetEndNode()->SetData<double>("Pressure", 1993);
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
		std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexagonalVesselNetwork.vtp");
		std::string progress_output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexagonalVesselNetwork_SAAProgress.dat");


		SimpleStructuralAdaptationSolver<2> solver;
		solver.SetWriteOutput(true);
		solver.SetOutputFileName(progress_output_filename);
		solver.SetTolerance(0.0001);
		solver.SetTimeIncrement(0.0001);
		solver.SetMaxIterations(1000);

		solver.Implement(vascular_network);

		// Write the network to file
		vascular_network->Write(output_filename);
	}

};

#endif
