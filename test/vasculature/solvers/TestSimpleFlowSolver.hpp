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

#ifndef TESTSIMPLEFLOWSOLVER_HPP_
#define TESTSIMPLEFLOWSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VasculatureGenerator.hpp"
#include "FakePetscSetup.hpp"
#include "SimpleFlowSolver.hpp"
#include "VasculatureData.hpp"

#include "Debug.hpp"

class TestSimpleFlowSolver : public CxxTest::TestSuite
{

	typedef boost::shared_ptr<VascularNode<2> > NodePtr2;
	typedef boost::shared_ptr<VascularNode<3> > NodePtr3;
	typedef boost::shared_ptr<CaVesselSegment<2> > SegmentPtr2;
	typedef boost::shared_ptr<CaVesselSegment<3> > SegmentPtr3;
	typedef boost::shared_ptr<CaVessel<2> > VesselPtr2;
	typedef boost::shared_ptr<CaVessel<3> > VesselPtr3;

public:



	void TestFlowThroughSingleSegment() throw(Exception)
	{

		// Make some nodes
		std::vector<ChastePoint<3> > points;
		points.push_back(ChastePoint<3>(0, 0, 0));
		points.push_back(ChastePoint<3>(5, 0, 0));

		std::vector<NodePtr3> nodes;
		for(unsigned i=0; i < points.size(); i++)
		{
			nodes.push_back(NodePtr3 (VascularNode<3>::Create(points[i])));
		}

		SegmentPtr3 p_segment(CaVesselSegment<3>::Create(nodes[0], nodes[1]));
		VesselPtr3 p_vessel(CaVessel<3>::Create(p_segment));

		// Generate the network
		boost::shared_ptr<CaVascularNetwork<3> > p_vascular_network(new CaVascularNetwork<3>());

		p_vascular_network->AddVessel(p_vessel);

		VasculatureData data;
		double impedance = 10.0;
		data.SetData("Impedance", impedance);
		p_vascular_network->SetSegmentData(data);

		VasculatureData node_data;
		bool is_input = false;
		node_data.SetData("Is Input", is_input);
		bool is_output = false;
		node_data.SetData("Is Output", is_output);
		p_vascular_network->SetNodeData(node_data);

		p_vessel->GetStartNode()->SetData<bool>("Is Input", true);
		p_vessel->GetStartNode()->SetData<double>("Pressure", 3393);

		p_vessel->GetEndNode()->SetData<bool>("Is Output", true);
		p_vessel->GetEndNode()->SetData<double>("Pressure", 1000.5);

		SimpleFlowSolver<3> solver;

		solver.Implement(p_vascular_network);

		TS_ASSERT_DELTA(p_vessel->GetStartNode()->GetData<double>("Pressure"),3393,1e-6);
		TS_ASSERT_DELTA(p_vessel->GetEndNode()->GetData<double>("Pressure"),1000.5,1e-6);

		TS_ASSERT_DELTA(p_vessel->GetData<double>("Flow Rate"),(3393-1000.5)/impedance,1e-6);
		TS_ASSERT_DELTA(p_segment->GetData<double>("Flow Rate"),(3393-1000.5)/impedance,1e-6);

		p_segment->SetData("Impedance",0.0);

		TS_ASSERT_THROWS_THIS(solver.Implement(p_vascular_network),"Impedance should be a positive number.");

	}

	void TestFlowThroughSingleVesselWithMultipleSegments() throw(Exception)
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
		boost::shared_ptr<CaVascularNetwork<3> > p_vascular_network(new CaVascularNetwork<3>());

		p_vascular_network->AddVessel(p_vessel);

		VasculatureData data;
		double impedance = 10.0;
		data.SetData("Impedance", impedance);
		p_vascular_network->SetSegmentData(data);

		VasculatureData node_data;
		bool is_input = false;
		node_data.SetData("Is Input", is_input);
		bool is_output = false;
		node_data.SetData("Is Output", is_output);
		p_vascular_network->SetNodeData(node_data);

		p_vessel->GetStartNode()->SetData<bool>("Is Input", true);
		p_vessel->GetStartNode()->SetData<double>("Pressure", 3393);

		p_vessel->GetEndNode()->SetData<bool>("Is Output", true);
		p_vessel->GetEndNode()->SetData<double>("Pressure", 1000.5);

		SimpleFlowSolver<3> solver;

		solver.Implement(p_vascular_network);

		for (unsigned i = 0; i < nodes.size(); i++)
		{
			TS_ASSERT_DELTA(nodes[i]->GetData<double>("Pressure"),3393 - (3393-1000.5)*i/(nodes.size()-1),1e-6);
		}

		TS_ASSERT_DELTA(p_vessel->GetStartNode()->GetData<double>("Pressure"),3393,1e-6);
		TS_ASSERT_DELTA(p_vessel->GetEndNode()->GetData<double>("Pressure"),1000.5,1e-6);

		TS_ASSERT_DELTA(p_vessel->GetData<double>("Flow Rate"),(3393-1000.5)/(segments.size()*impedance),1e-6);

		for (unsigned i = 0; i < segments.size(); i++)
		{
			TS_ASSERT_DELTA(segments[i]->GetData<double>("Flow Rate"),(3393-1000.5)/(segments.size()*impedance),1e-6);
		}


	}

	void TestFlowThroughBifurcation() throw(Exception)
	{

		// Make some nodes
		std::vector<ChastePoint<3> > points;
		points.push_back(ChastePoint<3>(0, 0, 0)); // input
		points.push_back(ChastePoint<3>(0, 1.0, 0)); // input
		points.push_back(ChastePoint<3>(1.0, 0.5, 0)); // bifurcation
		points.push_back(ChastePoint<3>(1.0, 1.0, 0)); // output

		std::vector<NodePtr3> nodes;
		for(unsigned i=0; i < points.size(); i++)
		{
			nodes.push_back(NodePtr3 (VascularNode<3>::Create(points[i])));
		}

		SegmentPtr3 p_segment1(CaVesselSegment<3>::Create(nodes[0], nodes[2]));
		SegmentPtr3 p_segment2(CaVesselSegment<3>::Create(nodes[1], nodes[2]));
		SegmentPtr3 p_segment3(CaVesselSegment<3>::Create(nodes[3], nodes[2]));

		VesselPtr3 p_vessel1(CaVessel<3>::Create(p_segment1));
		VesselPtr3 p_vessel2(CaVessel<3>::Create(p_segment2));
		VesselPtr3 p_vessel3(CaVessel<3>::Create(p_segment3));

		std::vector<VesselPtr3> vessels;
		vessels.push_back(p_vessel1); // lower input vessel
		vessels.push_back(p_vessel2); // upper input vessel
		vessels.push_back(p_vessel3); // output vessel

		// Generate the network
		boost::shared_ptr<CaVascularNetwork<3> > p_vascular_network(new CaVascularNetwork<3>());

		p_vascular_network->AddVessels(vessels);

		VasculatureData data;
		double impedance = 10.0;
		data.SetData("Impedance", impedance);
		p_vascular_network->SetSegmentData(data);

		VasculatureData node_data;
		bool is_input = false;
		node_data.SetData("Is Input", is_input);
		bool is_output = false;
		node_data.SetData("Is Output", is_output);
		p_vascular_network->SetNodeData(node_data);

		nodes[0]->SetData<bool>("Is Input", true);
		nodes[0]->SetData<double>("Pressure", 3393);

		nodes[1]->SetData<bool>("Is Input", true);
		nodes[1]->SetData<double>("Pressure", 3393);

		nodes[3]->SetData<bool>("Is Output", true);
		nodes[3]->SetData<double>("Pressure", 1000.5);

		SimpleFlowSolver<3> solver;

		solver.Implement(p_vascular_network);


		TS_ASSERT_DELTA(nodes[0]->GetData<double>("Pressure"),3393,1e-6);
		TS_ASSERT_DELTA(nodes[1]->GetData<double>("Pressure"),3393,1e-6);
		TS_ASSERT_DELTA(nodes[2]->GetData<double>("Pressure"),(2*3393 + 1000.5)/3,1e-6);
		TS_ASSERT_DELTA(nodes[3]->GetData<double>("Pressure"),1000.5,1e-6);

		TS_ASSERT_DELTA(vessels[0]->GetData<double>("Flow Rate"),(3393-nodes[2]->GetData<double>("Pressure"))/impedance,1e-6);
		TS_ASSERT_DELTA(vessels[1]->GetData<double>("Flow Rate"),(3393-nodes[2]->GetData<double>("Pressure"))/impedance,1e-6);
		TS_ASSERT_DELTA(vessels[2]->GetData<double>("Flow Rate"),-(nodes[2]->GetData<double>("Pressure")-1000.5)/impedance,1e-6);

		TS_ASSERT_DELTA(p_segment1->GetData<double>("Flow Rate"),(3393-nodes[2]->GetData<double>("Pressure"))/impedance,1e-6);
		TS_ASSERT_DELTA(p_segment2->GetData<double>("Flow Rate"),(3393-nodes[2]->GetData<double>("Pressure"))/impedance,1e-6);
		TS_ASSERT_DELTA(p_segment3->GetData<double>("Flow Rate"),-(nodes[2]->GetData<double>("Pressure")-1000.5)/impedance,1e-6);

		TS_ASSERT_DELTA(vessels[0]->GetData<double>("Absolute Flow Rate"),(3393-nodes[2]->GetData<double>("Pressure"))/impedance,1e-6);
		TS_ASSERT_DELTA(vessels[1]->GetData<double>("Absolute Flow Rate"),(3393-nodes[2]->GetData<double>("Pressure"))/impedance,1e-6);
		TS_ASSERT_DELTA(vessels[2]->GetData<double>("Absolute Flow Rate"),fabs((1000.5-nodes[2]->GetData<double>("Pressure")))/impedance,1e-6);

		TS_ASSERT_DELTA(p_segment1->GetData<double>("Absolute Flow Rate"),(3393-nodes[2]->GetData<double>("Pressure"))/impedance,1e-6);
		TS_ASSERT_DELTA(p_segment2->GetData<double>("Absolute Flow Rate"),(3393-nodes[2]->GetData<double>("Pressure"))/impedance,1e-6);
		TS_ASSERT_DELTA(p_segment3->GetData<double>("Absolute Flow Rate"),fabs((1000.5-nodes[2]->GetData<double>("Pressure")))/impedance,1e-6);

		double kirchoff_residual = vessels[0]->GetData<double>("Flow Rate") + vessels[1]->GetData<double>("Flow Rate") +
									vessels[2]->GetData<double>("Flow Rate");

		TS_ASSERT_DELTA(kirchoff_residual,0,1e-6);

	}

	void TestFlowThroughBifurcationHavingSwappedNodeLabels() throw(Exception)
	{

		// Make some nodes
		std::vector<ChastePoint<3> > points;
		points.push_back(ChastePoint<3>(0, 0, 0)); // input
		points.push_back(ChastePoint<3>(0, 1.0, 0)); // input
		points.push_back(ChastePoint<3>(1.0, 0.5, 0)); // bifurcation
		points.push_back(ChastePoint<3>(1.0, 1.0, 0)); // output

		std::vector<NodePtr3> nodes;
		for(unsigned i=0; i < points.size(); i++)
		{
			nodes.push_back(NodePtr3 (VascularNode<3>::Create(points[i])));
		}

		SegmentPtr3 p_segment1(CaVesselSegment<3>::Create(nodes[0], nodes[2]));
		SegmentPtr3 p_segment2(CaVesselSegment<3>::Create(nodes[1], nodes[2]));
		SegmentPtr3 p_segment3(CaVesselSegment<3>::Create(nodes[2], nodes[3]));

		VesselPtr3 p_vessel1(CaVessel<3>::Create(p_segment1));
		VesselPtr3 p_vessel2(CaVessel<3>::Create(p_segment2));
		VesselPtr3 p_vessel3(CaVessel<3>::Create(p_segment3));

		std::vector<VesselPtr3> vessels;
		vessels.push_back(p_vessel1); // lower input vessel
		vessels.push_back(p_vessel2); // upper input vessel
		vessels.push_back(p_vessel3); // output vessel

		// Generate the network
		boost::shared_ptr<CaVascularNetwork<3> > p_vascular_network(new CaVascularNetwork<3>());

		p_vascular_network->AddVessels(vessels);

		VasculatureData data;
		double impedance = 10.0;
		data.SetData("Impedance", impedance);
		p_vascular_network->SetSegmentData(data);

		VasculatureData node_data;
		bool is_input = false;
		node_data.SetData("Is Input", is_input);
		bool is_output = false;
		node_data.SetData("Is Output", is_output);
		p_vascular_network->SetNodeData(node_data);

		nodes[0]->SetData<bool>("Is Input", true);
		nodes[0]->SetData<double>("Pressure", 3393);

		nodes[1]->SetData<bool>("Is Input", true);
		nodes[1]->SetData<double>("Pressure", 3393);

		nodes[3]->SetData<bool>("Is Output", true);
		nodes[3]->SetData<double>("Pressure", 1000.5);

		SimpleFlowSolver<3> solver;

		solver.Implement(p_vascular_network);


		TS_ASSERT_DELTA(nodes[0]->GetData<double>("Pressure"),3393,1e-6);
		TS_ASSERT_DELTA(nodes[1]->GetData<double>("Pressure"),3393,1e-6);
		TS_ASSERT_DELTA(nodes[2]->GetData<double>("Pressure"),(2*3393 + 1000.5)/3,1e-6);
		TS_ASSERT_DELTA(nodes[3]->GetData<double>("Pressure"),1000.5,1e-6);

		TS_ASSERT_DELTA(vessels[0]->GetData<double>("Flow Rate"),(3393-nodes[2]->GetData<double>("Pressure"))/impedance,1e-6);
		TS_ASSERT_DELTA(vessels[1]->GetData<double>("Flow Rate"),(3393-nodes[2]->GetData<double>("Pressure"))/impedance,1e-6);
		TS_ASSERT_DELTA(vessels[2]->GetData<double>("Flow Rate"),(nodes[2]->GetData<double>("Pressure")-1000.5)/impedance,1e-6);

		TS_ASSERT_DELTA(p_segment1->GetData<double>("Flow Rate"),(3393-nodes[2]->GetData<double>("Pressure"))/impedance,1e-6);
		TS_ASSERT_DELTA(p_segment2->GetData<double>("Flow Rate"),(3393-nodes[2]->GetData<double>("Pressure"))/impedance,1e-6);
		TS_ASSERT_DELTA(p_segment3->GetData<double>("Flow Rate"),(nodes[2]->GetData<double>("Pressure")-1000.5)/impedance,1e-6);

		TS_ASSERT_DELTA(vessels[0]->GetData<double>("Absolute Flow Rate"),(3393-nodes[2]->GetData<double>("Pressure"))/impedance,1e-6);
		TS_ASSERT_DELTA(vessels[1]->GetData<double>("Absolute Flow Rate"),(3393-nodes[2]->GetData<double>("Pressure"))/impedance,1e-6);
		TS_ASSERT_DELTA(vessels[2]->GetData<double>("Absolute Flow Rate"),fabs((1000.5-nodes[2]->GetData<double>("Pressure")))/impedance,1e-6);

		TS_ASSERT_DELTA(p_segment1->GetData<double>("Absolute Flow Rate"),(3393-nodes[2]->GetData<double>("Pressure"))/impedance,1e-6);
		TS_ASSERT_DELTA(p_segment2->GetData<double>("Absolute Flow Rate"),(3393-nodes[2]->GetData<double>("Pressure"))/impedance,1e-6);
		TS_ASSERT_DELTA(p_segment3->GetData<double>("Absolute Flow Rate"),fabs((1000.5-nodes[2]->GetData<double>("Pressure")))/impedance,1e-6);

		double kirchoff_residual = vessels[0]->GetData<double>("Flow Rate") + vessels[1]->GetData<double>("Flow Rate") -
				vessels[2]->GetData<double>("Flow Rate");

		TS_ASSERT_DELTA(kirchoff_residual,0,1e-6);

	}

	void TestFlowThroughHexagonalNetwork() throw(Exception)
    {
		// Specify the network dimensions
		double vessel_length = 80.0;

		// Generate the network
		VasculatureGenerator<2> vascular_network_generator;
		boost::shared_ptr<CaVascularNetwork<2> > vascular_network = vascular_network_generator.GenerateHexagonalNetwork(1000,1000,vessel_length);

		VasculatureData data;
		double impedance = 10.0;
		data.SetData("Impedance", impedance);
		vascular_network->SetSegmentData(data);

		VasculatureData node_data;
		bool is_input = false;
		node_data.SetData("Is Input", is_input);
		bool is_output = false;
		node_data.SetData("Is Output", is_output);
		vascular_network->SetNodeData(node_data);

		std::vector<std::pair<double, double> > extents = vascular_network->GetExtents();
		double y_middle = (extents[1].first + extents[1].second) /2.0;

		typename std::vector<boost::shared_ptr<CaVessel<2> > >::iterator vessel_iterator;

		std::vector<boost::shared_ptr<CaVessel<2> > > vessels = vascular_network->GetVessels();

		for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
		{
			if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
			{
				if((*vessel_iterator)->GetStartNode()->GetLocation()[1] >  y_middle)
				{
					(*vessel_iterator)->GetStartNode()->SetData<bool>("Is Input", true);
					(*vessel_iterator)->GetStartNode()->SetData<double>("Pressure", 3393);
				}
			}
			if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
			{
				if((*vessel_iterator)->GetEndNode()->GetLocation()[1] >  y_middle)
				{
					(*vessel_iterator)->GetEndNode()->SetData<bool>("Is Input", true);
					(*vessel_iterator)->GetEndNode()->SetData<double>("Pressure", 3393);
				}
			}
			if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
			{
				if((*vessel_iterator)->GetStartNode()->GetLocation()[1] <=  y_middle)
				{
					(*vessel_iterator)->GetStartNode()->SetData<bool>("Is Output", true);
					(*vessel_iterator)->GetStartNode()->SetData<double>("Pressure", 1993);
				}
			}
			if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
			{
				if((*vessel_iterator)->GetEndNode()->GetLocation()[1] <=  y_middle)
				{
					(*vessel_iterator)->GetEndNode()->SetData<bool>("Is Output", true);
					(*vessel_iterator)->GetEndNode()->SetData<double>("Pressure", 1993);
				}
			}

		}

		SimpleFlowSolver<2> solver;
		solver.Implement(vascular_network);

		// Write the network to file
		OutputFileHandler output_file_handler("TestSimpleFlowSolver", false);
		std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexagonalVesselNetwork.vtp");
		vascular_network->Write(output_filename);
    }

};

#endif /*TESTSIMPLEFLOWSOLVER_HPP_*/
