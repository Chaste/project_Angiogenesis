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
#include "SimpleFlowSolver.hpp"
#include "VasculatureData.hpp"
#include "FakePetscSetup.hpp"

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

        double impedance = 1.e14;
		p_segment->GetFlowProperties()->SetImpedance(impedance);
		p_vascular_network->SetSegmentProperties(p_segment);

		p_vessel->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
		p_vessel->GetStartNode()->GetFlowProperties()->SetPressure(3393);

		p_vessel->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
		p_vessel->GetEndNode()->GetFlowProperties()->SetPressure(1000.5);

		SimpleFlowSolver<3> solver;
		solver.SetUp(p_vascular_network);
		solver.Implement();

		TS_ASSERT_DELTA(p_vessel->GetStartNode()->GetFlowProperties()->GetPressure(),3393,1e-6);
		TS_ASSERT_DELTA(p_vessel->GetEndNode()->GetFlowProperties()->GetPressure(),1000.5,1e-6);

		TS_ASSERT_DELTA(p_vessel->GetFlowRate(),(3393-1000.5)/impedance,1e-6);
		TS_ASSERT_DELTA(p_segment->GetFlowProperties()->GetFlowRate(),(3393-1000.5)/impedance,1e-6);

		p_segment->GetFlowProperties()->SetImpedance(-1.0);
		TS_ASSERT_THROWS_THIS(solver.UpdateImpedances(),"Impedance should be a positive number.");
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

        double impedance = 1.e14;
        p_segment1->GetFlowProperties()->SetImpedance(impedance);
        p_vascular_network->SetSegmentProperties(p_segment1);

		p_vessel->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
		p_vessel->GetStartNode()->GetFlowProperties()->SetPressure(3393);

		p_vessel->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
		p_vessel->GetEndNode()->GetFlowProperties()->SetPressure(1000.5);

		SimpleFlowSolver<3> solver;
		solver.SetUp(p_vascular_network);
		solver.Implement();

		for (unsigned i = 0; i < nodes.size(); i++)
		{
			TS_ASSERT_DELTA(nodes[i]->GetFlowProperties()->GetPressure(),3393 - (3393-1000.5)*i/(nodes.size()-1),1e-6);
		}

		TS_ASSERT_DELTA(p_vessel->GetStartNode()->GetFlowProperties()->GetPressure(),3393,1e-6);
		TS_ASSERT_DELTA(p_vessel->GetEndNode()->GetFlowProperties()->GetPressure(),1000.5,1e-6);

		TS_ASSERT_DELTA(p_vessel->GetFlowRate(),(3393-1000.5)/(segments.size()*impedance),1e-6);

		for (unsigned i = 0; i < segments.size(); i++)
		{
			TS_ASSERT_DELTA(segments[i]->GetFlowProperties()->GetFlowRate(),(3393-1000.5)/(segments.size()*impedance),1e-6);
		}

	}

    void TestFlowThroughMultipleVessels() throw(Exception)
    {

        // Make some nodes
        std::vector<ChastePoint<3> > points;
        points.push_back(ChastePoint<3>(1.0, 0, 0));
        points.push_back(ChastePoint<3>(2.0, 0, 0));
        points.push_back(ChastePoint<3>(3.0, 0, 0));

        std::vector<NodePtr3> nodes;
        for(unsigned i=0; i < points.size(); i++)
        {
            nodes.push_back(NodePtr3 (VascularNode<3>::Create(points[i])));
        }

        SegmentPtr3 p_segment1(CaVesselSegment<3>::Create(nodes[0], nodes[1]));
        SegmentPtr3 p_segment2(CaVesselSegment<3>::Create(nodes[1], nodes[2]));

        VesselPtr3 p_vessel1(CaVessel<3>::Create(p_segment1));
        VesselPtr3 p_vessel2(CaVessel<3>::Create(p_segment2));

        // Generate the network
        boost::shared_ptr<CaVascularNetwork<3> > p_vascular_network(new CaVascularNetwork<3>());

        p_vascular_network->AddVessel(p_vessel1);
        p_vascular_network->AddVessel(p_vessel2);

        double impedance = 1.e14;
        p_segment1->GetFlowProperties()->SetImpedance(impedance);
        p_vascular_network->SetSegmentProperties(p_segment1);

        nodes[0]->GetFlowProperties()->SetIsInputNode(true);
        nodes[0]->GetFlowProperties()->SetPressure(3393);
        nodes[2]->GetFlowProperties()->SetIsOutputNode(true);
        nodes[2]->GetFlowProperties()->SetPressure(1000.5);

        SimpleFlowSolver<3> solver;
        solver.SetUp(p_vascular_network);
        solver.Implement();

        TS_ASSERT_DELTA(p_vessel1->GetStartNode()->GetFlowProperties()->GetPressure(),3393,1e-6);
        TS_ASSERT_DELTA(p_vessel2->GetEndNode()->GetFlowProperties()->GetPressure(),1000.5,1e-6);
        TS_ASSERT_DELTA(nodes[1]->GetFlowProperties()->GetPressure(),(3393 + 1000.5) / 2.0, 1e-6);
        TS_ASSERT_DELTA(p_vessel1->GetFlowRate(),(3393-1000.5)/(2.0 * impedance),1e-6);

    }

	void TestFlowThroughBifurcation() throw(Exception)
	{

		// Make some nodes
		std::vector<ChastePoint<3> > points;
		points.push_back(ChastePoint<3>(0, 0, 0)); // input
		points.push_back(ChastePoint<3>(0, 1.0, 0)); // input
		points.push_back(ChastePoint<3>(1.0, 0.5, 0)); // bifurcation
		points.push_back(ChastePoint<3>(1.0, 1.0, 0));
		points.push_back(ChastePoint<3>(2.0, 1.0, 0));
		points.push_back(ChastePoint<3>(3.0, 1.0, 0));
		points.push_back(ChastePoint<3>(4.0, 1.0, 0)); // output

		std::vector<NodePtr3> nodes;
		for(unsigned i=0; i < points.size(); i++)
		{
			nodes.push_back(NodePtr3 (VascularNode<3>::Create(points[i])));
		}

		SegmentPtr3 p_segment1(CaVesselSegment<3>::Create(nodes[0], nodes[2]));
		SegmentPtr3 p_segment2(CaVesselSegment<3>::Create(nodes[1], nodes[2]));
		SegmentPtr3 p_segment3(CaVesselSegment<3>::Create(nodes[3], nodes[2]));
		SegmentPtr3 p_segment4(CaVesselSegment<3>::Create(nodes[4], nodes[3]));
		SegmentPtr3 p_segment5(CaVesselSegment<3>::Create(nodes[5], nodes[4]));
		SegmentPtr3 p_segment6(CaVesselSegment<3>::Create(nodes[6], nodes[5]));

		VesselPtr3 p_vessel1(CaVessel<3>::Create(p_segment1));
		VesselPtr3 p_vessel2(CaVessel<3>::Create(p_segment2));
		VesselPtr3 p_vessel3(CaVessel<3>::Create(p_segment3));
		VesselPtr3 p_vessel4(CaVessel<3>::Create(p_segment4));
		VesselPtr3 p_vessel5(CaVessel<3>::Create(p_segment5));
		VesselPtr3 p_vessel6(CaVessel<3>::Create(p_segment6));

		std::vector<VesselPtr3> vessels;
		vessels.push_back(p_vessel1); // lower input vessel
		vessels.push_back(p_vessel2); // upper input vessel
		vessels.push_back(p_vessel3);
		vessels.push_back(p_vessel4);
		vessels.push_back(p_vessel5);
		vessels.push_back(p_vessel6);

		// Generate the network
		boost::shared_ptr<CaVascularNetwork<3> > p_vascular_network(new CaVascularNetwork<3>());

		p_vascular_network->AddVessels(vessels);

        double impedance = 1.e14;;
        p_segment1->GetFlowProperties()->SetImpedance(impedance);
        p_vascular_network->SetSegmentProperties(p_segment1);

		nodes[0]->GetFlowProperties()->SetIsInputNode(true);
		nodes[0]->GetFlowProperties()->SetPressure(3393);

		nodes[1]->GetFlowProperties()->SetIsInputNode(true);
		nodes[1]->GetFlowProperties()->SetPressure(3393);

		nodes[6]->GetFlowProperties()->SetIsOutputNode(true);
		nodes[6]->GetFlowProperties()->SetPressure(1000.5);

		SimpleFlowSolver<3> solver;
		solver.SetUp(p_vascular_network);
		solver.Implement();

		TS_ASSERT_DELTA(nodes[0]->GetFlowProperties()->GetPressure(),3393,1e-6);
		TS_ASSERT_DELTA(nodes[1]->GetFlowProperties()->GetPressure(),3393,1e-6);
		TS_ASSERT_DELTA(nodes[2]->GetFlowProperties()->GetPressure(),(2.0*3393.0/10.0 + 1000.5/40.0)/(1.0/40.0 + 2.0/10.0),1e-6);
		TS_ASSERT_DELTA(nodes[6]->GetFlowProperties()->GetPressure(),1000.5,1e-6);

		TS_ASSERT_DELTA(vessels[0]->GetFlowRate(),(3393-nodes[2]->GetFlowProperties()->GetPressure())/impedance,1e-6);
		TS_ASSERT_DELTA(vessels[1]->GetFlowRate(),(3393-nodes[2]->GetFlowProperties()->GetPressure())/impedance,1e-6);
		TS_ASSERT_DELTA(vessels[5]->GetFlowRate(),-(nodes[2]->GetFlowProperties()->GetPressure()-1000.5)/(4.0 * impedance),1e-6);

		TS_ASSERT_DELTA(p_segment1->GetFlowProperties()->GetFlowRate(),(3393-nodes[2]->GetFlowProperties()->GetPressure())/impedance,1e-6);
		TS_ASSERT_DELTA(p_segment2->GetFlowProperties()->GetFlowRate(),(3393-nodes[2]->GetFlowProperties()->GetPressure())/impedance,1e-6);
		TS_ASSERT_DELTA(p_segment6->GetFlowProperties()->GetFlowRate(),-(nodes[2]->GetFlowProperties()->GetPressure()-1000.5)/(4.0 * impedance),1e-6);

		double kirchoff_residual = vessels[0]->GetFlowRate() + vessels[1]->GetFlowRate() +
									vessels[5]->GetFlowRate();

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

        double impedance = 1.e14;
        p_segment1->GetFlowProperties()->SetImpedance(impedance);
        p_vascular_network->SetSegmentProperties(p_segment1);

		nodes[0]->GetFlowProperties()->SetIsInputNode(true);
		nodes[0]->GetFlowProperties()->SetPressure(3393);

		nodes[1]->GetFlowProperties()->SetIsInputNode(true);
		nodes[1]->GetFlowProperties()->SetPressure(3393);

		nodes[3]->GetFlowProperties()->SetIsOutputNode(true);
		nodes[3]->GetFlowProperties()->SetPressure(1000.5);

		SimpleFlowSolver<3> solver;
		solver.SetUp(p_vascular_network);
		solver.Implement();

		TS_ASSERT_DELTA(nodes[0]->GetFlowProperties()->GetPressure(),3393,1e-6);
		TS_ASSERT_DELTA(nodes[1]->GetFlowProperties()->GetPressure(),3393,1e-6);
		TS_ASSERT_DELTA(nodes[2]->GetFlowProperties()->GetPressure(),(2*3393 + 1000.5)/3,1e-6);
		TS_ASSERT_DELTA(nodes[3]->GetFlowProperties()->GetPressure(),1000.5,1e-6);

		TS_ASSERT_DELTA(vessels[0]->GetFlowRate(),(3393-nodes[2]->GetFlowProperties()->GetPressure())/impedance,1e-6);
		TS_ASSERT_DELTA(vessels[1]->GetFlowRate(),(3393-nodes[2]->GetFlowProperties()->GetPressure())/impedance,1e-6);
		TS_ASSERT_DELTA(vessels[2]->GetFlowRate(),(nodes[2]->GetFlowProperties()->GetPressure()-1000.5)/impedance,1e-6);

		TS_ASSERT_DELTA(p_segment1->GetFlowProperties()->GetFlowRate(),(3393-nodes[2]->GetFlowProperties()->GetPressure())/impedance,1e-6);
		TS_ASSERT_DELTA(p_segment2->GetFlowProperties()->GetFlowRate(),(3393-nodes[2]->GetFlowProperties()->GetPressure())/impedance,1e-6);
		TS_ASSERT_DELTA(p_segment3->GetFlowProperties()->GetFlowRate(),(nodes[2]->GetFlowProperties()->GetPressure()-1000.5)/impedance,1e-6);

		double kirchoff_residual = vessels[0]->GetFlowRate() + vessels[1]->GetFlowRate() -
				vessels[2]->GetFlowRate();

		TS_ASSERT_DELTA(kirchoff_residual,0,1e-6);

	}

	void TestFlowThroughHexagonalNetwork() throw(Exception)
    {
		// Specify the network dimensions
		double vessel_length = 80.0;

		// Generate the network
		VasculatureGenerator<2> vascular_network_generator;
		boost::shared_ptr<CaVascularNetwork<2> > vascular_network = vascular_network_generator.GenerateHexagonalNetwork(1000,1000,vessel_length);

        // Make some nodes
        std::vector<ChastePoint<2> > points;
        points.push_back(ChastePoint<2>(0, 0, 0)); // input
        points.push_back(ChastePoint<2>(0, 1.0, 0)); // input

        std::vector<NodePtr2> nodes;
        for(unsigned i=0; i < points.size(); i++)
        {
            nodes.push_back(NodePtr2 (VascularNode<2>::Create(points[i])));
        }

		SegmentPtr2 p_segment1(CaVesselSegment<2>::Create(nodes[0], nodes[1]));

        double impedance = 1.e14;
        p_segment1->GetFlowProperties()->SetImpedance(impedance);
        vascular_network->SetSegmentProperties(p_segment1);

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
					(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
					(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3393);
				}
			}
			if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
			{
				if((*vessel_iterator)->GetEndNode()->GetLocation()[1] >  y_middle)
				{
					(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
					(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3393);
				}
			}
			if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
			{
				if((*vessel_iterator)->GetStartNode()->GetLocation()[1] <=  y_middle)
				{
					(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
					(*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(1993);
				}
			}
			if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
			{
				if((*vessel_iterator)->GetEndNode()->GetLocation()[1] <=  y_middle)
				{
					(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
					(*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(1993);
				}
			}

		}

		SimpleFlowSolver<2> solver;
		solver.SetUp(vascular_network);
		solver.Implement();

		// Write the network to file
		OutputFileHandler output_file_handler("TestSimpleFlowSolver", false);
		std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexagonalVesselNetwork.vtp");
		vascular_network->Write(output_filename);
    }

};

#endif /*TESTSIMPLEFLOWSOLVER_HPP_*/
