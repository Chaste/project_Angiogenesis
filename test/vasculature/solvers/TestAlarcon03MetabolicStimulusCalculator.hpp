//
//  TestAlarcon03MetabolicStimulusCalculator.hpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 03/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#ifndef VascularTumourGrowthModellingFramework_TestAlarcon03MetabolicStimulusCalculator_hpp
#define VascularTumourGrowthModellingFramework_TestAlarcon03MetabolicStimulusCalculator_hpp

#include <cxxtest/TestSuite.h>
#include "Alarcon03MetabolicStimulusCalculator.hpp"
#include <boost/shared_ptr.hpp>
#include "CaVascularNetwork.hpp"
#include <math.h>
#include <fstream>
#include "Debug.hpp"

class TestAlarcon03MetabolicStimulusCalculator : public CxxTest::TestSuite
{
	typedef boost::shared_ptr<VascularNode<2> > NodePtr2;
	typedef boost::shared_ptr<VascularNode<3> > NodePtr3;
	typedef boost::shared_ptr<CaVesselSegment<2> > SegmentPtr2;
	typedef boost::shared_ptr<CaVesselSegment<3> > SegmentPtr3;
	typedef boost::shared_ptr<CaVessel<2> > VesselPtr2;
	typedef boost::shared_ptr<CaVessel<3> > VesselPtr3;

public:

	void testAlarcon03MetabolicStimulusCalculator()
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

		VasculatureData segment_data;
		double flow_rate = 10.0;
		double haematocrit_level = 0.45;
		segment_data.SetData("Flow Rate", flow_rate);
		segment_data.SetData("Haematocrit Level",haematocrit_level);

		p_vascular_network->SetSegmentData(segment_data);

		boost::shared_ptr<Alarcon03MetabolicStimulusCalculator<3> > calculator(new Alarcon03MetabolicStimulusCalculator<3>());

		calculator->Calculate(p_vascular_network);

		double Q_ref = calculator->GetQRef();
		double k_m = calculator->GetKm();
		double MaxStimulus = calculator->GetMaxStimulus();

		double expected_metabolic_stimulus =  k_m*log10(Q_ref/(flow_rate*haematocrit_level) + 1.0);

		TS_ASSERT_DELTA(p_segment->template GetData<double>("Metabolic Stimulus"), expected_metabolic_stimulus,1e-6);

		p_segment->SetData("Flow Rate", 0.0);

		calculator->Calculate(p_vascular_network);

		expected_metabolic_stimulus = 0;

		TS_ASSERT_DELTA(p_segment->template GetData<double>("Metabolic Stimulus"), expected_metabolic_stimulus,1e-6);

		p_segment->SetData("Flow Rate", flow_rate);
		p_segment->SetData("Haematocrit Level",0.0);

		calculator->Calculate(p_vascular_network);

		expected_metabolic_stimulus = MaxStimulus;

		TS_ASSERT_DELTA(p_segment->template GetData<double>("Metabolic Stimulus"), expected_metabolic_stimulus,1e-6);
	}

};

#endif
