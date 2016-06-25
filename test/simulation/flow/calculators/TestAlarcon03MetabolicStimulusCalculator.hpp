//
//  TestAlarcon03MetabolicStimulusCalculator.hpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 03/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#ifndef TestAlarcon03MetabolicStimulusCalculator_hpp
#define TestAlarcon03MetabolicStimulusCalculator_hpp

#include <cxxtest/TestSuite.h>
#include <math.h>
#include <boost/shared_ptr.hpp>
#include "MetabolicStimulusCalculator.hpp"
#include "VesselNetwork.hpp"
#include "UnitCollection.hpp"
#include "OutputFileHandler.hpp"

class TestMetabolicStimulusCalculator : public CxxTest::TestSuite
{

public:

    void TestCalculator()
    {
        std::vector<boost::shared_ptr<VesselNode<3> > > nodes;
        nodes.push_back(VesselNode<3>::Create(0));
        nodes.push_back(VesselNode<3>::Create(100));
        double pressure = 3933.0;
        nodes[0]->GetFlowProperties()->SetPressure(pressure);
        nodes[1]->GetFlowProperties()->SetPressure(pressure);

        double flow_rate = 10.0;
        double haematocrit_level = 0.45;
        boost::shared_ptr<Vessel<3> > p_vessel(Vessel<3>::Create(VesselSegment<3>::Create(nodes[0], nodes[1])));
        boost::shared_ptr<VesselNetwork<3> > p_vascular_network = VesselNetwork<3>::Create();
        p_vascular_network->AddVessel(p_vessel);
        p_vessel->GetSegments()[0]->GetFlowProperties()->SetFlowRate(flow_rate);
        p_vessel->GetSegments()[0]->GetFlowProperties()->SetHaematocrit(haematocrit_level);

        boost::shared_ptr<MetabolicStimulusCalculator<3> > calculator(new MetabolicStimulusCalculator<3>());
        calculator->SetVesselNetwork(p_vascular_network);
        double Q_ref = calculator->GetQRef()/unit::unit_flow_rate;
        double k_m = calculator->GetKm()/unit::reciprocal_seconds;
        double MaxStimulus = calculator->GetMaxStimulus()/unit::reciprocal_seconds;
        double expected_metabolic_stimulus = k_m * log10(Q_ref / (flow_rate * haematocrit_level) + 1.0);
        TS_ASSERT_DELTA(p_vessel->GetSegments()[0]->GetFlowProperties()->GetGrowthStimulus(), expected_metabolic_stimulus, 1e-6);

        p_vessel->GetSegments()[0]->GetFlowProperties()->SetFlowRate(0.0);
        calculator->Calculate();
        expected_metabolic_stimulus = 0;
        TS_ASSERT_DELTA(p_vessel->GetSegments()[0]->GetFlowProperties()->GetGrowthStimulus(), expected_metabolic_stimulus, 1e-6);

        p_vessel->GetSegments()[0]->GetFlowProperties()->SetFlowRate(flow_rate);
        p_vessel->GetSegments()[0]->GetFlowProperties()->SetHaematocrit(0.0);
        calculator->Calculate();
        expected_metabolic_stimulus = MaxStimulus;
        TS_ASSERT_DELTA(p_vessel->GetSegments()[0]->GetFlowProperties()->GetGrowthStimulus(), expected_metabolic_stimulus, 1e-6);
    }
};

#endif
