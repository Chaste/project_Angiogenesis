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
#include "Alarcon03MetabolicStimulusCalculator.hpp"
#include "CaVascularNetwork.hpp"

class TestAlarcon03MetabolicStimulusCalculator : public CxxTest::TestSuite
{

public:

    void testAlarcon03MetabolicStimulusCalculator()
    {
        std::vector<boost::shared_ptr<VascularNode<3> > > nodes;
        nodes.push_back(VascularNode<3>::Create(0));
        nodes.push_back(VascularNode<3>::Create(100));
        double pressure = 3933.0;
        nodes[0]->GetFlowProperties()->SetPressure(pressure);
        nodes[1]->GetFlowProperties()->SetPressure(pressure);

        double flow_rate = 10.0;
        double haematocrit_level = 0.45;
        boost::shared_ptr<CaVessel<3> > p_vessel(CaVessel<3>::Create(CaVesselSegment<3>::Create(nodes[0], nodes[1])));
        boost::shared_ptr<CaVascularNetwork<3> > p_vascular_network = CaVascularNetwork<3>::Create();
        p_vascular_network->AddVessel(p_vessel);
        p_vessel->GetSegments()[0]->GetFlowProperties()->SetFlowRate(flow_rate);
        p_vessel->GetSegments()[0]->GetFlowProperties()->SetHaematocrit(haematocrit_level);

        boost::shared_ptr<Alarcon03MetabolicStimulusCalculator<3> > calculator(
                new Alarcon03MetabolicStimulusCalculator<3>());

        calculator->Calculate(p_vascular_network);
        double Q_ref = calculator->GetQRef();
        double k_m = calculator->GetKm();
        double MaxStimulus = calculator->GetMaxStimulus();
        double expected_metabolic_stimulus = k_m * log10(Q_ref / (flow_rate * haematocrit_level) + 1.0);
        TS_ASSERT_DELTA(p_vessel->GetSegments()[0]->GetFlowProperties()->GetMetabolicStimulus(), expected_metabolic_stimulus, 1e-6);

        p_vessel->GetSegments()[0]->GetFlowProperties()->SetFlowRate(0.0);
        calculator->Calculate(p_vascular_network);
        expected_metabolic_stimulus = 0;
        TS_ASSERT_DELTA(p_vessel->GetSegments()[0]->GetFlowProperties()->GetMetabolicStimulus(), expected_metabolic_stimulus, 1e-6);

        p_vessel->GetSegments()[0]->GetFlowProperties()->SetFlowRate(flow_rate);
        p_vessel->GetSegments()[0]->GetFlowProperties()->SetHaematocrit(0.0);
        calculator->Calculate(p_vascular_network);
        expected_metabolic_stimulus = MaxStimulus;
        TS_ASSERT_DELTA(p_vessel->GetSegments()[0]->GetFlowProperties()->GetMetabolicStimulus(), expected_metabolic_stimulus, 1e-6);
    }
};

#endif
