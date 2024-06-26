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

#ifndef TESTALARCON03METABOLICSTIMULUSCALCULATOR_HPP
#define TESTALARCON03METABOLICSTIMULUSCALCULATOR_HPP

#include <cxxtest/TestSuite.h>
#include <math.h>
#include <boost/shared_ptr.hpp>
#include "MetabolicStimulusCalculator.hpp"
#include "VesselNetwork.hpp"
#include "UnitCollection.hpp"
#include "OutputFileHandler.hpp"

#include "FakePetscSetup.hpp"

class TestMetabolicStimulusCalculator : public CxxTest::TestSuite
{

public:

    void TestCalculator()
    {
        std::vector<boost::shared_ptr<VesselNode<3> > > nodes;
        nodes.push_back(VesselNode<3>::Create(0));
        nodes.push_back(VesselNode<3>::Create(100));
        double pressure = 3933.0;
        nodes[0]->GetFlowProperties()->SetPressure(pressure*unit::pascals);
        nodes[1]->GetFlowProperties()->SetPressure(pressure*unit::pascals);

        double flow_rate = 10.0;
        double haematocrit_level = 0.45;
        boost::shared_ptr<Vessel<3> > p_vessel(Vessel<3>::Create(VesselSegment<3>::Create(nodes[0], nodes[1])));
        boost::shared_ptr<VesselNetwork<3> > p_vascular_network = VesselNetwork<3>::Create();
        p_vascular_network->AddVessel(p_vessel);
        p_vessel->GetSegments()[0]->GetFlowProperties()->SetFlowRate(flow_rate*unit::metre_cubed_per_second);
        p_vessel->GetSegments()[0]->GetFlowProperties()->SetHaematocrit(haematocrit_level);

        boost::shared_ptr<MetabolicStimulusCalculator<3> > calculator(new MetabolicStimulusCalculator<3>());
        calculator->SetVesselNetwork(p_vascular_network);
        double Q_ref = calculator->GetQRef()/unit::metre_cubed_per_second;
        double k_m = calculator->GetKm()/unit::per_second;
        double MaxStimulus = calculator->GetMaxStimulus()/unit::per_second;
        double expected_metabolic_stimulus = k_m * log10(Q_ref / (flow_rate * haematocrit_level) + 1.0);
        TS_ASSERT_DELTA(p_vessel->GetSegments()[0]->GetFlowProperties()->GetGrowthStimulus()/unit::per_second, expected_metabolic_stimulus, 1e-6);

        p_vessel->GetSegments()[0]->GetFlowProperties()->SetFlowRate(0.0 * unit::metre_cubed_per_second);
        calculator->Calculate();
        expected_metabolic_stimulus = 0;
        TS_ASSERT_DELTA(p_vessel->GetSegments()[0]->GetFlowProperties()->GetGrowthStimulus()/unit::per_second, expected_metabolic_stimulus, 1e-6);

        p_vessel->GetSegments()[0]->GetFlowProperties()->SetFlowRate(flow_rate * unit::metre_cubed_per_second);
        p_vessel->GetSegments()[0]->GetFlowProperties()->SetHaematocrit(0.0);
        calculator->Calculate();
        expected_metabolic_stimulus = MaxStimulus;
        TS_ASSERT_DELTA(p_vessel->GetSegments()[0]->GetFlowProperties()->GetGrowthStimulus()/unit::per_second, expected_metabolic_stimulus, 1e-6);
    }
};

#endif // TESTALARCONHAEMATOCRITSOLVER_HPP
