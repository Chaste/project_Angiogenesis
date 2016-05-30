//
//  TestAlarcon03MechanicalStimulusCalculator.hpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 03/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#ifndef TestAlarcon03MechanicalStimulusCalculator_hpp
#define TestAlarcon03MechanicalStimulusCalculator_hpp

#include <cxxtest/TestSuite.h>
#include <boost/shared_ptr.hpp>
#include <math.h>
#include "VascularNetwork.hpp"
#include "Alarcon03MechanicalStimulusCalculator.hpp"

class TestAlarcon03MechanicalStimulusCalculator : public CxxTest::TestSuite
{

public:

    void TestSingleVessel()
    {
        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > nodes;
        nodes.push_back(VascularNode<3>::Create(0));
        nodes.push_back(VascularNode<3>::Create(100));
        double pressure = 3933.0;
        nodes[0]->GetFlowProperties()->SetPressure(pressure);
        nodes[1]->GetFlowProperties()->SetPressure(pressure);

        boost::shared_ptr<Vessel<3> > p_vessel(Vessel<3>::Create(VesselSegment<3>::Create(nodes[0], nodes[1])));
        boost::shared_ptr<VascularNetwork<3> > p_vascular_network = VascularNetwork<3>::Create();
        p_vascular_network->AddVessel(p_vessel);
        double wall_shear_stress = 25.0;
        p_vessel->GetSegments()[0]->GetFlowProperties()->SetWallShearStress(wall_shear_stress);

        boost::shared_ptr<Alarcon03MechanicalStimulusCalculator<3> > calculator(new Alarcon03MechanicalStimulusCalculator<3>());
        calculator->Calculate(p_vascular_network);

        // convert pressure to mmhg
        double converted_pressure = pressure * 760 / (1.01 * pow(10.0, 5));
        double Tau_P = 0.1 * (100.0 - 86.0 * pow(exp(-5.0 * log10(log10(converted_pressure))), 5.4));
        double expected_mechanical_stimulus = log10((wall_shear_stress + calculator->GetTauRef()) / Tau_P);
        TS_ASSERT_DELTA(p_vessel->GetSegments()[0]->GetFlowProperties()->GetMechanicalStimulus(), expected_mechanical_stimulus, 1e-6);
    }

    void TestAlarcon03MechanicalStimulusVsPressure()
    {
        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > nodes;
        nodes.push_back(VascularNode<3>::Create(0));
        nodes.push_back(VascularNode<3>::Create(100));
        double pressure = 3933.0;
        nodes[0]->GetFlowProperties()->SetPressure(pressure);
        nodes[1]->GetFlowProperties()->SetPressure(pressure);

        boost::shared_ptr<Vessel<3> > p_vessel(Vessel<3>::Create(VesselSegment<3>::Create(nodes[0], nodes[1])));
        boost::shared_ptr<VascularNetwork<3> > p_vascular_network = VascularNetwork<3>::Create();
        p_vascular_network->AddVessel(p_vessel);
        double wall_shear_stress = 25.0;
        p_vessel->GetSegments()[0]->GetFlowProperties()->SetWallShearStress(wall_shear_stress);

        boost::shared_ptr<Alarcon03MechanicalStimulusCalculator<3> > calculator(new Alarcon03MechanicalStimulusCalculator<3>());
        calculator->Calculate(p_vascular_network);

        OutputFileHandler output_file(
                "Vasculature/StructuralAdaptationAlgorithm/TestAlarcon03MechanicalStimulusCalculator", true);

        std::string filename = output_file.GetOutputDirectoryFullPath() + "ShearStressVsPressure.dat";
        std::ofstream out(filename.c_str());
        out << "# Pressure (mmHg) # TauP (dyne/cm^2)\n";

        for (int mmHgPressure = 1; mmHgPressure < 101; mmHgPressure++)
        {
            nodes[0]->GetFlowProperties()->SetPressure(double(mmHgPressure) * (1.01 * pow(10.0, 5) / 760));
            nodes[1]->GetFlowProperties()->SetPressure(double(mmHgPressure) * (1.01 * pow(10.0, 5) / 760));
            calculator->Calculate(p_vascular_network);
            out << mmHgPressure << " " << calculator->GetTauP() / 0.1 << "\n";
        }
        out.close();
    }

};

#endif
