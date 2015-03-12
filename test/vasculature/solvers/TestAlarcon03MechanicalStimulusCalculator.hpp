//
//  TestAlarcon03MechanicalStimulusCalculator.hpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 03/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#ifndef VascularTumourGrowthModellingFramework_TestAlarcon03MechanicalStimulusCalculator_hpp
#define VascularTumourGrowthModellingFramework_TestAlarcon03MechanicalStimulusCalculator_hpp

#include <cxxtest/TestSuite.h>
#include "Alarcon03MechanicalStimulusCalculator.hpp"
#include <boost/shared_ptr.hpp>
#include "CaVascularNetwork.hpp"
#include <math.h>
#include <fstream>
#include "Debug.hpp"

class TestAlarcon03MechanicalStimulusCalculator : public CxxTest::TestSuite
{
	typedef boost::shared_ptr<VascularNode<2> > NodePtr2;
	typedef boost::shared_ptr<VascularNode<3> > NodePtr3;
	typedef boost::shared_ptr<CaVesselSegment<2> > SegmentPtr2;
	typedef boost::shared_ptr<CaVesselSegment<3> > SegmentPtr3;
	typedef boost::shared_ptr<CaVessel<2> > VesselPtr2;
	typedef boost::shared_ptr<CaVessel<3> > VesselPtr3;

public:
        
    void testAlarcon03MechanicalStimulusCalculator()
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

		VasculatureData node_data;
		double pressure = 3993.0;
		node_data.SetData("Pressure", pressure);

		p_vascular_network->SetNodeData(node_data);

		VasculatureData segment_data;
		double wall_shear_stress = 25.0;
		segment_data.SetData("Wall Shear Stress", wall_shear_stress);

		p_vascular_network->SetSegmentData(segment_data);

        boost::shared_ptr<Alarcon03MechanicalStimulusCalculator<3> > calculator(new Alarcon03MechanicalStimulusCalculator<3>());
      
        calculator->Calculate(p_vascular_network);

        // convert pressure to mmhg
        double converted_pressure = pressure * 760/(1.01*pow(10.0,5));

    	double Tau_P = 0.1*(100.0 - 86.0*pow(exp(-5.0*log10(log10(converted_pressure))), 5.4));

    	double expected_mechanical_stimulus = log10((wall_shear_stress + calculator->GetTauRef())/Tau_P);

        TS_ASSERT_DELTA(p_segment->template GetData<double>("Mechanical Stimulus"), expected_mechanical_stimulus,1e-6);
    }
    
    void testAlarcon03MechanicalStimulusVsPressure()
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

		VasculatureData node_data;
		double pressure = 3993.0;
		node_data.SetData("Pressure", pressure);

		p_vascular_network->SetNodeData(node_data);

		VasculatureData segment_data;
		double wall_shear_stress = 25.0;
		segment_data.SetData("Wall Shear Stress", wall_shear_stress);

		p_vascular_network->SetSegmentData(segment_data);

        boost::shared_ptr<Alarcon03MechanicalStimulusCalculator<3> > calculator(new Alarcon03MechanicalStimulusCalculator<3>());

        calculator->Calculate(p_vascular_network);

        OutputFileHandler output_file("Vasculature/StructuralAdaptationAlgorithm/TestAlarcon03MechanicalStimulusCalculator", true);

        std::string filename = output_file.GetOutputDirectoryFullPath() + "ShearStressVsPressure.dat";
        std::ofstream out(filename.c_str());
        out << "# Pressure (mmHg) # TauP (dyne/cm^2)\n";

        for (int mmHgPressure = 1; mmHgPressure < 101; mmHgPressure++)
        {
        	p_segment->GetNode(0)->SetData("Pressure", double(mmHgPressure)*(1.01*pow(10.0,5)/760));
        	p_segment->GetNode(1)->SetData("Pressure", double(mmHgPressure)*(1.01*pow(10.0,5)/760));


        	calculator->Calculate(p_vascular_network);

            out << mmHgPressure << " " << calculator->GetTauP()/0.1 << "\n";

        }

        out.close();
    }
    
};

#endif
