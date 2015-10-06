//
//  TestSimpleStructuralAdaptationSolver.hpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 03/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#ifndef TESTANASTA_HPP
#define TESTANASTA_HPP

#include <cxxtest/TestSuite.h>
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VasculatureGenerator.hpp"
#include "VascularNode.hpp"
#include "CaVesselSegment.hpp"
#include "CaVessel.hpp"
#include "CaVascularNetwork.hpp"
#include "Part.hpp"
#include "FiniteDifferenceSolver.hpp"
#include "AbstractAngiogenesisSolver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CaVesselSegment.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SimpleFlowSolver.hpp"
#include "OffLatticePrwGrowthDirectionModifier.hpp"
#include "AbstractSproutingRule.hpp"
#include "OffLatticeRandomNormalSproutingRule.hpp"
#include "OffLatticeTipAttractionGrowthDirectionModifier.hpp"
#include "RandomNumberGenerator.hpp"

class TestAnastamosis : public AbstractCellBasedTestSuite
{

public:

    void TestParallelConterflow() throw(Exception)
    {
        RandomNumberGenerator::Instance()->Reseed(123456);

        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        for(unsigned idx=0; idx<41; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 10.0, 0.0));
        }

        std::vector<boost::shared_ptr<VascularNode<3> > > top_nodes;
        for(unsigned idx=0; idx<41; idx++)
        {
            top_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 450.0, 0.0));
        }

        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVessel<3> > p_vessel2 = CaVessel<3>::Create(top_nodes);
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);
        p_network->SetSegmentRadii(10.0);

        for(unsigned idx=1; idx<40; idx+=2)
        {
            p_network->FormSprout(ChastePoint<3>(double(idx)*10, 10.0, 0.0), ChastePoint<3>(double(idx)*10, 20.0, 0.0));
            p_network->FormSprout(ChastePoint<3>(double(idx)*10, 450.0, 0.0), ChastePoint<3>(double(idx)*10, 440.0, 0.0));
        }

        boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> > p_grow_direction_modifier = OffLatticePrwGrowthDirectionModifier<3>::Create();
        boost::shared_ptr<OffLatticeTipAttractionGrowthDirectionModifier<3> > p_grow_direction_modifier2 = OffLatticeTipAttractionGrowthDirectionModifier<3>::Create();
        p_grow_direction_modifier2->SetNetwork(p_network);
        boost::shared_ptr<OffLatticeRandomNormalSproutingRule<3> > p_sprouting_rule = OffLatticeRandomNormalSproutingRule<3>::Create();
        p_sprouting_rule->SetSproutingProbability(0.005);

        OutputFileHandler output_file_handler("TestAnastamosis/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(40, 40);
        AbstractAngiogenesisSolver<3> angiogenesis_solver(p_network);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier2);
        angiogenesis_solver.SetSproutingRule(p_sprouting_rule);
        angiogenesis_solver.SetAnastamosisRadius(4.0);
        angiogenesis_solver.SetEndTime(40.0);

        angiogenesis_solver.SetOutputDirectory(output_directory);
        angiogenesis_solver.Run();
    }

};

#endif
