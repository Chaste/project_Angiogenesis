//
//  TestSimpleStructuralAdaptationSolver.hpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 03/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#ifndef TESTSPROUTINGRULES_HPP
#define TESTSPROUTINGRULES_HPP

#include <cxxtest/TestSuite.h>
#include "../../src/angiogenesis_solvers/growth_direction_modifiers/OffLatticePrwGrowthDirectionModifier.hpp"
#include "../../src/angiogenesis_solvers/growth_direction_modifiers/OffLatticeTipAttractionGrowthDirectionModifier.hpp"
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
#include "AngiogenesisSolver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CaVesselSegment.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SimpleFlowSolver.hpp"
#include "AbstractSproutingRule.hpp"
#include "OffLatticeRandomNormalSproutingRule.hpp"
#include "RandomNumberGenerator.hpp"

class TestPrwNetworkGenerator : public AbstractCellBasedTestSuite
{

public:

    void TestGenerateNetwork() throw(Exception)
    {
        RandomNumberGenerator::Instance()->Reseed(15565);

        // Make initial vessels on the forward and rear faces of a cube
        double cube_x = 800.0;
        double cube_y = 800.0;
        double cube_z = 800.0;
        boost::shared_ptr<Part<3> > p_part = Part<3>::Create();
        p_part->AddCuboid(cube_x, cube_y, cube_z);

        unsigned num_points = 20;
        std::vector<boost::shared_ptr<CaVessel<3> > > front_vessels;
        for(unsigned idx=0; idx<num_points; idx++)
        {
            double x_position = RandomNumberGenerator::Instance()->ranf() * cube_x;
            double y_position = RandomNumberGenerator::Instance()->ranf() * cube_y;
            double z_position = 0.0;

            boost::shared_ptr<VascularNode<3> > p_start_node = VascularNode<3>::Create(x_position, y_position, z_position);
            boost::shared_ptr<VascularNode<3> > p_end_node = VascularNode<3>::Create(x_position, y_position, z_position + 40.0);
            p_end_node->SetIsMigrating(true);
            front_vessels.push_back(CaVessel<3>::Create(p_start_node, p_end_node));
        }

        std::vector<boost::shared_ptr<CaVessel<3> > > back_vessels;
        for(unsigned idx=0; idx<num_points; idx++)
        {
            double x_position = RandomNumberGenerator::Instance()->ranf() * cube_x;
            double y_position = RandomNumberGenerator::Instance()->ranf() * cube_y;
            double z_position = cube_z;

            boost::shared_ptr<VascularNode<3> > p_start_node = VascularNode<3>::Create(x_position, y_position, z_position);
            boost::shared_ptr<VascularNode<3> > p_end_node = VascularNode<3>::Create(x_position, y_position, z_position - 40.0);
            p_end_node->SetIsMigrating(true);
            back_vessels.push_back(CaVessel<3>::Create(p_start_node, p_end_node));
        }

        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessels(front_vessels);
        p_network->AddVessels(back_vessels);
        p_network->SetSegmentRadii(10.0);

        boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> > p_grow_direction_modifier = OffLatticePrwGrowthDirectionModifier<3>::Create();
        boost::shared_ptr<OffLatticeTipAttractionGrowthDirectionModifier<3> > p_grow_direction_modifier2 = OffLatticeTipAttractionGrowthDirectionModifier<3>::Create();
        p_grow_direction_modifier2->SetNetwork(p_network);

//        boost::shared_ptr<OffLatticeRandomNormalSproutingRule<3> > p_sprouting_rule = OffLatticeRandomNormalSproutingRule<3>::Create();
//        p_sprouting_rule->SetSproutingProbability(0.005);

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(100, 100);
        AngiogenesisSolver<3> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetOutputDirectory("TestPrwNetworkGenerator");
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier2);
//        angiogenesis_solver.SetSproutingRule(p_sprouting_rule);
        angiogenesis_solver.SetAnastamosisRadius(5.0);
        angiogenesis_solver.SetEndTime(100.0);
        angiogenesis_solver.SetBoundingDomain(p_part);
        angiogenesis_solver.Run();
    }
};

#endif
