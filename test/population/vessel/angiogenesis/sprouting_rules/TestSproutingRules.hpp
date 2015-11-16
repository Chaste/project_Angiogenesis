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
#include "../../../../../src/population/vessel/angiogenesis/sprouting_rules/OffLatticeSproutingRule.hpp"
#include "OffLatticePrwGrowthDirectionModifier.hpp"
#include "OffLatticeTipAttractionGrowthDirectionModifier.hpp"
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
#include "FlowSolver.hpp"
#include "AbstractSproutingRule.hpp"
#include "RandomNumberGenerator.hpp"

class TestSproutingRules : public AbstractCellBasedTestSuite
{

public:

    void TestAbstractSprouting() throw(Exception)
    {
        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        for(unsigned idx=0; idx<9; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 10.0, 0.0));
        }
        bottom_nodes[0]->GetFlowProperties()->SetIsInputNode(true);
        bottom_nodes[0]->GetFlowProperties()->SetPressure(3000);
        bottom_nodes[8]->GetFlowProperties()->SetIsOutputNode(true);
        bottom_nodes[8]->GetFlowProperties()->SetPressure(1000);

        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);
        p_network->SetSegmentRadii(10.0);
        p_network->GetVesselSegments()[0]->GetFlowProperties()->SetViscosity(1.e-9);
        p_network->CopySegmentFlowProperties();

        boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> > p_grow_direction_modifier = OffLatticePrwGrowthDirectionModifier<3>::Create();
        boost::shared_ptr<AbstractSproutingRule<3> > p_sprouting_rule = AbstractSproutingRule<3>::Create();
        p_sprouting_rule->SetSproutingProbability(0.1);

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AngiogenesisSolver<3> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetOutputDirectory("TestSproutingRules/Abstract");
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.SetSproutingRule(p_sprouting_rule);
        angiogenesis_solver.Run();
    }

    void TestOffLatticeSprouting() throw(Exception)
    {
        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        for(unsigned idx=0; idx<9; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 10.0, 0.0));
        }
        bottom_nodes[0]->GetFlowProperties()->SetIsInputNode(true);
        bottom_nodes[0]->GetFlowProperties()->SetPressure(3000);
        bottom_nodes[8]->GetFlowProperties()->SetIsOutputNode(true);
        bottom_nodes[8]->GetFlowProperties()->SetPressure(1000);

        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);
        p_network->SetSegmentRadii(10.0);
        p_network->GetVesselSegments()[0]->GetFlowProperties()->SetViscosity(1.e-9);
        p_network->CopySegmentFlowProperties();

        boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> > p_grow_direction_modifier = OffLatticePrwGrowthDirectionModifier<3>::Create();
        boost::shared_ptr<OffLatticeRandomNormalSproutingRule<3> > p_sprouting_rule = OffLatticeRandomNormalSproutingRule<3>::Create();
        p_sprouting_rule->SetSproutingProbability(0.1);

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AngiogenesisSolver<3> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetOutputDirectory("TestSproutingRules/OffLattice");
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.SetSproutingRule(p_sprouting_rule);
        angiogenesis_solver.Run();
    }

    void TestOffLatticeSproutingWithAttraction() throw(Exception)
    {
        RandomNumberGenerator::Instance()->Reseed(15565);

        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        for(unsigned idx=0; idx<9; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*100, 10.0, 0.0));
        }
        bottom_nodes[0]->GetFlowProperties()->SetIsInputNode(true);
        bottom_nodes[0]->GetFlowProperties()->SetPressure(3000);
        bottom_nodes[8]->GetFlowProperties()->SetIsOutputNode(true);
        bottom_nodes[8]->GetFlowProperties()->SetPressure(1000);

        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);
        p_network->SetSegmentRadii(10.0);
        p_network->GetVesselSegments()[0]->GetFlowProperties()->SetViscosity(1.e-9);
        p_network->CopySegmentFlowProperties();

        boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> > p_grow_direction_modifier = OffLatticePrwGrowthDirectionModifier<3>::Create();
        boost::shared_ptr<OffLatticeTipAttractionGrowthDirectionModifier<3> > p_grow_direction_modifier2 = OffLatticeTipAttractionGrowthDirectionModifier<3>::Create();
        p_grow_direction_modifier2->SetNetwork(p_network);

        boost::shared_ptr<OffLatticeRandomNormalSproutingRule<3> > p_sprouting_rule = OffLatticeRandomNormalSproutingRule<3>::Create();
        p_sprouting_rule->SetSproutingProbability(0.1);

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(20, 20);
        AngiogenesisSolver<3> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetOutputDirectory("TestSproutingRules/OffLatticeAttraction");
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier2);
        angiogenesis_solver.SetSproutingRule(p_sprouting_rule);
        angiogenesis_solver.SetAnastamosisRadius(5.0);
        angiogenesis_solver.Run();
    }
};

#endif
