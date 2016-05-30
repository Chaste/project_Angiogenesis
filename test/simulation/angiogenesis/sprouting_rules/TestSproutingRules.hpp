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
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "Vessel.hpp"
#include "VascularNetwork.hpp"
#include "VascularNode.hpp"
#include "FunctionMap.hpp"
#include "LatticeBasedSproutingRule.hpp"
#include "Owen2011SproutingRule.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestSproutingRules : public AbstractCellBasedWithTimingsTestSuite
{

public:

    void TestLatticeBasedSproutingRuleSimpleNetwork() throw(Exception)
    {
        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        unsigned num_nodes = 100;
        double spacing = 1.0;
        for(unsigned idx=0; idx<num_nodes-1; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*spacing, 0.0, 0.0));
        }

        boost::shared_ptr<Vessel<3> > p_vessel = Vessel<3>::Create(bottom_nodes);
        boost::shared_ptr<VascularNetwork<3> > p_network = VascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel);

        // Set up a sprouting rule
        boost::shared_ptr<LatticeBasedSproutingRule<3> > p_sprouting_rule = LatticeBasedSproutingRule<3>::Create();
        p_sprouting_rule->SetSproutingProbability(0.2);
        p_sprouting_rule->SetVesselNetwork(p_network);

        // Test that we get some, but not all, sprouts
        RandomNumberGenerator::Instance()->Reseed(522525);
        std::vector<boost::shared_ptr<VascularNode<3> > > nodes = p_network->GetNodes();
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        for(unsigned idx=0; idx<5; idx++)
        {
            unsigned num_sprouts = p_sprouting_rule->GetSprouts(nodes).size();
            unsigned num_nodes = bottom_nodes.size();
            TS_ASSERT(num_sprouts>0)
            TS_ASSERT(num_sprouts<num_nodes)
        }
    }

    void TestOwen2011SproutingRuleSimpleNetwork() throw(Exception)
    {
        // Set up the grid
        boost::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
        double spacing = 40.0; //um
        p_grid->SetSpacing(spacing);
        std::vector<unsigned> extents(3, 1);
        extents[0] = 101; // num x
        extents[1] = 11; // num_y
        extents[2] = 1; // num_z
        p_grid->SetExtents(extents);

        // Make a network
        std::vector<boost::shared_ptr<VascularNode<2> > > bottom_nodes;
        unsigned num_nodes = 100;
        for(unsigned idx=0; idx<num_nodes-1; idx++)
        {
            bottom_nodes.push_back(VascularNode<2>::Create(double(idx)*spacing +spacing, 5.0 * spacing));
        }

        boost::shared_ptr<Vessel<2> > p_vessel = Vessel<2>::Create(bottom_nodes);
        boost::shared_ptr<VascularNetwork<2> > p_network = VascularNetwork<2>::Create();
        p_network->AddVessel(p_vessel);
        p_grid->SetVesselNetwork(p_network);

        // Set up a vegf field, 0.15 nM
        boost::shared_ptr<FunctionMap<2> > p_funciton_map = FunctionMap<2>::Create();
        p_funciton_map->SetGrid(p_grid);
        std::vector<double> vegf_field = std::vector<double>(extents[0]*extents[1], 0.15);
        p_funciton_map->UpdateSolution(vegf_field);

        // Set up a sprouting rule
        boost::shared_ptr<Owen2011SproutingRule<2> > p_sprouting_rule = Owen2011SproutingRule<2>::Create();
        p_sprouting_rule->SetSproutingProbability(0.2);
        p_sprouting_rule->SetGrid(p_grid);
        p_sprouting_rule->SetVesselNetwork(p_network);
        p_sprouting_rule->SetHybridSolver(p_funciton_map);

        // Test that we get some, but not all, sprouts
        RandomNumberGenerator::Instance()->Reseed(522525);
        std::vector<boost::shared_ptr<VascularNode<2> > > nodes = p_network->GetNodes();
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        for(unsigned idx=0; idx<5; idx++)
        {
            unsigned num_sprouts = p_sprouting_rule->GetSprouts(nodes).size();
            unsigned num_nodes = bottom_nodes.size();
            TS_ASSERT(num_sprouts>0)
            TS_ASSERT(num_sprouts<num_nodes)
        }
    }
};

#endif
