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
#include "LatticeBasedSproutingRule.hpp"
#include "Owen2011LatticeBasedSproutingRule.hpp"
#include "CaVessel.hpp"
#include "CaVascularNetwork.hpp"
#include "VascularNode.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FunctionMap.hpp"
#include "Debug.hpp"

class TestSproutingRules : public AbstractCellBasedTestSuite
{

public:

    void TestLatticeBasedSproutingRuleSimpleNetwork() throw(Exception)
    {
        // Set up the grid
        boost::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
        double spacing = 40.0; //um
        p_grid->SetSpacing(spacing);
        std::vector<unsigned> extents(3, 1);
        extents[0] = 11; // num x
        extents[1] = 11; // num_y
        extents[2] = 11; // num_z
        p_grid->SetExtents(extents);

        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        unsigned num_nodes = 10;
        for(unsigned idx=0; idx<num_nodes-1; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*spacing +spacing, 5.0 * spacing, 5.0 * spacing));
        }

        boost::shared_ptr<CaVessel<3> > p_vessel = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel);
        p_grid->SetVesselNetwork(p_network);

        // Set up a sprouting rule
        boost::shared_ptr<LatticeBasedSproutingRule<3> > p_sprouting_rule = LatticeBasedSproutingRule<3>::Create();
        p_sprouting_rule->SetSproutingProbability(0.1);
        p_sprouting_rule->SetGrid(p_grid);
        p_sprouting_rule->SetVesselNetwork(p_network);

        // Get the sprout directions
        RandomNumberGenerator::Instance()->Reseed(522525);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(40.0, 400);
        std::vector<boost::shared_ptr<VascularNode<3> > > nodes = p_network->GetNodes();

        // Repeat a few times
        for (unsigned idx=0; idx<1; idx++)
        {
            std::vector<c_vector<double, 3> > directions = p_sprouting_rule->GetSproutDirections(nodes);
            for (unsigned jdx=0; jdx<directions.size(); jdx++)
            {
                std::cout << directions[jdx] << std::endl;
            }
        }
    }

    void TestOwen2011SproutingRuleSimpleNetwork() throw(Exception)
    {
        // Set up the grid
        boost::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
        double spacing = 40.0; //um
        p_grid->SetSpacing(spacing);
        std::vector<unsigned> extents(3, 1);
        extents[0] = 11; // num x
        extents[1] = 11; // num_y
        extents[2] = 11; // num_z
        p_grid->SetExtents(extents);

        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        unsigned num_nodes = 10;
        for(unsigned idx=0; idx<num_nodes-1; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*spacing +spacing, 5.0 * spacing, 5.0 * spacing));
        }

        boost::shared_ptr<CaVessel<3> > p_vessel = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel);
        p_grid->SetVesselNetwork(p_network);

        // Set up a vegf field
        boost::shared_ptr<FunctionMap<3> > p_funciton_map = FunctionMap<3>::Create();
        p_funciton_map->SetGrid(p_grid);
        std::vector<double> vegf_field = std::vector<double>(extents[0]*extents[1]*extents[2], 0.0);
        for(unsigned idx=0; idx<extents[0]*extents[1]*extents[2]; idx++)
        {
            vegf_field[idx] = p_grid->GetLocationOf1dIndex(idx)[0] / (spacing * extents[0]); // 0 to 1 nM across the grid x extents
        }
        p_funciton_map->SetPointSolution(vegf_field);

        // Set up a sprouting rule
        boost::shared_ptr<Owen2011LatticeBasedSproutingRule<3> > p_sprouting_rule = Owen2011LatticeBasedSproutingRule<3>::Create();
        p_sprouting_rule->SetSproutingProbability(0.1);
        p_sprouting_rule->SetGrid(p_grid);
        p_sprouting_rule->SetVesselNetwork(p_network);
        p_sprouting_rule->SetHybridSolver(p_funciton_map);

        // Get the sprout directions
        RandomNumberGenerator::Instance()->Reseed(522525);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(40.0, 400);
        std::vector<boost::shared_ptr<VascularNode<3> > > nodes = p_network->GetNodes();

        // Repeat a few times
        for (unsigned idx=0; idx<1000; idx++)
        {
            std::vector<c_vector<double, 3> > directions = p_sprouting_rule->GetSproutDirections(nodes);
        }
    }
};

#endif
