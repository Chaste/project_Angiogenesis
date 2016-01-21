/*

Copyright (c) 2005-2014, University of Oxford.
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

#ifndef TestLatticeBasedMigrationRules_hpp
#define TestLatticeBasedMigrationRules_hpp

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VasculatureGenerator.hpp"
#include "VascularNode.hpp"
#include "CaVesselSegment.hpp"
#include "CaVessel.hpp"
#include "CaVascularNetwork.hpp"
#include "Part.hpp"
#include "AngiogenesisSolver.hpp"
#include "CaVesselSegment.hpp"
#include "LatticeBasedMigrationRule.hpp"
#include "FunctionMap.hpp"
#include "Owen2011MigrationRule.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestLatticeBasedMigrationRules : public AbstractCellBasedWithTimingsTestSuite
{

public:

    void TestSingleVessel() throw(Exception)
    {
        // Set the grid to move on
        boost::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
        double spacing = 100.0; //um
        p_grid->SetSpacing(spacing);
        std::vector<unsigned> extents(3, 1);
        extents[0] = 7; // num x
        extents[1] = 5; // num_y
        extents[2] = 1; // num_z
        p_grid->SetExtents(extents);

        // Make a vessel
        boost::shared_ptr<VascularNode<2> > p_node1 = VascularNode<2>::Create(0.0, 2.0*spacing);
        boost::shared_ptr<VascularNode<2> > p_node2 = VascularNode<2>::Create(spacing, 2.0*spacing);
        p_node2->SetIsMigrating(true);
        boost::shared_ptr<CaVessel<2> > p_vessel = CaVessel<2>::Create(p_node1, p_node2);
        boost::shared_ptr<CaVascularNetwork<2> > p_network = CaVascularNetwork<2>::Create();
        p_network->AddVessel(p_vessel);
        p_grid->SetVesselNetwork(p_network);

        // Set up the migration rule
        boost::shared_ptr<LatticeBasedMigrationRule<2> > p_migration_rule = LatticeBasedMigrationRule<2>::Create();
        p_migration_rule->SetGrid(p_grid);
        p_migration_rule->SetMovementProbability(0.1);
        p_migration_rule->SetNetwork(p_network);

        // Test that we move into the correct locations and that sometimes, but not always, we don't move
        RandomNumberGenerator::Instance()->Reseed(522525);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        unsigned not_moved = 0;
        for(unsigned idx=0; idx<100; idx++)
        {
            std::vector<int> indices = p_migration_rule->GetIndices(std::vector<boost::shared_ptr<VascularNode<2> > > (1, p_node2));
            if (indices[0] == -1)
            {
                not_moved++;
            }
            else
            {
                TS_ASSERT(indices[0] == 16 or indices[0] == 22 or indices[0] == 8)
                std::cout << indices[0];
            }
        }
        TS_ASSERT(not_moved>0)
        TS_ASSERT(not_moved<100)
    }

    void TestOwen2011SingleVessel() throw(Exception)
    {
        // Set the grid to move on
        boost::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
        double spacing = 100.0; //um
        p_grid->SetSpacing(spacing);
        std::vector<unsigned> extents(3, 1);
        extents[0] = 7; // num x
        extents[1] = 5; // num_y
        extents[2] = 1; // num_z
        p_grid->SetExtents(extents);

        // Make a vessel
        boost::shared_ptr<VascularNode<2> > p_node1 = VascularNode<2>::Create(0.0, 2.0*spacing);
        boost::shared_ptr<VascularNode<2> > p_node2 = VascularNode<2>::Create(spacing, 2.0*spacing);
        p_node2->SetIsMigrating(true);
        boost::shared_ptr<CaVessel<2> > p_vessel = CaVessel<2>::Create(p_node1, p_node2);
        boost::shared_ptr<CaVascularNetwork<2> > p_network = CaVascularNetwork<2>::Create();
        p_network->AddVessel(p_vessel);
        p_grid->SetVesselNetwork(p_network);

        // Set up a vegf field
        boost::shared_ptr<FunctionMap<2> > p_funciton_map = FunctionMap<2>::Create();
        p_funciton_map->SetGrid(p_grid);
        std::vector<double> vegf_field = std::vector<double>(extents[0]*extents[1], 0.0);
        double max_vegf = 0.2; //nM
        for(unsigned idx=0; idx<p_grid->GetNumberOfPoints(); idx++)
        {
            vegf_field[idx] = max_vegf * p_grid->GetLocationOf1dIndex(idx)[0] / (float(extents[0]) * spacing);
        }
        p_funciton_map->SetPointSolution(vegf_field);

        // Set up the migration rule
        boost::shared_ptr<Owen2011MigrationRule<2> > p_migration_rule = Owen2011MigrationRule<2>::Create();
        p_migration_rule->SetGrid(p_grid);
        p_migration_rule->SetMovementProbability(0.1);
        p_migration_rule->SetNetwork(p_network);
        p_migration_rule->SetHybridSolver(p_funciton_map);

        // Test that we move into the correct locations and that sometimes, but not always, we don't move.
        // Also check that we mostly move in the direction of the vegf gradient, but not always
        RandomNumberGenerator::Instance()->Reseed(522525);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        unsigned not_moved = 0;
        unsigned num_right = 0;
        for(unsigned idx=0; idx<100; idx++)
        {
            std::vector<int> indices = p_migration_rule->GetIndices(std::vector<boost::shared_ptr<VascularNode<2> > > (1, p_node2));
            if (indices[0] == -1)
            {
                not_moved++;
            }
            else
            {
                TS_ASSERT(indices[0] == 16 or indices[0] == 22 or indices[0] == 8)
                if(indices[0] == 16)
                {
                    num_right ++;
                }
            }
        }
        TS_ASSERT(not_moved>0)
        TS_ASSERT(not_moved<100)
        TS_ASSERT(num_right>0)
        TS_ASSERT(num_right<100)
    }

//    void DontTestSproutingWithFlow() throw(Exception)
//    {
//        // Make a network
//        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
//        for(unsigned idx=0; idx<9; idx++)
//        {
//            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 10.0, 0.0));
//        }
//        bottom_nodes[0]->GetFlowProperties()->SetIsInputNode(true);
//        bottom_nodes[0]->GetFlowProperties()->SetPressure(3000);
//        bottom_nodes[8]->GetFlowProperties()->SetIsOutputNode(true);
//        bottom_nodes[8]->GetFlowProperties()->SetPressure(1000);
//
//        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(bottom_nodes);
//        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
//        p_network->AddVessel(p_vessel1);
//        p_network->SetSegmentRadii(10.0);
//
//        for(unsigned idx=1; idx<6; idx+=2)
//        {
//            p_network->FormSprout(ChastePoint<3>(double(idx)*10, 10.0, 0.0), ChastePoint<3>(double(idx)*10, 20.0, 0.0));
//        }
//
//        p_network->UpdateSegments();
//        std::vector<boost::shared_ptr<CaVesselSegment<3> > > segments = p_network->GetVesselSegments();
//        for(unsigned idx=0; idx<segments.size(); idx++)
//        {
//            segments[idx]->GetFlowProperties()->SetViscosity(1.e-3);
//        }
//
//        boost::shared_ptr<OnLatticeRwGrowthDirectionModifier<3> > p_grow_direction_modifier =
//                boost::shared_ptr<OnLatticeRwGrowthDirectionModifier<3> >(new OnLatticeRwGrowthDirectionModifier<3>());
//
//        // Grow the vessel
//        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
//        AngiogenesisSolver<3> angiogenesis_solver;
//        angiogenesis_solver.SetVesselNetwork(p_network);
//        angiogenesis_solver.SetOutputDirectory("TestOnLatticeRwGrowthDirectionModifier/SproutingFlow/");
//        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
//        angiogenesis_solver.Run();
//    }
};

#endif // TestLatticeBasedMigrationRules_hpp
