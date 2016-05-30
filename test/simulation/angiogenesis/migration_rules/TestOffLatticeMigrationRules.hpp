//
//  TestSimpleStructuralAdaptationSolver.hpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 03/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#ifndef TestOffLatticePrwGrowthDirectionModifier_hpp
#define TestOffLatticePrwGrowthDirectionModifier_hpp

#include <cxxtest/TestSuite.h>
#include "OffLatticeMigrationRule.hpp"
#include "OffLatticeSproutingRule.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "FunctionMap.hpp"
#include "VasculatureGenerator.hpp"
#include "VascularNode.hpp"
#include "VesselSegment.hpp"
#include "Vessel.hpp"
#include "VascularNetwork.hpp"
#include "Part.hpp"
#include "FiniteDifferenceSolver.hpp"
#include "AngiogenesisSolver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "VesselSegment.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FlowSolver.hpp"

class TestOffLatticeMigrationRules : public AbstractCellBasedTestSuite
{

public:

    void Test2dMigration() throw(Exception)
    {
        MAKE_PTR_ARGS(OutputFileHandler, p_handler, ("TestOffLatticeMigrationRules2d"));

        // Set up the grid
        boost::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
        double spacing = 40.0; //um
        p_grid->SetSpacing(spacing);

        std::vector<unsigned> extents(3, 1);
        extents[0] = 25; // num x
        extents[1] = 25; // num_y
        extents[2] = 25; // num_z
        p_grid->SetExtents(extents);

        // Prescribe a linearly increasing vegf field using a function map
        boost::shared_ptr<FunctionMap<3> > p_funciton_map = FunctionMap<3>::Create();
        p_funciton_map->SetGrid(p_grid);
        std::vector<double> vegf_field = std::vector<double>(extents[0] * extents[1] * extents[2], 0.0);
        for (unsigned idx = 0; idx < extents[0] * extents[1] * extents[2]; idx++)
        {
            vegf_field[idx] = 0.2*p_grid->GetLocationOf1dIndex(idx)[0] / (spacing * extents[0]);
        }

        p_grid->Write(p_handler);
        p_funciton_map->SetFileHandler(p_handler);
        p_funciton_map->SetFileName("Function.vti");
        p_funciton_map->Setup();
        p_funciton_map->UpdateSolution(vegf_field);
        p_funciton_map->Write();

        //Set up the limbal vessel
        VasculatureGenerator<3> generator;
        c_vector<double, 3> start_point;
        start_point[0] = 2.0 * spacing; // two lattice points to the right
        start_point[1] = 2.0 * spacing;
        start_point[2] = 10.0*spacing;

        double length = spacing * (extents[1] - 3); // full domain in y direction
        unsigned divisions = extents[1] - 2; // divide the vessel to coincide with grid
        unsigned alignment_axis = 1; // pointing y direction
        boost::shared_ptr<VascularNetwork<3> > p_network = generator.GenerateSingleVessel(length, start_point,
                                                                                            divisions, alignment_axis);

        boost::shared_ptr<OffLatticeMigrationRule<3> > p_migration_rule = OffLatticeMigrationRule<3>::Create();
        p_migration_rule->SetHybridSolver(p_funciton_map); // This contains the vegf field
        p_migration_rule->SetNetwork(p_network);

        boost::shared_ptr<OffLatticeSproutingRule<3> > p_sprouting_rule = OffLatticeSproutingRule<3>::Create();
        p_sprouting_rule->SetHybridSolver(p_funciton_map); // This contains the vegf field
        p_sprouting_rule->SetSproutingProbability(0.01);
        p_sprouting_rule->SetVesselNetwork(p_network);

        AngiogenesisSolver<3> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetMigrationRule(p_migration_rule);
        angiogenesis_solver.SetSproutingRule(p_sprouting_rule);
        angiogenesis_solver.SetOutputFileHandler(p_handler);

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(100.0, 100);
        angiogenesis_solver.Run(true);
    }

};

#endif
