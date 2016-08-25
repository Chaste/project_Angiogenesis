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
#include "VesselNode.hpp"
#include "VesselSegment.hpp"
#include "Vessel.hpp"
#include "VesselNetwork.hpp"
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
        p_grid->SetSpacing(spacing*1.e-6*unit::metres);

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
        boost::shared_ptr<VesselNetwork<3> > p_network = generator.GenerateSingleVessel(length, start_point,
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
