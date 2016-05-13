/*

 Copyright (c) 2005-2015, University of Oxford.
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

#ifndef TESTSNAILTRAIL_HPP_
#define TESTSNAILTRAIL_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>

#include "Owen2011SproutingRule.hpp"
#include "Owen2011MigrationRule.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "RegularGrid.hpp"
#include "VasculatureGenerator.hpp"
#include "CaVascularNetwork.hpp"
#include "AngiogenesisSolver.hpp"
#include "FunctionMap.hpp"
#include "FlowSolver.hpp"
#include "VascularTumourSolver.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestSnailTrail : public AbstractCellBasedWithTimingsTestSuite
{

public:

    void DontTest2dSnailTrailWithPrescribedVegf() throw (Exception)
    {
        MAKE_PTR_ARGS(OutputFileHandler, p_handler, ("TestSnailTrailModel"));

        // Set up the grid
        boost::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
        double spacing = 40.0; //um
        p_grid->SetSpacing(spacing);

        std::vector<unsigned> extents(3, 1);
        extents[0] = 25; // num x
        extents[1] = 25; // num_y
        p_grid->SetExtents(extents);

        // Prescribe a linearly increasing vegf field using a function map
        std::vector<double> vegf_field = std::vector<double>(extents[0] * extents[1], 0.0);
        for (unsigned idx = 0; idx < extents[0] * extents[1]; idx++)
        {
            vegf_field[idx] = 0.2*p_grid->GetLocationOf1dIndex(idx)[0] / (spacing * extents[0]);
        }

        boost::shared_ptr<FunctionMap<3> > p_funciton_map = FunctionMap<3>::Create();
        p_funciton_map->SetGrid(p_grid);
        p_funciton_map->SetFileHandler(p_handler);
        p_funciton_map->SetFileName("Function.vti");
        p_funciton_map->Setup();
        p_funciton_map->UpdateSolution(vegf_field);
        p_funciton_map->Write();

        //Set up the limbal vessel
        VasculatureGenerator<3> generator;
        c_vector<double, 3> start_point;
        start_point[0] = 2.0 * spacing; // two lattice points to the right
        start_point[1] = 0.0;
        start_point[2] = 0.0;

        double length = spacing * (extents[1] - 1); // full domain in y direction
        unsigned divisions = extents[1] - 2; // divide the vessel to coincide with grid
        unsigned alignment_axis = 1; // pointing y direction
        boost::shared_ptr<CaVascularNetwork<3> > p_network = generator.GenerateSingleVessel(length, start_point,
                                                                                            divisions, alignment_axis);

        boost::shared_ptr<Owen2011MigrationRule<3> > p_migration_rule = Owen2011MigrationRule<3>::Create();
        p_migration_rule->SetGrid(p_grid);
        p_migration_rule->SetHybridSolver(p_funciton_map); // This contains the vegf field
        p_migration_rule->SetCellMotilityParameter(100.0);
        p_migration_rule->SetCellChemotacticParameter(80000.0);
        p_migration_rule->SetNetwork(p_network);

        boost::shared_ptr<Owen2011SproutingRule<3> > p_sprouting_rule = Owen2011SproutingRule<3>::Create();
        p_sprouting_rule->SetHybridSolver(p_funciton_map); // This contains the vegf field
        p_sprouting_rule->SetSproutingProbability(0.5);
        p_sprouting_rule->SetGrid(p_grid);
        p_sprouting_rule->SetVesselNetwork(p_network);

        AngiogenesisSolver<3> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetVesselGrid(p_grid);
        angiogenesis_solver.SetMigrationRule(p_migration_rule);
        angiogenesis_solver.SetSproutingRule(p_sprouting_rule);
        angiogenesis_solver.SetOutputFileHandler(p_handler);

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(100.0, 100);
        angiogenesis_solver.Run(true);
    }

    void DontTest3dSnailTrailWithPrescribedVegf() throw (Exception)
    {
        MAKE_PTR_ARGS(OutputFileHandler, p_handler, ("TestSnailTrailModel3d"));

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
        std::vector<double> vegf_field = std::vector<double>(extents[0] * extents[1]* extents[2], 0.0);
        for (unsigned idx = 0; idx < extents[0] * extents[1]* extents[2]; idx++)
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
        start_point[1] = 0.0;
        start_point[2] = 10.0 * spacing; // roughly midway

        double length = spacing * (extents[1] - 1); // full domain in y direction
        unsigned divisions = extents[1] - 2; // divide the vessel to coincide with grid
        unsigned alignment_axis = 1; // pointing y direction
        boost::shared_ptr<CaVascularNetwork<3> > p_network = generator.GenerateSingleVessel(length, start_point,
                                                                                            divisions, alignment_axis);

        boost::shared_ptr<Owen2011MigrationRule<3> > p_migration_rule = Owen2011MigrationRule<3>::Create();
        p_migration_rule->SetGrid(p_grid);
        p_migration_rule->SetHybridSolver(p_funciton_map); // This contains the vegf field
        p_migration_rule->SetCellMotilityParameter(100.0);
        p_migration_rule->SetCellChemotacticParameter(80000.0);
        p_migration_rule->SetNetwork(p_network);

        boost::shared_ptr<Owen2011SproutingRule<3> > p_sprouting_rule = Owen2011SproutingRule<3>::Create();
        p_sprouting_rule->SetHybridSolver(p_funciton_map); // This contains the vegf field
        p_sprouting_rule->SetSproutingProbability(0.05); // reduce for 3d
        p_sprouting_rule->SetGrid(p_grid);
        p_sprouting_rule->SetVesselNetwork(p_network);

        AngiogenesisSolver<3> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetVesselGrid(p_grid);
        angiogenesis_solver.SetMigrationRule(p_migration_rule);
        angiogenesis_solver.SetSproutingRule(p_sprouting_rule);
        angiogenesis_solver.SetOutputFileHandler(p_handler);

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(100.0, 100);
        angiogenesis_solver.Run(true);
    }

    void Test2dSnailTrailWithPrescribedVegfAndFlow() throw (Exception)
    {
        MAKE_PTR_ARGS(OutputFileHandler, p_handler, ("TestSnailTrailModel/Flow"));

        // Set up the grid
        boost::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
        double spacing = 40.0; //um
        p_grid->SetSpacing(spacing);

        std::vector<unsigned> extents(3, 1);
        extents[0] = 25; // num x
        extents[1] = 25; // num_y
        p_grid->SetExtents(extents);

        // Prescribe a linearly increasing vegf field using a function map
        boost::shared_ptr<FunctionMap<3> > p_funciton_map = FunctionMap<3>::Create();
        p_funciton_map->SetGrid(p_grid);
        std::vector<double> vegf_field = std::vector<double>(extents[0] * extents[1], 0.0);
        for (unsigned idx = 0; idx < extents[0] * extents[1]; idx++)
        {
            vegf_field[idx] = 0.2 + p_grid->GetLocationOf1dIndex(idx)[0] / (spacing * extents[0]); // 0.1 to 1.1 nM across the grid x extents
        }

        p_funciton_map->SetFileHandler(p_handler);
        p_funciton_map->SetFileName("Function.vti");
        p_funciton_map->Setup();
        p_funciton_map->UpdateSolution(vegf_field);
        p_funciton_map->Write();

        //Set up the limbal vessel
        VasculatureGenerator<3> generator;
        c_vector<double, 3> start_point;
        start_point[0] = 2.0 * spacing; // two lattice points to the right
        start_point[1] = 0.0;
        start_point[2] = 0.0;

        double length = spacing * (extents[1] - 1); // full domain in y direction
        unsigned divisions = extents[1] - 2; // divide the vessel to coincide with grid
        unsigned alignment_axis = 1; // pointing y direction
        boost::shared_ptr<CaVascularNetwork<3> > p_network = generator.GenerateSingleVessel(length, start_point,
                                                                                            divisions, alignment_axis);

        p_network->SetNodeRadii(10.0);
        p_network->SetSegmentRadii(10.0);
        p_network->GetVesselSegments()[0]->GetFlowProperties()->SetViscosity(1.e-9);
        p_network->GetVessels()[0]->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
        p_network->GetVessels()[0]->GetStartNode()->GetFlowProperties()->SetPressure(0.005);
        p_network->GetVessels()[0]->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
        p_network->GetVessels()[0]->GetEndNode()->GetFlowProperties()->SetPressure(0.003);

        boost::shared_ptr<Owen2011MigrationRule<3> > p_migration_rule = Owen2011MigrationRule<3>::Create();
        p_migration_rule->SetGrid(p_grid);
        p_migration_rule->SetCellMotilityParameter(0.001);
        p_migration_rule->SetNetwork(p_network);
        p_migration_rule->SetHybridSolver(p_funciton_map); // This contains the vegf field

        boost::shared_ptr<Owen2011SproutingRule<3> > p_sprouting_rule = Owen2011SproutingRule<3>::Create();
        p_sprouting_rule->SetHybridSolver(p_funciton_map); // This contains the vegf field
        p_sprouting_rule->SetSproutingProbability(0.3);
        p_sprouting_rule->SetGrid(p_grid);
        p_sprouting_rule->SetVesselNetwork(p_network);

        boost::shared_ptr<AngiogenesisSolver<3> > p_angiogenesis_solver = AngiogenesisSolver<3>::Create();
        p_angiogenesis_solver->SetVesselNetwork(p_network);
        p_angiogenesis_solver->SetVesselGrid(p_grid);
        p_angiogenesis_solver->SetMigrationRule(p_migration_rule);
        p_angiogenesis_solver->SetSproutingRule(p_sprouting_rule);
        p_angiogenesis_solver->SetOutputFileHandler(p_handler);
        p_angiogenesis_solver->SetVesselGrid(p_grid);

        boost::shared_ptr<FlowSolver<3> > p_flow_solver = FlowSolver<3>::Create();
        std::vector<boost::shared_ptr<CaVesselSegment<3> > > segments = p_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(1.e-3);
        }

        boost::shared_ptr<VascularTumourSolver<3> > p_vascular_tumour_solver = VascularTumourSolver<3>::Create();
        p_vascular_tumour_solver->SetVesselNetwork(p_network);
        p_vascular_tumour_solver->SetOutputFrequency(1);
        p_vascular_tumour_solver->SetAngiogenesisSolver(p_angiogenesis_solver);
        p_vascular_tumour_solver->SetFlowSolver(p_flow_solver);
        p_vascular_tumour_solver->SetOutputFileHandler(p_handler);

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(100.0, 400);
        p_vascular_tumour_solver->Run();
    }
};

#endif /*TESTSNAILTRAIL_HPP_*/
