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

#include "LatticeBasedSproutingRule.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "RegularGrid.hpp"
#include "VasculatureGenerator.hpp"
#include "CaVascularNetwork.hpp"
#include "AngiogenesisSolver.hpp"
#include "FunctionMap.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestSnailTrail : public AbstractCellBasedWithTimingsTestSuite
{

public:

    void Test2dSnailTrailWithPrescribedVegf() throw (Exception)
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
//        boost::shared_ptr<FunctionMap<3> > p_funciton_map = FunctionMap<3>::Create();
//        p_funciton_map->SetGrid(p_grid);
//        std::vector<double> vegf_field = std::vector<double>(extents[0]*extents[1], 0.0);
//        for(unsigned idx=0; idx<extents[0]*extents[1]; idx++)
//        {
//            vegf_field[idx] = p_grid->GetLocationOf1dIndex(idx)[0] / (spacing * extents[0]); // 0 to 1 nM across the grid x extents
//        }
//        p_funciton_map->SetPointSolution(vegf_field);
//
//        std::map<std::string, std::vector<double> > data;
//        data["Vegf"] = vegf_field;
//        p_grid->Write(p_handler);
//        p_funciton_map->SetFileHandler(p_handler);
//        p_funciton_map->SetFileName("Function.vti");
//        p_funciton_map->Setup();
//        p_funciton_map->UpdateSolution(data);
//        p_funciton_map->Write();

        //Set up the limbal vessel
        VasculatureGenerator<3> generator;
        c_vector<double, 3> start_point;
        start_point[0] = 2.0 * spacing; // two lattice points to the right
        start_point[1] = 0.0;
        start_point[2] = 0.0;

        double length = spacing*(extents[1]-1); // full domain in y direction
        unsigned divisions = extents[1] - 2; // divide the vessel to coincide with grid
        unsigned alignment_axis = 1; // pointing y direction
        boost::shared_ptr<CaVascularNetwork<3> > p_network = generator.GenerateSingleVessel(length, start_point, divisions, alignment_axis);

        // Set up the angiogenesis solver
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(40.0, 800);

        AngiogenesisSolver<3> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetVesselGrid(p_grid);

//        boost::shared_ptr<OnLatticeRwGrowthDirectionModifier<3> > p_random_walk_modifier =
//                boost::shared_ptr<OnLatticeRwGrowthDirectionModifier<3> >(new OnLatticeRwGrowthDirectionModifier<3>());
//
//        boost::shared_ptr<OffLatticeSolutionDependentGrowthDirectionModifier<3> > p_solution_dependent_modifier =
//                boost::shared_ptr<OffLatticeSolutionDependentGrowthDirectionModifier<3> >(new OffLatticeSolutionDependentGrowthDirectionModifier<3>());
//
//        boost::shared_ptr<Owen2011LatticeBasedSproutingRule<3> > p_sprouting_rule = Owen2011LatticeBasedSproutingRule<3>::Create();
//        p_sprouting_rule->SetHybridSolver(p_funciton_map); // This contains the vegf field

        boost::shared_ptr<LatticeBasedSproutingRule<3> > p_sprouting_rule = LatticeBasedSproutingRule<3>::Create();

//        angiogenesis_solver.AddGrowthDirectionModifier(p_random_walk_modifier);
//        angiogenesis_solver.AddGrowthDirectionModifier(p_solution_dependent_modifier);
        angiogenesis_solver.SetSproutingRule(p_sprouting_rule);
        angiogenesis_solver.SetOutputFileHandler(p_handler);
        angiogenesis_solver.Run(true);
    }
};

#endif /*TESTSNAILTRAIL_HPP_*/
