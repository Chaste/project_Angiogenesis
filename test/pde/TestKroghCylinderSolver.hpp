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
 * Redistributions in binary form must reproduce the abovea copyright notice,
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

#ifndef TESTDISCRETESOURCE_HPP_
#define TESTDISCRETESOURCE_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include <string>
#include "SmartPointers.hpp"
#include "Part.hpp"
#include "HybridLinearEllipticPde.hpp"
#include "HybridNonLinearEllipticPde.hpp"
#include "FiniteDifferenceSolver.hpp"
#include "VascularNetwork.hpp"
#include "VasculatureGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "RegularGrid.hpp"
#include "DiscreteSource.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestDiscreteSource : public CxxTest::TestSuite
{

public:

    void TestWithVessels() throw(Exception)
    {
        // Set up the vessel network
        double vessel_length = 100;
        VasculatureGenerator<3> generator;
        boost::shared_ptr<VascularNetwork<3> > p_network = generator.GenerateSingleVessel(vessel_length);

        // Set up the grid
        boost::shared_ptr<Part<3> > p_domain = Part<3>::Create();
        p_domain->AddCuboid(vessel_length, vessel_length, vessel_length);
        c_vector<double, 3> translation_vector;
        translation_vector[0] = -vessel_length/2.0;
        translation_vector[1] = -vessel_length/2.0;
        translation_vector[2] = 0.0;

        p_domain->Translate(translation_vector);
        boost::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
        p_grid->GenerateFromPart(p_domain, 10.0);

        // Choose the PDE
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_pde = HybridLinearEllipticPde<3>::Create();
        p_pde->SetIsotropicDiffusionConstant(0.0033);
        p_pde->SetContinuumConstantInUTerm(-2.e-6);

        // Set up the discrete source
        boost::shared_ptr<DiscreteSource<3> > p_vessel_source_lin = DiscreteSource<3>::Create();
        p_vessel_source_lin->SetValue(-1.e3);
        p_vessel_source_lin->SetType(SourceType::VESSEL);
        p_vessel_source_lin->SetSource(SourceStrength::PRESCRIBED);
        p_vessel_source_lin->SetIsLinearInSolution(true);

        boost::shared_ptr<DiscreteSource<3> > p_vessel_source_const = DiscreteSource<3>::Create();
        p_vessel_source_const->SetValue(40.e3);
        p_vessel_source_const->SetType(SourceType::VESSEL);
        p_vessel_source_const->SetSource(SourceStrength::PRESCRIBED);
        p_vessel_source_const->SetIsLinearInSolution(false);

        p_pde->AddDiscreteSource(p_vessel_source_lin);
        p_pde->AddDiscreteSource(p_vessel_source_const);

        // Set up and run the solver
        FiniteDifferenceSolver<3> solver;
        solver.SetGrid(p_grid);
        solver.SetPde(p_pde);
        solver.SetVesselNetwork(p_network);

        MAKE_PTR_ARGS(OutputFileHandler, p_output_file_handler, ("TestDiscreteSource/WithVessels", false));
        solver.SetFileHandler(p_output_file_handler);
        solver.SetWriteSolution(true);
        solver.Solve();
    }

    void TestNonLinearWithVessels() throw(Exception)
    {
        // Set up the vessel network
        double vessel_length = 100;
        VasculatureGenerator<3> generator;
        boost::shared_ptr<VascularNetwork<3> > p_network = generator.GenerateSingleVessel(vessel_length);

        // Set up the grid
        boost::shared_ptr<Part<3> > p_domain = Part<3>::Create();
        p_domain->AddCuboid(vessel_length, vessel_length, vessel_length);
        c_vector<double, 3> translation_vector;
        translation_vector[0] = -vessel_length/2.0;
        translation_vector[1] = -vessel_length/2.0;
        translation_vector[2] = 0.0;

        p_domain->Translate(translation_vector);
        boost::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
        p_grid->GenerateFromPart(p_domain, 10.0);

        // Choose the PDE
        boost::shared_ptr<HybridNonLinearEllipticPde<3> > p_pde = HybridNonLinearEllipticPde<3>::Create();
        p_pde->SetIsotropicDiffusionConstant(0.0033);
        p_pde->SetContinuumConstantInUTerm(-2.e-6);
        p_pde->SetThreshold(2.5);

        // Set up the discrete source
        boost::shared_ptr<DiscreteSource<3> > p_vessel_source_lin = DiscreteSource<3>::Create();
        p_vessel_source_lin->SetValue(-1.e3);
        p_vessel_source_lin->SetType(SourceType::VESSEL);
        p_vessel_source_lin->SetSource(SourceStrength::PRESCRIBED);
        p_vessel_source_lin->SetIsLinearInSolution(true);

        boost::shared_ptr<DiscreteSource<3> > p_vessel_source_const = DiscreteSource<3>::Create();
        p_vessel_source_const->SetValue(40.e3);
        p_vessel_source_const->SetType(SourceType::VESSEL);
        p_vessel_source_const->SetSource(SourceStrength::PRESCRIBED);
        p_vessel_source_const->SetIsLinearInSolution(false);

        p_pde->AddDiscreteSource(p_vessel_source_lin);
        p_pde->AddDiscreteSource(p_vessel_source_const);

        // Set up and run the solver
        FiniteDifferenceSolver<3> solver;
        solver.SetGrid(p_grid);
        solver.SetNonLinearPde(p_pde);
        solver.SetVesselNetwork(p_network);

        MAKE_PTR_ARGS(OutputFileHandler, p_output_file_handler, ("TestDiscreteSource/NonLinearWithVessels", false));
        solver.SetFileHandler(p_output_file_handler);
        solver.SetWriteSolution(true);
        solver.Solve();
    }
};

#endif /*TESTDISCRETESOURCE_HPP_*/
