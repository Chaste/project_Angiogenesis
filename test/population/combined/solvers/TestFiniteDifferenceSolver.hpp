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

#ifndef TESTFINITEDIFFERENCESOLVER_HPP_
#define TESTFINITEDIFFERENCESOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include "UblasIncludes.hpp"
#include "Part.hpp"
#include "HybridLinearEllipticPde.hpp"
#include "FiniteDifferenceSolver.hpp"
#include "CaVascularNetwork.hpp"
#include "VasculatureGenerator.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "OutputFileHandler.hpp"
#include "VasculatureGenerator.hpp"
//#include "SimpleCell.hpp"
//#include "SimpleCellPopulation.hpp"

class TestFiniteDifferenceSolver : public CxxTest::TestSuite
{

public:

    void Test3dKroghCylinderNetwork()
    {
        // Set up the vessel network
        double vessel_length = 100;
        VasculatureGenerator<3> generator;
        boost::shared_ptr<CaVascularNetwork<3> > p_network = generator.GenerateSingleVessel(vessel_length);

        // Set up the PDE domain
        boost::shared_ptr<Part<3> > p_domain = Part<3>::Create();
        p_domain->AddCuboid(vessel_length, vessel_length, vessel_length);

        // Choose the PDE
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_pde = HybridLinearEllipticPde<3>::Create();
        p_pde->SetDiffusionConstant(0.0033);
        p_pde->SetLinearInUTerm(-2.e-7);

        // Set up and run the solver
        FiniteDifferenceSolver<3> solver;
        solver.SetVesselNetwork(p_network);
        solver.SetExtents(p_domain, 10.0);
        solver.SetPde(p_pde);

        OutputFileHandler output_file_handler("TestFiniteDifferenceSolver/KroghCylinder3d", false);
        solver.SetWorkingDirectory(output_file_handler.GetOutputDirectoryFullPath());
        solver.Solve(true);
    }

//    void Test3dKroghCylinderNetworkWithCells()
//    {
//        // Set up the vessel network
//        double vessel_length = 100;
//        VasculatureGenerator<3> generator;
//        boost::shared_ptr<CaVascularNetwork<3> > p_network = generator.GenerateSingleVessel(vessel_length);
//
//        // Set up the PDE domain
//        boost::shared_ptr<Part<3> > p_domain = Part<3>::Create();
//        p_domain->AddCuboid(vessel_length, vessel_length, vessel_length);
//
//        // Set up the cells
//        boost::shared_ptr<SimpleCellPopulation<3> > p_population = SimpleCellPopulation<3>::Create();
//        double spacing = 10;
//        unsigned num_x = unsigned(vessel_length/spacing) + 1;
//        p_population->GenerateCellsOnGrid(num_x, num_x, num_x, spacing);
//        p_population->BooleanWithVesselNetwork(p_network);
//
//        // Choose the PDE
//        boost::shared_ptr<HybridLinearEllipticPde<3> > p_pde = HybridLinearEllipticPde<3>::Create();
//        p_pde->SetDiffusionConstant(0.0033);
//        p_pde->SetLinearInUTerm(-2.e-7);
//
//        // Set up and run the solver
//        FiniteDifferenceSolver<3> solver;
//        solver.SetVesselNetwork(p_network);
//        solver.SetCellPopulation(p_population);
//        solver.SetExtents(p_domain, 10.0);
//        solver.SetPde(p_pde);
//
//        OutputFileHandler output_file_handler("TestFiniteDifferenceSolver/KroghCylinder3dCells", false);
//        solver.SetWorkingDirectory(output_file_handler.GetOutputDirectoryFullPath());
//        solver.Solve(true);
//    }

    void Test2dBifurcationNetwork()
    {
        // Set up the vessel network
        double vessel_length = 100;
        VasculatureGenerator<3> generator;
        boost::shared_ptr<CaVascularNetwork<3> > p_network = generator.GenerateBifurcationUnit(vessel_length);

        // Set up the tissue domain
        boost::shared_ptr<Part<3> > p_domain = Part<3>::Create();
        p_domain->AddRectangle(4.0 * vessel_length, 2.0 * vessel_length);

        // Define the PDE
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_pde = HybridLinearEllipticPde<3>::Create();
        p_pde->SetDiffusionConstant(0.0033);
        p_pde->SetLinearInUTerm(-2.e-7);

        // Set up and run the simulation
        FiniteDifferenceSolver<3> solver;
        solver.SetVesselNetwork(p_network);
        solver.SetExtents(p_domain, 10.0);
        solver.SetPde(p_pde);

        OutputFileHandler output_file_handler("TestFiniteDifferenceSolver/Bifurcation2d", false);
        solver.SetWorkingDirectory(output_file_handler.GetOutputDirectoryFullPath());
        solver.Solve(true);
    }

    void Test3dBifurcationNetwork()
    {
        // Set up the vessel network
        double vessel_length = 100;
        VasculatureGenerator<3> generator;
        boost::shared_ptr<CaVascularNetwork<3> > p_network = generator.GenerateBifurcationUnit(vessel_length);

        // Set up the tissue domain
        boost::shared_ptr<Part<3> > p_domain = Part<3>::Create();
        p_domain->AddCuboid(4.0 * vessel_length, 2.0 * vessel_length, 2.0 * vessel_length);

        // Define the PDE
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_pde = HybridLinearEllipticPde<3>::Create();
        p_pde->SetDiffusionConstant(0.0033);
        p_pde->SetLinearInUTerm(-2.e-7);

        // Set up and run the simulation
        FiniteDifferenceSolver<3> solver;
        solver.SetVesselNetwork(p_network);
        solver.SetExtents(p_domain, 10.0);
        solver.SetPde(p_pde);

        OutputFileHandler output_file_handler("TestFiniteDifferenceSolver/Bifurcation3d", false);
        solver.SetWorkingDirectory(output_file_handler.GetOutputDirectoryFullPath());
        solver.Solve(true);
    }
};

#endif /*TESTFINITEDIFFERENCESOLVER_HPP_*/
