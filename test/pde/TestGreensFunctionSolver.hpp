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

#ifndef TESTGREENSFUNCTIONSOLVER_HPP_
#define TESTGREENSFUNCTIONSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "GreensFunctionSolver.hpp"
#include "VascularNetwork.hpp"
#include "VasculatureGenerator.hpp"
#include "HybridLinearEllipticPde.hpp"
#include "RegularGrid.hpp"
#include "UnitCollections.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestGreenFunctionSolver : public CxxTest::TestSuite
{

public:

    void TestSingleVessel3d()
    {
        double vessel_length = 2.0;
        VasculatureGenerator<3> generator;
        c_vector<double,3> centre = zero_vector<double>(3);
        centre[0] = 0.5;
        centre[1] = 0.5;
        centre[2] = 0.0;

        boost::shared_ptr<VascularNetwork<3> > p_network = generator.GenerateSingleVessel(vessel_length, centre, 14.0);
        p_network->SetSegmentRadii(0.05*unit::metres);

        // Set up the grid
        boost::shared_ptr<Part<3> > p_domain = Part<3>::Create();
        p_domain->AddCuboid(1.0, 1.0, 2.0);
        boost::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
        p_grid->GenerateFromPart(p_domain, 0.1);

        // Choose the PDE
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_pde = HybridLinearEllipticPde<3>::Create();
        p_pde->SetIsotropicDiffusionConstant(1);
        p_pde->SetContinuumConstantInUTerm(-2.0);

        GreensFunctionSolver<3> solver;
        solver.SetVesselNetwork(p_network);
        solver.SetGrid(p_grid);
        solver.SetPde(p_pde);

        MAKE_PTR_ARGS(OutputFileHandler, p_output_file_handler, ("TestGreensFunctionSolver/KroghCylinder3d", false));
        solver.SetFileHandler(p_output_file_handler);
        solver.Setup();
        solver.SetWriteSolution(true);
        solver.Solve();
    }
};

#endif /*TESTGREENSFUNCTIONSOLVER_HPP_*/
