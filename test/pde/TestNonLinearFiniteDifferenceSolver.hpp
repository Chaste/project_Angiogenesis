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

#ifndef TESTFINITEELEMENTSOLVER_HPP_
#define TESTFINITEELEMENTSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include <string>
#include <boost/lexical_cast.hpp>
#include "FiniteDifferenceSolver.hpp"
#include "UblasIncludes.hpp"
#include "Part.hpp"
#include "Vertex.hpp"
#include "HybridLinearEllipticPde.hpp"
#include "HybridNonLinearEllipticPde.hpp"
#include "VesselNetwork.hpp"
#include "VasculatureGenerator.hpp"
#include "SmartPointers.hpp"
#include "HybridBoundaryCondition.hpp"
#include "HybridMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "Debug.hpp"

class TestNonLinearFiniteDifferenceSolver : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestBox() throw(Exception)
    {
        // Set up the mesh
        boost::shared_ptr<Part<3> > p_domain = Part<3>::Create();
        p_domain->AddCuboid(4.0, 4.0, 4.0);

        boost::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
        p_grid->GenerateFromPart(p_domain, 0.5*1.e-6*unit::metres);

        // Choose the PDE
        boost::shared_ptr<HybridNonLinearEllipticPde<3> > p_non_linear_pde = HybridNonLinearEllipticPde<3>::Create();
        p_non_linear_pde->SetIsotropicDiffusionConstant(1.0);
        p_non_linear_pde->SetContinuumConstantInUTerm(-2.0);
        p_non_linear_pde->SetThreshold(0.0625);

        // Choose the Boundary conditions
        boost::shared_ptr<HybridBoundaryCondition<3> > p_outer_boundary_condition = HybridBoundaryCondition<3>::Create();
        p_outer_boundary_condition->SetValue(1.0);

        FiniteDifferenceSolver<3> solver;
        solver.SetGrid(p_grid);
        solver.SetNonLinearPde(p_non_linear_pde);
        solver.AddBoundaryCondition(p_outer_boundary_condition);

        MAKE_PTR_ARGS(OutputFileHandler, p_output_file_handler, ("TestNonLinearFiniteDifferenceSolver/Box", false));
        solver.SetFileHandler(p_output_file_handler);
        solver.SetFileName("output_nl_fd.vti");
        solver.SetWriteSolution(true);
        solver.Solve();
    }
};

#endif /*TESTFINITEELEMENTSOLVER_HPP_*/
