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
#include "FiniteElementSolver.hpp"
#include "UblasIncludes.hpp"
#include "Part.hpp"
#include "Vertex.hpp"
#include "HybridLinearEllipticPde.hpp"
#include "HybridNonLinearEllipticPde.hpp"
#include "CaVascularNetwork.hpp"
#include "VasculatureGenerator.hpp"
#include "SmartPointers.hpp"
#include "VasculatureData.hpp"
#include "HybridBoundaryCondition.hpp"
#include "HybridMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "Debug.hpp"

class TestNonLinearFiniteElementSolver : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestBox() throw(Exception)
    {
        // Set up the mesh
        boost::shared_ptr<Part<3> > p_domain = Part<3>::Create();
        p_domain->AddCuboid(200.0, 200.0, 200.0);
        boost::shared_ptr<HybridMesh<3, 3> > p_mesh = HybridMesh<3, 3>::Create();
        p_mesh->GenerateFromPart(p_domain, 800);

        // Choose the PDE
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_linear_pde = HybridLinearEllipticPde<3>::Create();
        p_linear_pde->SetIsotropicDiffusionConstant(0.0033);
        p_linear_pde->SetContinuumConstantInUTerm(-8.e-5);

        boost::shared_ptr<HybridNonLinearEllipticPde<3> > p_non_linear_pde = HybridNonLinearEllipticPde<3>::Create();
        p_non_linear_pde->SetIsotropicDiffusionConstant(0.0033);
        p_non_linear_pde->SetContinuumConstantInUTerm(-8.e-5);

        // Choose the Boundary conditions
        boost::shared_ptr<HybridBoundaryCondition<3> > p_outer_boundary_condition = HybridBoundaryCondition<3>::Create();
        p_outer_boundary_condition->SetValue(40.0);

        FiniteElementSolver<3> solver;
        solver.SetMesh(p_mesh);
        solver.SetPde(p_linear_pde);
        solver.SetNonLinearPde(p_non_linear_pde);
        solver.AddBoundaryCondition(p_outer_boundary_condition);

        MAKE_PTR_ARGS(OutputFileHandler, p_output_file_handler, ("TestNonLinearFiniteElementSolver/Box", false));
        solver.SetFileHandler(p_output_file_handler);
        solver.SetFileName("output_nl");
        solver.SetWriteSolution(true);
        solver.SetUseSimpleNetonSolver(true);
        solver.Solve();
    }
};

#endif /*TESTFINITEELEMENTSOLVER_HPP_*/
