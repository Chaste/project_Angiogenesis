/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTSOLVINGPDES_HPP_
#define TESTSOLVINGPDES_HPP_

/* = Solving PDEs in the Angiogenesis Project =
 * This tutorial demonstrates methods for solving PDEs in the Angiogenesis Project. It is noted
 * that the way to set up PDEs differs from that of Cell Based Chaste, although the same solver
 * can be used behind the scenes.
 *
 * The following is covered:
 * * Solving a linear reaction-diffusion PDE with finite differences.
 * * Solving a linear reaction-diffusion PDE with finite finite elements and discrete sinks and sources.
 * * Solving a non-linear reaction-diffusion PDE with finite differences.
 * * Solving a linear reaction-diffusion PDE using Green's functions.
 * * Interacting with regular grids and finite element meshes.
 *
 * = The Test =
 * Start by introducing the necessary header files. The first contain functionality for setting up unit tests,
 * smart pointer tools and output management.
 */
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "SmartPointers.hpp"
#include "OutputFileHandler.hpp"
#include "FileFinder.hpp"
/*
 * Dimensional analysis.
 */
#include "DimensionalChastePoint.hpp"
#include "UnitCollection.hpp"
#include "Owen11Parameters.hpp"
#include "GenericParameters.hpp"
#include "ParameterCollection.hpp"
#include "BaseUnits.hpp"
/*
 * Geometry tools.
 */
#include "MappableGridGenerator.hpp"
#include "Part.hpp"
/*
 * Vessel networks.
 */
#include "VesselNetwork.hpp"
#include "VesselNetworkGenerator.hpp"
/*
 * Grids and PDEs.
 */
#include "DiscreteContinuumMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "FiniteElementSolver.hpp"
#include "FiniteDifferenceSolver.hpp"
#include "DiscreteSource.hpp"
#include "VesselBasedDiscreteSource.hpp"
#include "DiscreteContinuumBoundaryCondition.hpp"
#include "LinearSteadyStateDiffusionReactionPde.hpp"
#include "MichaelisMentenSteadyStateDiffusionReactionPde.hpp"
/*
 * This should appear last.
 */
#include "PetscSetupAndFinalize.hpp"
class TestSolvingPdesLiteratePaper : public AbstractCellBasedWithTimingsTestSuite
{
public:
    /*
     * = Test 1 - Linear Reaction Diffusion With Finite Differences =
     * In the first example we will solve a steady-state linear reaction diffusion
     * PDE with finite differences.
     */
    void TestLinearReactionDiffusionPdeWithFiniteDifferences() throw(Exception)
    {
        MAKE_PTR_ARGS(OutputFileHandler, p_handler, ("TestSolvingPdesLiteratePaper/TestLinearReactionDiffusionPdeWithFiniteDifferences"));
        /*
         * We will work in microns
         */
        units::quantity<unit::length> reference_length(1.0 * unit::microns);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);
        /*
         * Set up a simulation domain, which will be a cuboid.
         */
        units::quantity<unit::length> domain_width(100.0 * 1.e-6 * unit::microns);
        units::quantity<unit::length> domain_height(100.0 * 1.e-6 * unit::microns);
        units::quantity<unit::length> domain_depth(20.0 * 1.e-6 * unit::microns);
        boost::shared_ptr<Part<3> > p_domain = Part<3>::Create();
        p_domain->AddCuboid(domain_width, domain_height, domain_depth, DimensionalChastePoint<3>(0.0, 0.0, 0.0));
        /*
         * Make a regular grid on the domain
         */
        boost::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
        p_grid->GenerateFromPart(p_domain, 10.0*reference_length);
        /*
         * Set up a PDE, we will model oxygen diffusion.
         */
        boost::shared_ptr<LinearSteadyStateDiffusionReactionPde<3> > p_oxygen_pde = LinearSteadyStateDiffusionReactionPde<3>::Create();
        units::quantity<unit::diffusivity> oxygen_diffusivity(1.e-6*unit::metre_squared_per_second);
        p_oxygen_pde->SetIsotropicDiffusionConstant(oxygen_diffusivity);
        /*
         * Add continuum sink term for cells
         */
        units::quantity<unit::rate> oxygen_consumption_rate(1.e-6*unit::per_second);
        p_oxygen_pde->SetContinuumLinearInUTerm(-oxygen_consumption_rate);
        /*
         * Add a Dirichlet boundary condition on the left face of the domain.
         */
        p_domain->GetFacet(DimensionalChastePoint<3>(0.0,
                                                     domain_height/(2.0*reference_length),
                                                     domain_depth/(2.0*reference_length)))->SetLabel("boundary_1");
        boost::shared_ptr<DiscreteContinuumBoundaryCondition<3> > p_left_face_boundary = DiscreteContinuumBoundaryCondition<3>::Create();
        p_left_face_boundary->SetType(BoundaryConditionType::FACET);
        p_left_face_boundary->SetDomain(p_domain);
        p_left_face_boundary->SetValue(10.0*unit::mole_per_metre_cubed);
        p_left_face_boundary->SetLabelName("boundary_1");
        /*
         * Set up the PDE solvers for the oxygen problem
         */
        boost::shared_ptr<FiniteDifferenceSolver<3> > p_oxygen_solver = FiniteDifferenceSolver<3>::Create();
        p_oxygen_solver->SetPde(p_oxygen_pde);
        p_oxygen_solver->SetGrid(p_grid);
        p_oxygen_solver->AddBoundaryCondition(p_left_face_boundary);
        p_oxygen_solver->SetLabel("oxygen");
        p_oxygen_solver->SetFileHandler(p_handler);
        p_oxygen_solver->SetFileName("fd_solution.vti");
        p_oxygen_solver->SetWriteSolution(true);
        p_oxygen_solver->Setup();
        p_oxygen_solver->Solve();
    }
};

#endif /*TESTSOLVINGPDES_HPP_*/
