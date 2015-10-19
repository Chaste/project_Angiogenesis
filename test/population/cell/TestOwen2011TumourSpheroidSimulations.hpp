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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */

#ifndef TESTOwen2011TUMOURSPHEROIDSIMULATIONS_HPP_
#define TESTOwen2011TUMOURSPHEROIDSIMULATIONS_HPP_

/*
 * = An example showing how to run tumour spheroid simulations =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to simulate a growing cell monolayer culture or
 * multicellular tumour spheroid. Like the crypt simulations, tumour spheroid simulations
 * include cell-cycle models and force laws to determine how cells divide and
 * move. In tumour spheroid simulations, however, these are also coupled to a
 * system of partial differential equations (PDEs) that determine the concentration
 * of specified nutrients (e.g. oxygen) throughout the cell population. Also, unlike
 * in a crypt simulation (for example), the cell population may grow substantially as the simulation
 * progresses.
 *
 * In summary, the main difference between this tutorial and the other cell-based simulation
 * tutorials is that a PDE is defined, which is used in the simulation.
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * As in the other cell-based simulation tutorials, we begin by including the necessary header files. We have
 * encountered some of these files already. Recall that often {{{CheckpointArchiveTypes.hpp}}}
 * or {{{CellBasedSimulationArchiver.hpp}}} must be included as the first Chaste header.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "RandomNumberGenerator.hpp"
#include "SmartPointers.hpp"
#include "ApoptoticCellKiller.hpp"
#include "AveragedSourcePde.hpp"
/*
 * The {{{SimpleOxygenBasedCellCycleModel}}} header file defines a cell-cycle model in which
 * a cell's rate of progress through G1 phase changes over time in a simple manner, according
 * to the local oxygen concentration. We also include the {{{WildTypeCellMutationState}}}
 * header file, which defines a wild type cell mutation state that we will use to construct
 * cells. A cell mutation state is always required when constructing a cell, however
 * in earlier simulation tutorial we used a helper classes (({{{CellsGenerator}}} and {{{CryptCellsGenerator}}}) that
 * allowed us to avoid having to construct cells directly.
 */
#include "Owen2011OxygenBasedCellCycleModel.hpp"
#include "Owen2011OxygenBasedCellCycleModelWithoutOde.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CancerCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"

/*
 * The next three header files define: a PDE that describes how oxygen is transported via through the
 * domain via diffusion and is consumed by live cells; a constant-valued boundary condition to
 * associate with the PDE; and a PDE handler class, which is passed to the simulation object and
 * handles the numerical solution of any PDEs.
 */
#include "CellwiseSourcePde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "CellBasedPdeHandler.hpp"

/*
 * We use an {{{OffLatticeSimulation}}}.
 */
#include "OffLatticeSimulation.hpp"
#include "VoronoiDataWriter.hpp"

#include "Owen2011TrackingModifier.hpp"

/*
 * The header file {{{PetscSetupAndFinalize.hpp}}} must be included in all tests which use Petsc. This is
 * a suite of data structures and routines that are used in the finite element
 * PDE solvers, which is how we solve the oxygen transport PDE.
 */
#include "PetscSetupAndFinalize.hpp"
#include "Debug.hpp"
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

/*
 * Having included all the necessary header files, we proceed by defining the test class.
 */
class TestOwen2011TumourSpheroidSimulations : public AbstractCellBasedTestSuite
{
public:

    void TestSimpleSpheroidOwen2011OxygenBasedCellCycleModel()
    {

        Timer::Reset();

        /*
         * First we want to create a '''non-periodic''' 'honeycomb' mesh.
         * We use the honeycomb mesh generator, as before, saying 10 cells wide
         * and 10 cells high. Note that the thickness of the ghost nodes layer is
         * 0, i.e. there are no ghost nodes, and the {{{false}}} indicates that the
         * returned mesh is '''not''' cylindrical. In contrast to the crypt simulation
         * tutorial, here we call {{{GetMesh()}}} on the {{{HoneycombMeshGenerator}}}
         * object to return the mesh, which is of type {{{MutableMesh}}}.
         */
        HoneycombMeshGenerator generator(20, 20, 0);
        MutableMesh<2,2>* p_mesh = generator.GetCircularMesh(2);

        /*
         * Next, we need to create some cells. Unlike in the the crypt simulation
         * tutorial, we don't just use a {{{CellsGenerator}}} class, but do it manually,
         * in a loop. First, we define a {{{std::vector}}} of cell pointers.
         */
        std::vector<CellPtr> cells;

        /*
         * This line defines a mutation state to be used for all cells, of type
         * `WildTypeCellMutationState` (i.e. 'healthy'):
         */
        MAKE_PTR(CancerCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        // Set up oxygen_concentration
        double oxygen_concentration = 20.0;

        /*
         * Now we loop over the nodes...
         */
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {

            /*
             * ...then create a cell, giving it a {{{SimpleOxygenBasedCellCycleModel}}}.
             * The spatial dimension (1, 2 or 3) needs to be set on the cell-cycle model before it is passed to the cell.
             */
            Owen2011OxygenBasedCellCycleModelWithoutOde* const p_model = new Owen2011OxygenBasedCellCycleModelWithoutOde;
            p_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->SetApoptosisTime(30);

            /*
             * We also alter the default cell-cycle times.
             */
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            /*
             * We now define a random birth time, chosen from [-T,0], where
             * T = t,,1,, + t,,2,,, where t,,1,, is a parameter representing the G,,1,, duration
             * of a 'stem' cell, and t,,2,, is the basic S+G,,2,,+M phases duration...
             */
            double birth_time = RandomNumberGenerator::Instance()->ranf() * 26.67; // time needed for cancer cell to complete cycle when well oxygenated according to Owen2011
            /*
             * ...then we set the birth time and push the cell back into the vector
             * of cells.
             */
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);

            p_cell->GetCellData()->SetItem("oxygen", oxygen_concentration);

        }

        /*
         * Now that we have defined the cells, we can define the {{{CellPopulation}}}. We use a
         * {{{MeshBasedCellPopulation}}} since although the cell population is mesh-based, it does
         * not include any ghost nodes. The constructor takes in the mesh and the cells vector.
         */
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputResultsForChasteVisualizer(false);
        cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        //        cell_population.SetWriteVtkAsPoints(false);
        //        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        AveragedSourcePde<2> pde(cell_population, -0.0359);
        ConstBoundaryCondition<2> bc(oxygen_concentration);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("oxygen");
        pde_handler.AddPdeAndBc(&pde_and_bc);


        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(10.0, 10.0);
        ChasteCuboid<2> cuboid(lower, upper);
        pde_handler.UseCoarsePdeMesh(2.0, cuboid, true);
        pde_handler.SetImposeBcsOnCoarseBoundary(true);

        /*
         * We are now in a position to construct an {{{OffLatticeSimulationWithPdes}}} object,
         * using the cell population. We then pass the PDE handler object to the simulation.
         */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetCellBasedPdeHandler(&pde_handler);

        /*
         * We next set the output directory and end time.
         */
        std::string resultsDirectoryName = "TestOwen2011TumourSpheroidGrowth";
        simulator.SetOutputDirectory(resultsDirectoryName);
        simulator.SetSamplingTimestepMultiple(60); // 60 saves a point every 30min
        simulator.SetEndTime(50.0);

        /*
         * Create cell killer to remove apoptotic cell from simulation
         */
        boost::shared_ptr<ApoptoticCellKiller<2> > apoptotic_cell_killer(new ApoptoticCellKiller<2>(&cell_population));
        simulator.AddCellKiller(apoptotic_cell_killer);


        // Create a Cell Concentration tracking modifier and add it to the simulation
        MAKE_PTR(Owen2011TrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        /*
         * We must now create one or more force laws, which determine the mechanics of
         * the cell population. As in the crypt simulation tutorial, we assume that a cell
         * experiences a force from each neighbour that can be represented as a linear overdamped
         * spring, so we use a {{{GeneralisedLinearSpringForce}}} object.
         * Note that we have called the method {{{SetCutOffLength}}} on the
         * {{{GeneralisedLinearSpringForce}}} before passing it to the simulator: this call
         * modifies the force law so that two neighbouring cells do not impose
         * a force on each other if they are located more than 3 units (=3 cell widths)
         * away from each other. This modification is necessary when no ghost nodes are used,
         * for example to avoid artificially large forces between cells that lie close together
         * on the spheroid boundary.
         */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(3);
        simulator.AddForce(p_linear_force);

        /*
         * We call {{{Solve()}}} on the simulator to run the simulation.
         */
        simulator.Solve();
        Timer::Print("Elapsed time = ");
    }



};

#endif /*TESTOwen2011TUMOURSPHEROIDSIMULATIONS_HPP_*/
