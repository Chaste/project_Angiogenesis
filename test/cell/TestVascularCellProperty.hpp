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

#ifndef TESTVASCULARCELLPROPERTY_HPP_
#define TESTVASCULARCELLPROPERTY_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
/* Most Chaste code uses PETSc to solve linear algebra problems.  This involves starting PETSc at the beginning of a test-suite
 * and closing it at the end.  (If you never run code in parallel then it is safe to replace PetscSetupAndFinalize.hpp with FakePetscSetup.hpp)
 */
#include "VascularCellProperty.hpp"
/* The remaining header files define classes that will be used in the cell-based
 * simulation test. We have encountered each of these header files in previous cell-based
 * Chaste tutorials. */
#include "AbstractForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * @file
 *
 * This is an example of a CxxTest test suite, used to test the source
 * code, and also used to run simulations (as it provides a handy
 * shortcut to compile and link against the correct libraries using scons).
 *
 * You can #include any of the files in the project 'src' folder.
 * For example here we #include "Hello.hpp"
 *
 * You can utilise any of the code in the main the Chaste trunk
 * in exactly the same way.
 * NOTE: you will have to alter the project SConscript file lines 41-44
 * to enable #including of code from the 'heart', 'cell_based' or 'crypt'
 * components of Chaste.
 */

class TestCreatingAndUsingANewCellPropertyTutorial : public AbstractCellBasedTestSuite
{
public:

    /*
     * === Testing the cell property ===
     *
     * We begin by testing that our new cell property is implemented correctly.
     */
    void TestVascularCellProperty() throw(Exception)
    {
        /* We begin by testing that some of the base class methods work correctly.
         * We typically use shared pointers to create and access a cell property
         * like {{{VascularCellProperty}}}, for which it makes sense for all cells
         * that have the same mutation to share a pointer to the same cell property
         * object (although strictly speaking, they are not required to). Observe that
         * in this case we have provided a value for the member variable {{{mColour}}}
         * in the {{{VascularCellProperty}}} constructor.*/
        MAKE_PTR_ARGS(VascularCellProperty, p_property, (8));

        /* Each cell property has a member variable, {{{mCellCount}}}, which
         * stores the number of cells with this cell property. We can test whether
         * {{{mCellCount}}} is being updated correctly by our cell property, as follows. */
        TS_ASSERT_EQUALS(p_property->GetCellCount(), 0u);
        p_property->IncrementCellCount();
        TS_ASSERT_EQUALS(p_property->GetCellCount(), 1u);
        p_property->DecrementCellCount();
        TS_ASSERT_EQUALS(p_property->GetCellCount(), 0u);
        TS_ASSERT_THROWS_THIS(p_property->DecrementCellCount(),
                "Cannot decrement cell count: no cells have this cell property");

        /* We can also test whether our cell property is of a given type, as follows. */
        TS_ASSERT_EQUALS(p_property->IsType<WildTypeCellMutationState>(), false);
        TS_ASSERT_EQUALS(p_property->IsType<VascularCellProperty>(), true);

        /* We can also test that archiving is implemented correctly for our cell
         * property, as follows (further details on how to implement and
         * test archiving can be found at ChasteGuides/BoostSerialization).  */
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "property.arch";

        {
            AbstractCellProperty* const p_const_property = new VascularCellProperty(7);
            p_const_property->IncrementCellCount();

            TS_ASSERT_EQUALS(p_const_property->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(dynamic_cast<VascularCellProperty*>(p_const_property)->GetColour(), 7u);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_const_property;

            delete p_const_property;
        }

        {
            AbstractCellProperty* p_arch_property;

            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_arch_property;

            TS_ASSERT_EQUALS(p_arch_property->GetCellCount(), 1u);

            VascularCellProperty* p_real_property = dynamic_cast<VascularCellProperty*>(p_arch_property);
            TS_ASSERT(p_real_property != NULL);
            TS_ASSERT_EQUALS(p_real_property->GetColour(), 7u);

            delete p_arch_property;
        }
    }

//    /*
//     * === Using the cell property in a cell-based simulation ===
//     *
//     * We conclude with a brief test demonstrating how {{{VascularCellProperty}}} can be used
//     * in a cell-based simulation.
//     */
//    void TestOffLatticeSimulationWithVascularCellProperty() throw(Exception)
//    {
//        /* Note that HoneycombMeshGenerator, used in this test, is not
//         *  yet implemented in parallel. */
//
//        /* We use the {{{HoneycombMeshGenerator}}} to create a honeycomb mesh covering a
//         * circular domain of given radius, and use this to generate a {{{NodesOnlyMesh}}}
//         * as follows. */
//        HoneycombMeshGenerator generator(10, 10);
//        MutableMesh<2,2>* p_generating_mesh = generator.GetCircularMesh(5);
//
//        NodesOnlyMesh<2> mesh;
//        /* We construct the mesh using the generating mesh and a cut-off 1.5 which defines the
//         * connectivity in the mesh.
//         */
//        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);
//
//        /* We now create a shared pointer to our new property, as follows. */
//        MAKE_PTR(VascularCellProperty, p_Vascular);
//        /*
//         * Also create a shared pointer to a cell label so we can visualize the
//         * different cell types. Note that this is also a {{{CellProperty}}}.
//         */
//        MAKE_PTR(CellLabel, p_label);
//
//        /* Next, we create some cells. We don't use a Generator as we want to give some cells the new cell property, therefore
//         * we create the cells in a loop, as follows.*/
//        MAKE_PTR(WildTypeCellMutationState, p_state);
//        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
//        std::vector<CellPtr> cells;
//        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
//        {
//            /* For each node we create a cell with our cell-cycle model and the wild-type cell mutation state.
//             * We then add the property {{{VascularCellProperty}}} to a random selection of the cells, as follows. */
//            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
//
//            CellPropertyCollection collection;
//            if (RandomNumberGenerator::Instance()->ranf() < 0.2)
//            {
//                collection.AddProperty(p_Vascular);
//                collection.AddProperty(p_label);
//            }
//
//            CellPtr p_cell(new Cell(p_state, p_model, false, collection));
//            p_cell->SetCellProliferativeType(p_diff_type);
//
//            /* Now, we define a random birth time, chosen from [-T,0], where
//             * T = t,,1,, + t,,2,,, where t,,1,, is a parameter representing the G,,1,, duration
//             * of a stem cell, and t,,2,, is the basic S+G,,2,,+M phases duration.
//             */
//            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
//                                    (p_model->GetStemCellG1Duration()
//                                        + p_model->GetSG2MDuration());
//
//            /* Finally, we set the birth time and push the cell back into the vector of cells. */
//            p_cell->SetBirthTime(birth_time);
//            cells.push_back(p_cell);
//        }
//
//        /* Now that we have defined the mesh and cells, we can define the cell population. The constructor
//         * takes in the mesh and the cells vector. */
//        NodeBasedCellPopulation<2> cell_population(mesh, cells);
//
//        /* In order to visualize labelled cells we need to use the following command.*/
//        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
//
//        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
//         * and set the output directory, output multiple, and end time. */
//        OffLatticeSimulation<2> simulator(cell_population);
//        simulator.SetOutputDirectory("TestOffLatticeSimulationWithVascularCellProperty");
//        simulator.SetSamplingTimestepMultiple(12);
//        simulator.SetEndTime(10.0);
//
//        /* We create a force law and pass it to the {{{OffLatticeSimulation}}}. */
//        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
//        p_linear_force->SetCutOffLength(1.5);
//        simulator.AddForce(p_linear_force);
//
//        /* Now create a {{{MotlieForce}}} and pass it to the {{{OffLatticeSimulation}}}. */
//        MAKE_PTR(MyMotiveForce, p_motive_force);
//        simulator.AddForce(p_motive_force);
//
//        /* To run the simulation, we call {{{Solve()}}}. */
//        simulator.Solve();
//    }
};
/*
 * When you visualize the results with
 *
 * {{{java Visualize2dCentreCells /tmp/$USER/testoutput/TestOffLatticeSimulationWithVascularCellProperty/results_from_time_0}}}
 *
 * you should see a collection of cells with the {{{VascularCellProperty}}} (labelled dark blue) moving towards the origin.
 */

#endif /*TESTVASCULARCELLPROPERTY_HPP_*/
