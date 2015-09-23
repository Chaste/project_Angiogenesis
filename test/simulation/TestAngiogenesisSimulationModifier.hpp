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

#ifndef TESTSPHEROIDWITHANGIOGENESIS_HPP_
#define TESTSPHEROIDWITHANGIOGENESIS_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PottsMeshGenerator.hpp"
#include "Cell.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "CaBasedCellPopulation.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "OnLatticeSimulation.hpp"
#include "Node.hpp"
#include "NodesOnlyMesh.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AngiogenesisModifier.hpp"
#include "FiniteDifferenceSolver.hpp"
#include "HybridLinearEllipticPde.hpp"
#include "VasculatureGenerator.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "CellLabelWriter.hpp"

class TestSpheroidWithAngiogenesis : public AbstractCellBasedTestSuite
{
public:

    void TestCaBasedSpheroid() throw (Exception)
    {
        // Create the domain
        double domain_x = 800.0;
        double domain_y = 800.0;
        double domain_z = 400.0;
        boost::shared_ptr<Part<3> > p_domain = Part<3> ::Create();
        p_domain->AddCuboid(domain_x, domain_y, domain_z);

        // Create the vessels
        VasculatureGenerator<3> vessel_generator;
        c_vector<double, 3> vessel_start_point;
        vessel_start_point[0] = 400.0;
        vessel_start_point[1] = 400.0;
        vessel_start_point[2] = 0.0;
        boost::shared_ptr<CaVascularNetwork<3> > p_network = vessel_generator.GenerateSingleVessel(400.0, vessel_start_point);

        // Create a Potts mesh
        double spacing = 40.0;
        unsigned num_x = unsigned(domain_x/spacing) + 1;
        unsigned num_y = unsigned(domain_y/spacing) + 1;
        unsigned num_z = unsigned(domain_z/spacing) + 1;
        PottsMeshGenerator<3> generator(num_x, 0, 0, num_y, 0, 0, num_z, 0, 0);
        PottsMesh<3>* p_mesh = generator.GetMesh();
        p_mesh->Scale(spacing, spacing, spacing);

        // Create cells in a cylinder
        c_vector<double,3> centre;
        centre[0] = domain_x/2.0;
        centre[1] = domain_y/2.0;
        centre[2] = 0.0;
        double radius = 100.0;
        std::vector<unsigned> location_indices;
        for(unsigned kdx=0; kdx<num_z; kdx++)
        {
            for(unsigned jdx=0; jdx<num_y; jdx++)
            {
                for(unsigned idx=0; idx<num_x; idx++)
                {
                    unsigned location_index = idx + num_x * jdx + num_x * num_y * kdx;
                    c_vector<double,3> location;
                    location[0] = double(idx) * spacing;
                    location[1] = double(jdx) * spacing;
                    location[2] = double(kdx) * spacing;

                    if(norm_2(location - centre) <= radius && location[2] >= 40.0 && location[2] <= 80.0)
                    {
                        location_indices.push_back(location_index);
                    }
                }
            }
        }

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size());

        // Create cell population
        CaBasedCellPopulation<3> cell_population(*p_mesh, cells, location_indices);
        cell_population.SetCellAncestorsToLocationIndices();

        // Create a pde and solver
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_pde = HybridLinearEllipticPde<3>::Create();
        p_pde->SetDiffusionConstant(0.002);
        double consumptionRate = 0.0;
        p_pde->SetLinearInUTerm(consumptionRate);

        boost::shared_ptr<FiniteDifferenceSolver<3> > p_solver = boost::shared_ptr<FiniteDifferenceSolver<3> >(new FiniteDifferenceSolver<3>());
        p_solver->SetExtents(p_domain, 40.0);
        p_solver->SetPde(p_pde);
        p_solver->SetVesselNetwork(p_network);
        p_solver->SetBoundaryConditionType(BoundaryConditionType::LINE);
        p_solver->SetBoundaryConditionSource(BoundaryConditionSource::USER);
        p_solver->SetBoundaryConditionValue(40.0);
        p_solver->SetBoundaryConditionName("Oxygen");

        boost::shared_ptr<AngiogenesisModifier<3> > p_modifier = boost::shared_ptr<AngiogenesisModifier<3> >(new AngiogenesisModifier<3>);
        p_modifier->SetSolver(p_solver);

        OnLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestSpheroidAngiogenesis/CaBased");
        simulator.SetDt(1.0);
        simulator.SetEndTime(10.0);
        simulator.AddSimulationModifier(p_modifier);
        simulator.Solve();
    }

    void DOntTestNodeBasedSpheroid() throw (Exception)
    {
        // Create the domain
        double domain_x = 800.0;
        double domain_y = 800.0;
        double domain_z = 400.0;
        boost::shared_ptr<Part<3> > p_domain = Part<3> ::Create();
        p_domain->AddCuboid(domain_x, domain_y, domain_z);

        // Create the vessels
        VasculatureGenerator<3> generator;
        c_vector<double, 3> vessel_start_point;
        vessel_start_point[0] = 400.0;
        vessel_start_point[1] = 400.0;
        vessel_start_point[2] = 0.0;
        boost::shared_ptr<CaVascularNetwork<3> > p_network = generator.GenerateSingleVessel(400.0, vessel_start_point);

        // Create cells in a cylinder
        double spacing = 40.0;
        double pde_spacing = 40.0;
        unsigned num_x = unsigned(domain_x/spacing) + 1;
        unsigned num_y = unsigned(domain_y/spacing) + 1;
        unsigned num_z = unsigned(domain_z/spacing) + 1;
        std::vector<Node<3>*> nodes;
        c_vector<double,3> centre;
        centre[0] = domain_x/2.0;
        centre[1] = domain_y/2.0;
        centre[2] = 0.0;
        double radius = 100.0;

        unsigned counter = 0;
        for(unsigned kdx=0; kdx<num_z; kdx++)
        {
            for(unsigned jdx=0; jdx<num_y; jdx++)
            {
                for(unsigned idx=0; idx<num_x; idx++)
                {
                    c_vector<double,3> location;
                    location[0] = double(idx) * spacing;
                    location[1] = double(jdx) * spacing;
                    location[2] = 0.0;

                    double position_z = double(kdx) * spacing;
                    if(norm_2(location - centre) <= radius && position_z >= 160.0 && position_z <= 240.0)
                    {
                        nodes.push_back(new Node<3>(counter,  false, location[0], location[1], position_z));
                        counter++;
                    }
                }
            }
        }

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5 * spacing);
        std::vector<CellPtr> cells;

        CellsGenerator<SimpleOxygenBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(2.0 * spacing);
        cell_population.AddCellWriter<CellLabelWriter>();

        // Create a pde and solver
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_pde = HybridLinearEllipticPde<3>::Create();
        p_pde->SetDiffusionConstant(0.002);
        double consumptionRate = 0.0;
        p_pde->SetLinearInUTerm(consumptionRate);

        boost::shared_ptr<FiniteDifferenceSolver<3> > p_solver = boost::shared_ptr<FiniteDifferenceSolver<3> >(new FiniteDifferenceSolver<3>());
        p_solver->SetExtents(p_domain, pde_spacing);
        p_solver->SetPde(p_pde);
        p_solver->SetVesselNetwork(p_network);
        p_solver->SetBoundaryConditionType(BoundaryConditionType::LINE);
        p_solver->SetBoundaryConditionSource(BoundaryConditionSource::USER);
        p_solver->SetBoundaryConditionValue(40.0);
        p_solver->SetBoundaryConditionName("Oxygen");

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestSpheroidAngiogenesis/NodeBased");
        simulator.SetDt(1.0);
        simulator.SetEndTime(100.0);

        boost::shared_ptr<AngiogenesisModifier<3> > p_modifier = boost::shared_ptr<AngiogenesisModifier<3> >(new AngiogenesisModifier<3>);
        p_modifier->SetSolver(p_solver);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        simulator.AddForce(p_force);
        simulator.Solve();

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
};

#endif /* TESTSPHEROIDWITHANGIOGENESIS_HPP_ */
