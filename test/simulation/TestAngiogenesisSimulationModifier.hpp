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
#include "Polygon.hpp"
#include "NodesOnlyMesh.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "FakePetscSetup.hpp"
#include "AngiogenesisModifier.hpp"
#include "FiniteDifferenceSolver.hpp"
#include "HybridLinearEllipticPde.hpp"
#include "VasculatureGenerator.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "DirichletBoundaryCondition.hpp"
#include "DiscreteSource.hpp"
#include "CellLabelWriter.hpp"
#include "CaVessel.hpp"
#include "VascularNode.hpp"
#include "GeometryTools.hpp"
#include "OffLatticePrwGrowthDirectionModifier.hpp"
#include "OffLatticeTipAttractionGrowthDirectionModifier.hpp"
#include "OffLatticeSolutionDependentGrowthDirectionModifier.hpp"
#include "OffLatticeRandomNormalSproutingRule.hpp"
#include "RandomNumberGenerator.hpp"

class TestSpheroidWithAngiogenesis : public AbstractCellBasedTestSuite
{

    boost::shared_ptr<Part<3> > GetSimulationDomain()
    {
        double domain_x = 800.0;
        double domain_y = 800.0;
        double domain_z = 200.0;
        boost::shared_ptr<Part<3> > p_domain = Part<3> ::Create();
        p_domain->AddCuboid(domain_x, domain_y, domain_z);
        return p_domain;
    }

    boost::shared_ptr<CaVascularNetwork<3> > GetVesselNetwork()
    {
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();

        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        for(unsigned idx=0; idx<81; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 50.0, 100.0));
        }
        boost::shared_ptr<CaVessel<3> > p_vessel_1 = CaVessel<3>::Create(bottom_nodes);

        std::vector<boost::shared_ptr<VascularNode<3> > > top_nodes;
        for(unsigned idx=0; idx<81; idx++)
        {
            top_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 750.0, 100.0));
        }
        boost::shared_ptr<CaVessel<3> > p_vessel_2 = CaVessel<3>::Create(top_nodes);

        // Set up flow properties
        p_network->AddVessel(p_vessel_1);
        p_network->AddVessel(p_vessel_2);
        p_network->GetVessel(0)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
        p_network->GetVessel(0)->GetStartNode()->GetFlowProperties()->SetPressure(3000.0);
        p_network->GetVessel(0)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
        p_network->GetVessel(0)->GetEndNode()->GetFlowProperties()->SetPressure(1000.0);

        p_network->GetVessel(1)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
        p_network->GetVessel(1)->GetStartNode()->GetFlowProperties()->SetPressure(3000.0);
        p_network->GetVessel(1)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
        p_network->GetVessel(1)->GetEndNode()->GetFlowProperties()->SetPressure(1000.0);

        p_network->UpdateSegments();
        p_network->SetSegmentRadii(10.0);
        std::vector<boost::shared_ptr<CaVesselSegment<3> > > segments = p_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(1.e-9);
        }
        return p_network;
    }

    boost::shared_ptr<Part<3> > GetInitialTumourCellRegion()
    {
        double radius = 100.0;
        double depth = 200.0;
        c_vector<double, 3> origin;
        origin[0] = 400.0;
        origin[1] = 400.0;
        origin[2] = 0.0;

        boost::shared_ptr<Part<3> > p_domain = Part<3> ::Create();
        boost::shared_ptr<Polygon> circle = p_domain->AddCircle(radius, origin);
        p_domain->Extrude(circle, depth);
        return p_domain;
    }

    boost::shared_ptr<FiniteDifferenceSolver<3> > GetOxygenSolver(boost::shared_ptr<Part<3> > p_domain,
                                                                  boost::shared_ptr<CaVascularNetwork<3> > p_network)
    {
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_oxygen_pde = HybridLinearEllipticPde<3>::Create();
        p_oxygen_pde->SetDiffusionConstant(0.0033);
        p_oxygen_pde->SetVariableName("oxygen");

        boost::shared_ptr<DiscreteSource<3> > p_cell_oxygen_sink = DiscreteSource<3>::Create();
        p_cell_oxygen_sink->SetType(SourceType::MULTI_POINT);
        p_cell_oxygen_sink->SetSource(SourceStrength::PRESCRIBED);
        p_cell_oxygen_sink->SetValue(1.e-6);
        p_cell_oxygen_sink->SetIsLinearInSolution(true);
        p_oxygen_pde->AddDiscreteSource(p_cell_oxygen_sink);

        boost::shared_ptr<DirichletBoundaryCondition<3> > p_vessel_ox_boundary_condition = DirichletBoundaryCondition<3>::Create();
        p_vessel_ox_boundary_condition->SetValue(40.0);
        p_vessel_ox_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
        p_vessel_ox_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

        boost::shared_ptr<FiniteDifferenceSolver<3> > p_oxygen_solver = FiniteDifferenceSolver<3>::Create();
        p_oxygen_solver->SetExtents(p_domain, 40.0);
        p_oxygen_solver->SetPde(p_oxygen_pde);
        p_oxygen_solver->SetVesselNetwork(p_network);
        p_oxygen_solver->AddDirichletBoundaryCondition(p_vessel_ox_boundary_condition);

        return p_oxygen_solver;
    }

    boost::shared_ptr<FiniteDifferenceSolver<3> > GetVegfSolver(boost::shared_ptr<Part<3> > p_domain,
                                                                  boost::shared_ptr<CaVascularNetwork<3> > p_network)
    {
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_vegf_pde = HybridLinearEllipticPde<3>::Create();
        p_vegf_pde->SetDiffusionConstant(0.0033);
        p_vegf_pde->SetVariableName("vegf");
        p_vegf_pde->SetLinearInUTerm(-1.e-7);

        boost::shared_ptr<DiscreteSource<3> > p_cell_vegf_source = DiscreteSource<3>::Create();
        p_cell_vegf_source->SetType(SourceType::MULTI_POINT);
        p_cell_vegf_source->SetSource(SourceStrength::PRESCRIBED);
        p_cell_vegf_source->SetValue(-1.e-4);
        p_cell_vegf_source->SetIsLinearInSolution(false);
        p_vegf_pde->AddDiscreteSource(p_cell_vegf_source);

        boost::shared_ptr<FiniteDifferenceSolver<3> > p_vegf_solver = FiniteDifferenceSolver<3>::Create();
        p_vegf_solver->SetExtents(p_domain, 40.0);
        p_vegf_solver->SetPde(p_vegf_pde);
        p_vegf_solver->SetVesselNetwork(p_network);

        return p_vegf_solver;
    }

public:

    void DOntTestCaBasedSpheroid() throw (Exception)
    {
        // Create the simulation domain
        boost::shared_ptr<Part<3> > p_domain = GetSimulationDomain();

        // Create a lattice for the cell population
        double spacing = 40.0;
        unsigned num_x = unsigned(p_domain->GetBoundingBox()[1]/spacing) + 1;
        unsigned num_y = unsigned(p_domain->GetBoundingBox()[3]/spacing) + 1;
        unsigned num_z = unsigned(p_domain->GetBoundingBox()[5]/spacing) + 1;
        PottsMeshGenerator<3> generator(num_x, 0, 0, num_y, 0, 0, num_z, 0, 0);
        PottsMesh<3>* p_mesh = generator.GetMesh();
        p_mesh->Scale(spacing, spacing, spacing);

        // Create a tumour cells in a cylinder in the middle of the domain
        boost::shared_ptr<Part<3> > p_tumour_cell_region = GetInitialTumourCellRegion();
        std::vector<unsigned> location_indices = p_tumour_cell_region->GetContainingGridIndices(num_x, num_y, num_z, spacing);

        std::vector<CellPtr> cells;
        CellsGenerator<SimpleOxygenBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size());

        // Create cell population
        CaBasedCellPopulation<3> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddCellWriter<CellLabelWriter>();

        // Create the vessel network
        boost::shared_ptr<CaVascularNetwork<3> > p_network = GetVesselNetwork();

        // Create the oxygen pde solver
        boost::shared_ptr<FiniteDifferenceSolver<3> > p_oxygen_solver = GetOxygenSolver(p_domain, p_network);

        // Create the vegf pde solver
        boost::shared_ptr<FiniteDifferenceSolver<3> > p_vegf_solver = GetVegfSolver(p_domain, p_network);

        // Create the angiogenesis solver
        boost::shared_ptr<AbstractAngiogenesisSolver<3> > p_angiogenesis_solver = AbstractAngiogenesisSolver<3>::Create(p_network);
        p_angiogenesis_solver->AddPdeSolver(p_oxygen_solver);
        p_angiogenesis_solver->AddPdeSolver(p_vegf_solver);

        boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> > p_grow_direction_modifier = OffLatticePrwGrowthDirectionModifier<3>::Create();
        boost::shared_ptr<OffLatticeTipAttractionGrowthDirectionModifier<3> > p_grow_direction_modifier2 = OffLatticeTipAttractionGrowthDirectionModifier<3>::Create();
        p_grow_direction_modifier2->SetNetwork(p_network);
        boost::shared_ptr<OffLatticeRandomNormalSproutingRule<3> > p_sprouting_rule = OffLatticeRandomNormalSproutingRule<3>::Create();
        p_sprouting_rule->SetSproutingProbability(0.05);

        p_angiogenesis_solver->AddGrowthDirectionModifier(p_grow_direction_modifier);
        p_angiogenesis_solver->AddGrowthDirectionModifier(p_grow_direction_modifier2);
        p_angiogenesis_solver->SetSproutingRule(p_sprouting_rule);
        p_angiogenesis_solver->SetAnastamosisRadius(5.0);

        boost::shared_ptr<AngiogenesisModifier<3> > p_simulation_modifier = boost::shared_ptr<AngiogenesisModifier<3> >(new AngiogenesisModifier<3>);
        p_simulation_modifier->SetAngiogenesisSolver(p_angiogenesis_solver);
        OnLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestAngiogenesisSimulationModifier/CaBased");
        simulator.SetDt(1.0);
        simulator.SetEndTime(10.0);
        simulator.AddSimulationModifier(p_simulation_modifier);
        simulator.Solve();
    }

    void TestNodeBasedSpheroid() throw (Exception)
    {
        //



        // Create the domain
        boost::shared_ptr<Part<3> > p_domain = GetSimulationDomain();

        // Create nodes corresponding to cell positions
        double spacing = 40.0;
        unsigned num_x = unsigned(p_domain->GetBoundingBox()[1]/spacing) + 1;
        unsigned num_y = unsigned(p_domain->GetBoundingBox()[3]/spacing) + 1;
        unsigned num_z = unsigned(p_domain->GetBoundingBox()[5]/spacing) + 1;

        // Create a tumour cells in a cylinder in the middle of the domain
        boost::shared_ptr<Part<3> > p_tumour_cell_region = GetInitialTumourCellRegion();
        std::vector<unsigned> location_indices = p_tumour_cell_region->GetContainingGridIndices(num_x, num_y, num_z, spacing);

        std::vector<Node<3>*> nodes;
        for(unsigned idx=0; idx<location_indices.size(); idx++)
        {
            c_vector<double, 3> location = Grid::GetLocationOf1dIndex(location_indices[idx], num_x, num_y, spacing);
            nodes.push_back(new Node<3>(idx, location, false));
        }
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5 * spacing);
        std::vector<CellPtr> cells;
        CellsGenerator<SimpleOxygenBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(2.0 * spacing);
        cell_population.AddCellWriter<CellLabelWriter>();

        // Create the vessel network
        boost::shared_ptr<CaVascularNetwork<3> > p_network = GetVesselNetwork();

        // Create the oxygen pde solver
        boost::shared_ptr<FiniteDifferenceSolver<3> > p_oxygen_solver = GetOxygenSolver(p_domain, p_network);

        // Create the vegf pde solver
        boost::shared_ptr<FiniteDifferenceSolver<3> > p_vegf_solver = GetVegfSolver(p_domain, p_network);

        // Create the angiogenesis solver
        boost::shared_ptr<AbstractAngiogenesisSolver<3> > p_angiogenesis_solver = AbstractAngiogenesisSolver<3>::Create(p_network);
        p_angiogenesis_solver->AddPdeSolver(p_oxygen_solver);
        p_angiogenesis_solver->AddPdeSolver(p_vegf_solver);

        boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> > p_grow_direction_modifier = OffLatticePrwGrowthDirectionModifier<3>::Create();
        boost::shared_ptr<OffLatticeTipAttractionGrowthDirectionModifier<3> > p_grow_direction_modifier2 = OffLatticeTipAttractionGrowthDirectionModifier<3>::Create();
        boost::shared_ptr<OffLatticeSolutionDependentGrowthDirectionModifier<3> > p_grow_direction_modifier3 = OffLatticeSolutionDependentGrowthDirectionModifier<3>::Create();
        p_grow_direction_modifier2->SetNetwork(p_network);
        p_grow_direction_modifier3->SetSolver(p_vegf_solver);
        p_grow_direction_modifier3->SetStrength(0.2);
        boost::shared_ptr<OffLatticeRandomNormalSproutingRule<3> > p_sprouting_rule = OffLatticeRandomNormalSproutingRule<3>::Create();
        p_sprouting_rule->SetSproutingProbability(0.005);

        p_angiogenesis_solver->AddGrowthDirectionModifier(p_grow_direction_modifier);
        p_angiogenesis_solver->AddGrowthDirectionModifier(p_grow_direction_modifier2);
        p_angiogenesis_solver->AddGrowthDirectionModifier(p_grow_direction_modifier3);
        p_angiogenesis_solver->SetSproutingRule(p_sprouting_rule);
        p_angiogenesis_solver->SetAnastamosisRadius(5.0);

        boost::shared_ptr<AngiogenesisModifier<3> > p_simulation_modifier = boost::shared_ptr<AngiogenesisModifier<3> >(new AngiogenesisModifier<3>);
        p_simulation_modifier->SetAngiogenesisSolver(p_angiogenesis_solver);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestAngiogenesisSimulationModifier/NodeBased");
        simulator.SetDt(1.0);
        simulator.SetEndTime(100.0);
        simulator.AddSimulationModifier(p_simulation_modifier);

        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        simulator.AddForce(p_force);
        simulator.Solve();

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
};

#endif /* TESTSPHEROIDWITHANGIOGENESIS_HPP_ */