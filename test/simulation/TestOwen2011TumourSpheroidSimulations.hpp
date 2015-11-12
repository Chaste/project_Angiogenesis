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

#ifndef TESTOwen2011TUMOURSPHEROIDSIMULATIONS_HPP_
#define TESTOwen2011TUMOURSPHEROIDSIMULATIONS_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "RandomNumberGenerator.hpp"
#include "SmartPointers.hpp"
#include "ApoptoticCellKiller.hpp"
#include "AveragedSourcePde.hpp"
#include "Owen2011OxygenBasedCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CancerCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellLabelWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellwiseSourcePde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "CellBasedPdeHandler.hpp"
#include "OffLatticeSimulation.hpp"
#include "VoronoiDataWriter.hpp"
#include "Owen2011TrackingModifier.hpp"
#include "AngiogenesisModifier.hpp"
#include "FiniteDifferenceSolver.hpp"
#include "HybridLinearEllipticPde.hpp"
#include "DirichletBoundaryCondition.hpp"
#include "DiscreteSource.hpp"
#include "QuiescentCancerCellMutationState.hpp"
#include "CancerCellMutationState.hpp"
#include "ApoptoticCellProperty.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Debug.hpp"
#include "RegularGrid.hpp"
#include "CaBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "PottsMesh.hpp"
#include "OnLatticeSimulation.hpp"
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

class TestOwen2011TumourSpheroidSimulations : public AbstractCellBasedTestSuite
{
public:



    void TestSimpleOnLatticeSpheroidOwen2011OxygenBasedCellCycleModel() throw (Exception)
    {
        Timer::Reset();

        // set up simulation domain
        double domain_x = 40.0;
        double domain_y = 40.0;
        boost::shared_ptr<Part<2> > p_domain = Part<2>::Create();
        p_domain->AddRectangle(domain_x, domain_y);

        // Create a lattice for the cell population
        double spacing = 1.0;
        unsigned num_x = unsigned(p_domain->GetBoundingBox()[1]/spacing);
        unsigned num_y = unsigned(p_domain->GetBoundingBox()[3]/spacing);
        PottsMeshGenerator<2> generator(num_x, 0, 0, num_y, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();
        p_mesh->Scale(spacing, spacing);

        // get initital tumour cell region
        double radius = 3.0;
        c_vector<double, 2> origin;
        origin[0] = 18.0;
        origin[1] = 18.0;
        boost::shared_ptr<Part<2> > p_sub_domain = Part<2>::Create();
        boost::shared_ptr<Polygon> circle = p_sub_domain->AddCircle(radius, origin);
        std::vector<unsigned> location_indices;

        for(unsigned ind = 0; ind < p_mesh->GetNumNodes(); ind++)
        {
            if(p_sub_domain->IsPointInPart(p_mesh->GetNode(ind)->rGetLocation()))
            {
                location_indices.push_back(ind);
            }
        }

        // Make the cells and set up the state and type
        std::vector<CellPtr> cells;
        MAKE_PTR(CancerCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        // Set up oxygen_concentration (mmHg)
        double oxygen_concentration = 30.0;

        for (unsigned i=0; i<location_indices.size(); i++)
        {
            // Assign an oxygen based cell cycle model, which requires a dimension to be set.
            Owen2011OxygenBasedCellCycleModel* const p_model = new Owen2011OxygenBasedCellCycleModel;
            p_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);

            p_cell->SetApoptosisTime(30);
            cells.push_back(p_cell);

            // Start all cells with the specified oxygen concentration
            p_cell->GetCellData()->SetItem("oxygen", oxygen_concentration);
        }

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.SetOutputResultsForChasteVisualizer(false);
        cell_population.AddCellWriter<CellLabelWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();

        // Create a grid to solve PDEs on
        boost::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
        p_grid->SetSpacing(4.0);
        std::vector<unsigned> extents;
        extents.push_back(10); // num_x
        extents.push_back(10); // num_y
        extents.push_back(1); // num_z
        p_grid->SetExtents(extents);

        //        c_vector<double,2> origin; // Grid bottom left corner
        //        origin[0]= -20.0;
        //        origin[1]= -20.0;
        //        p_grid->SetOrigin(origin);

        // Create the oxygen pde, discrete sources and boundary condition
        boost::shared_ptr<HybridLinearEllipticPde<2> > p_oxygen_pde = HybridLinearEllipticPde<2>::Create();
        p_oxygen_pde->SetDiffusionConstant(0.0033/400.0); // assume cell width is 20 microns
        p_oxygen_pde->SetVariableName("oxygen");

        boost::shared_ptr<DiscreteSource<2> > p_cell_oxygen_sink = DiscreteSource<2>::Create();
        p_cell_oxygen_sink->SetType(SourceType::CELL); // cell population is added automatically in AngiogenesisModifier
        p_cell_oxygen_sink->SetSource(SourceStrength::PRESCRIBED);
        p_cell_oxygen_sink->SetValue(2.e-8);
        p_cell_oxygen_sink->SetIsLinearInSolution(true);
        p_oxygen_pde->AddDiscreteSource(p_cell_oxygen_sink);

        ApoptoticCellProperty apoptotic_property;
        QuiescentCancerCellMutationState quiescent_property;
        WildTypeCellMutationState normal_cell_state;
        std::map<unsigned, double > mutationSpecificConsumptionRateMap;
        mutationSpecificConsumptionRateMap[apoptotic_property.GetColour()] =  0.0;
        mutationSpecificConsumptionRateMap[p_state.get()->GetColour()] =  2.5e-8;
        mutationSpecificConsumptionRateMap[quiescent_property.GetColour()] =  2.5e-8;
        mutationSpecificConsumptionRateMap[normal_cell_state.GetColour()] =  1e-8;
        p_cell_oxygen_sink->SetMutationSpecificConsumptionRateMap(mutationSpecificConsumptionRateMap);

        boost::shared_ptr<DirichletBoundaryCondition<2> > p_domain_ox_boundary_condition = DirichletBoundaryCondition<2>::Create();
        p_domain_ox_boundary_condition->SetValue(oxygen_concentration);
        p_domain_ox_boundary_condition->SetType(BoundaryConditionType::OUTER);
        p_domain_ox_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

        // Create the pde solver
        boost::shared_ptr<FiniteDifferenceSolver<2> > p_oxygen_solver = FiniteDifferenceSolver<2>::Create();
        p_oxygen_solver->SetGrid(p_grid);
        p_oxygen_solver->SetPde(p_oxygen_pde);
        p_oxygen_solver->AddDirichletBoundaryCondition(p_domain_ox_boundary_condition);

        // Create the vegf pde, discrete sources and boundary condition
        boost::shared_ptr<HybridLinearEllipticPde<2> > p_vegf_pde = HybridLinearEllipticPde<2>::Create();
        p_vegf_pde->SetDiffusionConstant(0.0033/(400.0*145.0)); // assume cell width is 20 microns and vegf D is oxygen D/145.0
        p_vegf_pde->SetVariableName("VEGF");
        p_vegf_pde->SetLinearInUTerm(-0.8);

        // VEGF for normal cells
        boost::shared_ptr<DiscreteSource<2> > p_cell_vegf_source = DiscreteSource<2>::Create();
        p_cell_vegf_source->SetType(SourceType::CELL); // cell population is added automatically in AngiogenesisModifier
        p_cell_vegf_source->SetSource(SourceStrength::PRESCRIBED);
        p_cell_vegf_source->SetIsLinearInSolution(false);
        p_vegf_pde->AddDiscreteSource(p_cell_vegf_source);
        std::map<unsigned, double> vegf_cell_color_source_rates;
        std::map<unsigned, double > vegf_cell_color_source_thresholds;
        vegf_cell_color_source_rates[normal_cell_state.GetColour()] = 0.6; // Normal cell mutation state
        p_cell_vegf_source->SetMutationSpecificConsumptionRateMap(vegf_cell_color_source_rates);
        vegf_cell_color_source_thresholds[normal_cell_state.GetColour()] = 0.27;
        p_cell_vegf_source->SetLabelName("VEGF");
        p_cell_vegf_source->SetMutationSpecificConsumptionRateThresholdMap(vegf_cell_color_source_thresholds);

        // VEGF for cancer cells
        boost::shared_ptr<DiscreteSource<2> > p_cancer_cell_vegf_source = DiscreteSource<2>::Create();
        p_cancer_cell_vegf_source->SetType(SourceType::CELL); // cell population is added automatically in AngiogenesisModifier
        p_cancer_cell_vegf_source->SetSource(SourceStrength::PRESCRIBED);
        p_cancer_cell_vegf_source->SetIsLinearInSolution(false);
        std::map<unsigned, double> vegf_cancer_cell_color_source_rates;
        vegf_cancer_cell_color_source_rates[quiescent_property.GetColour()] = 0.6; // Quiescent cancer cell mutation state
        p_cell_vegf_source->SetMutationSpecificConsumptionRateMap(vegf_cancer_cell_color_source_rates);
        p_vegf_pde->AddDiscreteSource(p_cancer_cell_vegf_source);

        boost::shared_ptr<DirichletBoundaryCondition<2> > p_domain_vegf_boundary_condition = DirichletBoundaryCondition<2>::Create();
        p_domain_vegf_boundary_condition->SetValue(0.0);
        p_domain_vegf_boundary_condition->SetType(BoundaryConditionType::OUTER);
        p_domain_vegf_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

        // Create the pde solver
        boost::shared_ptr<FiniteDifferenceSolver<2> > p_vegf_solver = FiniteDifferenceSolver<2>::Create();
        p_vegf_solver->SetGrid(p_grid);
        p_vegf_solver->SetPde(p_vegf_pde);
        p_vegf_solver->AddDirichletBoundaryCondition(p_domain_vegf_boundary_condition);

        boost::shared_ptr<AngiogenesisSolver<2> > p_angiogenesis_solver = AngiogenesisSolver<2>::Create();
        p_angiogenesis_solver->AddPdeSolver(p_oxygen_solver);
        p_angiogenesis_solver->AddPdeSolver(p_vegf_solver);
        p_angiogenesis_solver->SetOutputFrequency(1);

        boost::shared_ptr<AngiogenesisModifier<2> > p_simulation_modifier = boost::shared_ptr<AngiogenesisModifier<2> >(new AngiogenesisModifier<2>);
        p_simulation_modifier->SetAngiogenesisSolver(p_angiogenesis_solver);

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.AddSimulationModifier(p_simulation_modifier);

        /*
         * We next set the output directory and end time.
         */
        std::string resultsDirectoryName = "TestOwen2011OnLatticeTumourSpheroidGrowthWithODEWithHybridSolver";
        simulator.SetOutputDirectory(resultsDirectoryName);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetDt(1);
        simulator.SetEndTime(300);

        /*
         * Create cell killer to remove apoptotic cell from simulation
         */
        //  boost::shared_ptr<ApoptoticCellKiller<2> > apoptotic_cell_killer(new ApoptoticCellKiller<2>(&cell_population));
        //  simulator.AddCellKiller(apoptotic_cell_killer);

        // Create a Cell Concentration tracking modifier and add it to the simulation
        MAKE_PTR(Owen2011TrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        simulator.Solve();

        Timer::Print("Elapsed time = ");
    }

    void dontTestSimpleOffLatticeSpheroidOwen2011OxygenBasedCellCycleModel() throw (Exception)
    {
        Timer::Reset();

        // Set up a mesh to define cells on
        HoneycombMeshGenerator generator(20, 20, 0);
        MutableMesh<2,2>* p_mesh = generator.GetCircularMesh(3);

        // Make the cells and set up the state and type
        std::vector<CellPtr> cells;
        MAKE_PTR(CancerCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        // Set up oxygen_concentration (mmHg)
        double oxygen_concentration = 30.0;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            // Assign an oxygen based cell cycle model, which requires a dimension to be set.
            Owen2011OxygenBasedCellCycleModel* const p_model = new Owen2011OxygenBasedCellCycleModel;
            p_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);

            p_cell->SetApoptosisTime(30);
            cells.push_back(p_cell);

            // Start all cells with the specified oxygen concentration
            p_cell->GetCellData()->SetItem("oxygen", oxygen_concentration);
        }

        // Create the cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputResultsForChasteVisualizer(false);
        cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();

        // Create a grid to solve PDEs on
        boost::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
        p_grid->SetSpacing(4.0);
        std::vector<unsigned> extents;
        extents.push_back(10); // num_x
        extents.push_back(10); // num_y
        extents.push_back(1); // num_z
        p_grid->SetExtents(extents);

        c_vector<double,2> origin; // Grid bottom left corner
        origin[0]= -18.0;
        origin[1]= -18.0;
        p_grid->SetOrigin(origin);

        // Create the oxygen pde, discrete sources and boundary condition
        boost::shared_ptr<HybridLinearEllipticPde<2> > p_oxygen_pde = HybridLinearEllipticPde<2>::Create();
        p_oxygen_pde->SetDiffusionConstant(0.0033/400.0); // assume cell width is 20 microns
        p_oxygen_pde->SetVariableName("oxygen");

        boost::shared_ptr<DiscreteSource<2> > p_cell_oxygen_sink = DiscreteSource<2>::Create();
        p_cell_oxygen_sink->SetType(SourceType::CELL); // cell population is added automatically in AngiogenesisModifier
        p_cell_oxygen_sink->SetSource(SourceStrength::PRESCRIBED);
        p_cell_oxygen_sink->SetValue(2.e-8);
        p_cell_oxygen_sink->SetIsLinearInSolution(true);
        p_oxygen_pde->AddDiscreteSource(p_cell_oxygen_sink);

        ApoptoticCellProperty apoptotic_property;
        QuiescentCancerCellMutationState quiescent_property;
        WildTypeCellMutationState normal_cell_state;
        std::map<unsigned, double > mutationSpecificConsumptionRateMap;
        mutationSpecificConsumptionRateMap[apoptotic_property.GetColour()] =  0.0;
        mutationSpecificConsumptionRateMap[p_state.get()->GetColour()] =  2.5e-8;
        mutationSpecificConsumptionRateMap[quiescent_property.GetColour()] =  2.5e-8;
        mutationSpecificConsumptionRateMap[normal_cell_state.GetColour()] =  1e-8;
        p_cell_oxygen_sink->SetMutationSpecificConsumptionRateMap(mutationSpecificConsumptionRateMap);

        boost::shared_ptr<DirichletBoundaryCondition<2> > p_domain_ox_boundary_condition = DirichletBoundaryCondition<2>::Create();
        p_domain_ox_boundary_condition->SetValue(oxygen_concentration);
        p_domain_ox_boundary_condition->SetType(BoundaryConditionType::OUTER);
        p_domain_ox_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

        // Create the pde solver
        boost::shared_ptr<FiniteDifferenceSolver<2> > p_oxygen_solver = FiniteDifferenceSolver<2>::Create();
        p_oxygen_solver->SetGrid(p_grid);
        p_oxygen_solver->SetPde(p_oxygen_pde);
        p_oxygen_solver->AddDirichletBoundaryCondition(p_domain_ox_boundary_condition);

        // Create the vegf pde, discrete sources and boundary condition
        boost::shared_ptr<HybridLinearEllipticPde<2> > p_vegf_pde = HybridLinearEllipticPde<2>::Create();
        p_vegf_pde->SetDiffusionConstant(0.0033/(400.0*145.0)); // assume cell width is 20 microns and vegf D is oxygen D/145.0
        p_vegf_pde->SetVariableName("VEGF");
        p_vegf_pde->SetLinearInUTerm(0.8);

        // VEGF for normal cells
        boost::shared_ptr<DiscreteSource<2> > p_cell_vegf_source = DiscreteSource<2>::Create();
        p_cell_vegf_source->SetType(SourceType::CELL); // cell population is added automatically in AngiogenesisModifier
        p_cell_vegf_source->SetSource(SourceStrength::PRESCRIBED);
        p_cell_vegf_source->SetIsLinearInSolution(false);
        p_vegf_pde->AddDiscreteSource(p_cell_vegf_source);
        std::map<unsigned, double> vegf_cell_color_source_rates;
        std::map<unsigned, double > vegf_cell_color_source_thresholds;
        vegf_cell_color_source_rates[normal_cell_state.GetColour()] = -0.6; // Normal cell mutation state
        p_cell_vegf_source->SetMutationSpecificConsumptionRateMap(vegf_cell_color_source_rates);
        vegf_cell_color_source_thresholds[normal_cell_state.GetColour()] = 0.27;
        p_cell_vegf_source->SetLabelName("VEGF");
        p_cell_vegf_source->SetMutationSpecificConsumptionRateThresholdMap(vegf_cell_color_source_thresholds);

        // VEGF for cancer cells
        boost::shared_ptr<DiscreteSource<2> > p_cancer_cell_vegf_source = DiscreteSource<2>::Create();
        p_cancer_cell_vegf_source->SetType(SourceType::CELL); // cell population is added automatically in AngiogenesisModifier
        p_cancer_cell_vegf_source->SetSource(SourceStrength::PRESCRIBED);
        p_cancer_cell_vegf_source->SetIsLinearInSolution(false);
        std::map<unsigned, double> vegf_cancer_cell_color_source_rates;
        vegf_cancer_cell_color_source_rates[quiescent_property.GetColour()] = 0.6; // Quiescent cancer cell mutation state
        p_cell_vegf_source->SetMutationSpecificConsumptionRateMap(vegf_cancer_cell_color_source_rates);
        p_vegf_pde->AddDiscreteSource(p_cancer_cell_vegf_source);

        boost::shared_ptr<DirichletBoundaryCondition<2> > p_domain_vegf_boundary_condition = DirichletBoundaryCondition<2>::Create();
        p_domain_vegf_boundary_condition->SetValue(0.0);
        p_domain_vegf_boundary_condition->SetType(BoundaryConditionType::OUTER);
        p_domain_vegf_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

        // Create the pde solver
        boost::shared_ptr<FiniteDifferenceSolver<2> > p_vegf_solver = FiniteDifferenceSolver<2>::Create();
        p_vegf_solver->SetGrid(p_grid);
        p_vegf_solver->SetPde(p_vegf_pde);
        p_vegf_solver->AddDirichletBoundaryCondition(p_domain_vegf_boundary_condition);

        boost::shared_ptr<AngiogenesisSolver<2> > p_angiogenesis_solver = AngiogenesisSolver<2>::Create();
        p_angiogenesis_solver->AddPdeSolver(p_oxygen_solver);
        p_angiogenesis_solver->AddPdeSolver(p_vegf_solver);
        p_angiogenesis_solver->SetOutputFrequency(100);

        boost::shared_ptr<AngiogenesisModifier<2> > p_simulation_modifier = boost::shared_ptr<AngiogenesisModifier<2> >(new AngiogenesisModifier<2>);
        p_simulation_modifier->SetAngiogenesisSolver(p_angiogenesis_solver);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.AddSimulationModifier(p_simulation_modifier);

        /*
         * We next set the output directory and end time.
         */
        std::string resultsDirectoryName = "TestOwen2011OffLatticeTumourSpheroidGrowthWithODEWithHybridSolver";
        simulator.SetOutputDirectory(resultsDirectoryName);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetDt(0.01);
        simulator.SetEndTime(150);

        /*
         * Create cell killer to remove apoptotic cell from simulation
         */
        //  boost::shared_ptr<ApoptoticCellKiller<2> > apoptotic_cell_killer(new ApoptoticCellKiller<2>(&cell_population));
        //  simulator.AddCellKiller(apoptotic_cell_killer);


        // Create a Cell Concentration tracking modifier and add it to the simulation
        MAKE_PTR(Owen2011TrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(3);
        simulator.AddForce(p_linear_force);

        simulator.Solve();

        Timer::Print("Elapsed time = ");
    }


};

#endif /*TESTOwen2011TUMOURSPHEROIDSIMULATIONS_HPP_*/
