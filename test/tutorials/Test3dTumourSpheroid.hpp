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

#ifndef TEST3DTUMOURSPHEROIDSIMULATIONS_HPP_
#define TEST3DTUMOURSPHEROIDSIMULATIONS_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "Owen2011OxygenBasedCellCycleModel.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellLabelWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "OffLatticeSimulation.hpp"
#include "VoronoiDataWriter.hpp"
#include "Owen2011TrackingModifier.hpp"
#include "QuiescentCancerCellMutationState.hpp"
#include "CancerCellMutationState.hpp"
#include "ApoptoticCellProperty.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CancerCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "RegularGrid.hpp"
#include "FiniteDifferenceSolver.hpp"
#include "HybridLinearEllipticPde.hpp"
#include "HybridBoundaryCondition.hpp"
#include "CellStateDependentDiscreteSource.hpp"
#include "DiscreteSource.hpp"
#include "VascularTumourModifier.hpp"
#include "CaBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "PottsMesh.hpp"
#include "OnLatticeSimulation.hpp"

#include "PetscSetupAndFinalize.hpp"

class Test3dTumourSpheroid : public AbstractCellBasedWithTimingsTestSuite
{
	boost::shared_ptr<HybridLinearEllipticPde<3> > GetOxygenPde()
	{
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_pde = HybridLinearEllipticPde<3>::Create();
        p_pde->SetIsotropicDiffusionConstant(0.0033/400.0); // assume cell width is 20 microns

        // Add a cell state specific discrete source for cells consuming oxygen
        boost::shared_ptr<CellStateDependentDiscreteSource<3> > p_cell_sink = CellStateDependentDiscreteSource<3>::Create();
        std::map<unsigned, double > state_specific_rates;

		MAKE_PTR(WildTypeCellMutationState, p_normal_cell_state);
        MAKE_PTR(CancerCellMutationState, p_cancer_state);
        MAKE_PTR(QuiescentCancerCellMutationState, p_quiescent_cancer_state);
        MAKE_PTR(ApoptoticCellProperty, p_apoptotic_property);

        state_specific_rates[p_apoptotic_property->GetColour()] =  0.0;
        state_specific_rates[p_cancer_state->GetColour()] =  -3e-8;
        state_specific_rates[p_quiescent_cancer_state->GetColour()] =  -3e-8;
        state_specific_rates[p_normal_cell_state->GetColour()] =  -2e-8;

        p_cell_sink->SetStateRateMap(state_specific_rates);
        p_pde->AddDiscreteSource(p_cell_sink);

        return p_pde;
	}

	boost::shared_ptr<HybridLinearEllipticPde<3> > GetVegfPde()
	{
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_pde = HybridLinearEllipticPde<3>::Create();
        p_pde->SetIsotropicDiffusionConstant(0.0033/(400.0*145.0)); // assume cell width is 20 microns and vegf D is oxygen D/145.0
        p_pde->SetContinuumLinearInUTerm(0.8); //Vegf decay

        // VEGF release for normal cells, release only when intracellular vegf reaches a certain value
        boost::shared_ptr<CellStateDependentDiscreteSource<3> > p_normal_cell_source = CellStateDependentDiscreteSource<3>::Create();
        p_normal_cell_source->SetIsLinearInSolution(false); // constant vegf release rate

        // Set mutation specific source strengths and thresholds
        std::map<unsigned, double> normal_cell_rates;
        std::map<unsigned, double> normal_cell_rate_thresholds;

		MAKE_PTR(WildTypeCellMutationState, p_normal_cell_state);
		normal_cell_rates[p_normal_cell_state->GetColour()] = 0.6;
		normal_cell_rate_thresholds[p_normal_cell_state->GetColour()] = 0.27;

		p_normal_cell_source->SetStateRateMap(normal_cell_rates);
		p_normal_cell_source->SetLabelName("VEGF");
		p_normal_cell_source->SetStateRateThresholdMap(normal_cell_rate_thresholds);
		p_pde->AddDiscreteSource(p_normal_cell_source);

        // VEGF release for cancer cells, now there is no threshold so we use a different source
        boost::shared_ptr<CellStateDependentDiscreteSource<3> > p_cancer_cell_source = CellStateDependentDiscreteSource<3>::Create();
        p_cancer_cell_source->SetIsLinearInSolution(false); // constant vegf release rate
        std::map<unsigned, double> cancer_cell_rates;

        MAKE_PTR(QuiescentCancerCellMutationState, p_quiescent_cancer_state);
        cancer_cell_rates[p_quiescent_cancer_state->GetColour()] = 0.6; // Quiescent cancer cell mutation state
        p_cancer_cell_source->SetStateRateMap(cancer_cell_rates);
        p_pde->AddDiscreteSource(p_cancer_cell_source);

        return p_pde;
	}

public:

    void TestOnLattice3dTumourSpheroid() throw (Exception)
	{
	   // Set up simulation domain
	   double domain_x = 40.0;
	   double domain_y = 40.0;
	   double domain_z = 5.0;
	   boost::shared_ptr<Part<3> > p_domain = Part<3>::Create();
	   p_domain->AddCuboid(domain_x, domain_y, domain_z);

	   // Create a lattice for the cell population
	   double spacing = 1.0;
	   unsigned num_x = unsigned(p_domain->GetBoundingBox()[1]/spacing);
	   unsigned num_y = unsigned(p_domain->GetBoundingBox()[3]/spacing);
	   unsigned num_z = unsigned(p_domain->GetBoundingBox()[5]/spacing);
	   PottsMeshGenerator<3> generator(num_x, 0, 0, num_y, 0, 0,num_z ,0,0);
	   PottsMesh<3>* p_mesh = generator.GetMesh();
	   p_mesh->Scale(spacing, spacing);

	   // Get initital tumour cell region
	   double radius = 3.0;
	   c_vector<double, 3> origin;
	   origin[0] = 18.0;
	   origin[1] = 18.0;
	   origin[2] = 4.0;
	   boost::shared_ptr<Part<3> > p_sub_domain = Part<3>::Create();
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

	   double oxygen_concentration = 30.0;
	   for (unsigned i=0; i<location_indices.size(); i++)
	   {
		   // Assign an oxygen based cell cycle model, which requires a dimension to be set.
		   Owen2011OxygenBasedCellCycleModel* const p_model = new Owen2011OxygenBasedCellCycleModel;
		   p_model->SetDimension(3);
		   CellPtr p_cell(new Cell(p_state, p_model));
		   p_cell->SetCellProliferativeType(p_stem_type);
		   p_cell->SetApoptosisTime(30);
		   cells.push_back(p_cell);
		   p_cell->GetCellData()->SetItem("oxygen", oxygen_concentration);
	   }

	   // Create cell population
	   CaBasedCellPopulation<3> cell_population(*p_mesh, cells, location_indices);
	   cell_population.SetOutputResultsForChasteVisualizer(false);
	   cell_population.AddCellWriter<CellLabelWriter>();
	   cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
	   cell_population.AddCellWriter<CellMutationStatesWriter>();
	   cell_population.AddCellWriter<CellProliferativeTypesWriter>();
	   cell_population.AddCellWriter<CellProliferativePhasesWriter>();

	   // Create a grid to solve PDEs on
	   boost::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
	   p_grid->SetSpacing(2.0);
	   std::vector<unsigned> extents;
	   extents.push_back(40); // num_x
	   extents.push_back(40); // num_y
	   extents.push_back(5); // num_z
	   p_grid->SetExtents(extents);

       // Create the oxygen pde solver
       boost::shared_ptr<FiniteDifferenceSolver<3> > p_oxygen_solver = FiniteDifferenceSolver<3>::Create();
       p_oxygen_solver->SetGrid(p_grid);
       p_oxygen_solver->SetPde(GetOxygenPde());
       p_oxygen_solver->SetLabel("oxygen");

       // Add a dirichlet boundary condition for oxygen on the outer walls of the domain
       boost::shared_ptr<HybridBoundaryCondition<3> > p_domain_ox_boundary_condition = HybridBoundaryCondition<3>::Create();
       p_domain_ox_boundary_condition->SetValue(oxygen_concentration);
       p_oxygen_solver->AddBoundaryCondition(p_domain_ox_boundary_condition);

       // Create the vegf pde solver
       boost::shared_ptr<FiniteDifferenceSolver<3> > p_vegf_solver = FiniteDifferenceSolver<3>::Create();
       p_vegf_solver->SetGrid(p_grid);
       p_vegf_solver->SetPde(GetVegfPde());

       // Add a 0.0 dirichlet boundary condition for vegf on the outer walls of the domain
       boost::shared_ptr<HybridBoundaryCondition<3> > p_domain_vegf_boundary_condition = HybridBoundaryCondition<3>::Create();
       p_vegf_solver->AddBoundaryCondition(p_domain_vegf_boundary_condition);
       p_vegf_solver->SetLabel("VEGF");

       // Create the vascular tumour solver, which manages all pde solves
       boost::shared_ptr<VascularTumourSolver<3> > p_vascular_tumour_solver = VascularTumourSolver<3>::Create();
       p_vascular_tumour_solver->AddHybridSolver(p_oxygen_solver);
       p_vascular_tumour_solver->AddHybridSolver(p_vegf_solver);
       p_vascular_tumour_solver->SetOutputFrequency(10);

	   OnLatticeSimulation<3> simulator(cell_population);

       // Create the vascular tumour modifier which integrates with cell based Chaste
       boost::shared_ptr<VascularTumourModifier<3> > p_vascular_tumour_modifier = VascularTumourModifier<3>::Create();
       p_vascular_tumour_modifier->SetVascularTumourSolver(p_vascular_tumour_solver);
       simulator.AddSimulationModifier(p_vascular_tumour_modifier);

	   // Create a Cell Concentration tracking modifier and add it to the simulation
	   MAKE_PTR(Owen2011TrackingModifier<3>, p_modifier);
	   simulator.AddSimulationModifier(p_modifier);

	   std::string resultsDirectoryName = "Test3dTumourSpheroid/OnLattice";
	   simulator.SetOutputDirectory(resultsDirectoryName);
	   simulator.SetSamplingTimestepMultiple(10);
	   simulator.SetDt(1);
	   simulator.SetEndTime(300);

	   simulator.Solve();
	}
};

#endif /*TEST3DTUMOURSPHEROIDSIMULATIONS_HPP_*/
