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

#ifndef TEST2DVASCULARTUMOURGROWTH_HPP_
#define TEST2DVASCULARTUMOURGROWTH_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "Owen2011OxygenBasedCellCycleModel.hpp"
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
#include "ApoptoticCellKiller.hpp"

#include "VasculatureGenerator.hpp"
#include "CaVascularNetwork.hpp"
#include "CaBasedCellPopulationWithVessels.hpp"
#include "CaBasedCellPopulationWithVesselsGenerator.hpp"
#include "AngiogenesisSolver.hpp"
#include "AngiogenesisSolverUsingCellPopulationWithVessels.hpp"
#include "FLowSolver.hpp"

#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

class Test2dVascularTumourGrowth : public AbstractCellBasedWithTimingsTestSuite
{
    boost::shared_ptr<HybridLinearEllipticPde<2> > GetOxygenPde()
                    {
        boost::shared_ptr<HybridLinearEllipticPde<2> > p_pde = HybridLinearEllipticPde<2>::Create();
        p_pde->SetIsotropicDiffusionConstant(8700000 / (400.0)); // assume cell width is 20 microns

        // Add a cell state specific discrete source for cells consuming oxygen
        boost::shared_ptr<CellStateDependentDiscreteSource<2> > p_cell_sink =
                CellStateDependentDiscreteSource<2>::Create();
        p_cell_sink->SetIsLinearInSolution(true);
        std::map<unsigned, double> state_specific_rates;

        MAKE_PTR(WildTypeCellMutationState, p_normal_cell_state);
        MAKE_PTR(CancerCellMutationState, p_cancer_state);
        MAKE_PTR(QuiescentCancerCellMutationState, p_quiescent_cancer_state);
        MAKE_PTR(ApoptoticCellProperty, p_apoptotic_property);
        MAKE_PTR(TipCellMutationState, p_tip_state);
        MAKE_PTR(StalkCellMutationState, p_stalk_state);

        state_specific_rates[p_apoptotic_property->GetColour()] = 0.0;
        state_specific_rates[p_cancer_state->GetColour()] = -7800;
        state_specific_rates[p_quiescent_cancer_state->GetColour()] = -7800;
        state_specific_rates[p_tip_state->GetColour()] = -5000;
        state_specific_rates[p_stalk_state->GetColour()] = -5000;
        state_specific_rates[p_normal_cell_state->GetColour()] = -5000;

        p_cell_sink->SetStateRateMap(state_specific_rates);
        p_pde->AddDiscreteSource(p_cell_sink);

        // todo this needs to be updated so that source strength is proportional to haematocrit level
        boost::shared_ptr<CellStateDependentDiscreteSource<2> > p_cell_source =
                CellStateDependentDiscreteSource<2>::Create();
        p_cell_source->SetIsLinearInSolution(false); // constant oxygen release rate
        std::map<unsigned, double> state_specific_rates2;

        double segment_radius = 0.5;
        double segment_length = 1.0;
        double permeability = 20000;
        double inter_vessel_O2_level = 20;
        state_specific_rates2[p_stalk_state->GetColour()] = 2*M_PI*segment_radius*segment_length*permeability*inter_vessel_O2_level/pow(segment_length,3.0);
        state_specific_rates2[p_tip_state->GetColour()] = 2*M_PI*segment_radius*segment_length*permeability*inter_vessel_O2_level/pow(segment_length,3.0);
        p_cell_source->SetStateRateMap(state_specific_rates2);
        p_pde->AddDiscreteSource(p_cell_source);

        p_pde->SetVariableName("Oxygen");

        return p_pde;
                    }

    // todo need to check parameters in sink/source terms in here
    boost::shared_ptr<HybridLinearEllipticPde<2> > GetVegfPde()
                    {
        boost::shared_ptr<HybridLinearEllipticPde<2> > p_pde = HybridLinearEllipticPde<2>::Create();
        p_pde->SetIsotropicDiffusionConstant(60000 / (400.0)); // assume cell width is 20 microns
        p_pde->SetContinuumLinearInUTerm(-0.8); //Vegf decay

        // VEGF release for normal cells and quiescent cancer cells:
        // normal cells release only when intracellular vegf reaches a certain value
        // quiescent cancer cells release always release vegf
        boost::shared_ptr<CellStateDependentDiscreteSource<2> > p_normal_and_quiescent_cell_source = CellStateDependentDiscreteSource<
                2>::Create();
        p_normal_and_quiescent_cell_source->SetIsLinearInSolution(false); // constant vegf release rate

        // Set mutation specific source strengths and thresholds
        std::map<unsigned, double> normal_and_quiescent_cell_rates;
        std::map<unsigned, double> normal_and_quiescent_cell_rate_thresholds;

        MAKE_PTR(QuiescentCancerCellMutationState, p_quiescent_cancer_state);
        MAKE_PTR(WildTypeCellMutationState, p_normal_cell_state);
        normal_and_quiescent_cell_rates[p_normal_cell_state->GetColour()] = 0.6;
        normal_and_quiescent_cell_rate_thresholds[p_normal_cell_state->GetColour()] = 0.27;
        normal_and_quiescent_cell_rates[p_quiescent_cancer_state->GetColour()] = 0.6;
        normal_and_quiescent_cell_rate_thresholds[p_quiescent_cancer_state->GetColour()] = 0;

        p_normal_and_quiescent_cell_source->SetStateRateMap(normal_and_quiescent_cell_rates);
        p_normal_and_quiescent_cell_source->SetLabelName("VEGF");
        p_normal_and_quiescent_cell_source->SetStateRateThresholdMap(normal_and_quiescent_cell_rate_thresholds);
        p_pde->AddDiscreteSource(p_normal_and_quiescent_cell_source);

        // VEGF release for cancer cells and tip cells, now there is no threshold so we use a different source
        boost::shared_ptr<CellStateDependentDiscreteSource<2> > p_other_cell_sinks = CellStateDependentDiscreteSource<
                2>::Create();
        p_other_cell_sinks->SetIsLinearInSolution(true); // linear vegf uptake rate
        std::map<unsigned, double> other_cell_rates;

        double segment_radius = 0.5;
        double segment_length = 1.0;
        double permeability = 15;
        MAKE_PTR(TipCellMutationState, p_tip_state);
        MAKE_PTR(StalkCellMutationState, p_stalk_state);
        other_cell_rates[p_tip_state->GetColour()] = -2*M_PI*segment_radius*segment_length*permeability/pow(segment_length,3.0); // tip cell mutation state
        other_cell_rates[p_stalk_state->GetColour()] = -2*M_PI*segment_radius*segment_length*permeability/pow(segment_length,3.0); // stalk cell mutation state
        p_other_cell_sinks->SetStateRateMap(other_cell_rates);
        p_pde->AddDiscreteSource(p_other_cell_sinks);

        p_pde->SetVariableName("VEGF");

        return p_pde;
                    }

    boost::shared_ptr<CaVascularNetwork<2> > GetHexagonalNetwork(double domain_x, double domain_y)
                {

        VasculatureGenerator<2> network_generator;
        boost::shared_ptr<CaVascularNetwork<2> > p_network = network_generator.GenerateHexagonalNetwork(domain_x, domain_y, 7);

        std::vector<ChastePoint<2> > points;
        points.push_back(ChastePoint<2>(0, 0));
        points.push_back(ChastePoint<2>(5, 0));

        std::vector<boost::shared_ptr<VascularNode<2> > > nodes;
        for(unsigned i=0; i < points.size(); i++)
        {
            nodes.push_back(boost::shared_ptr<VascularNode<2> > (VascularNode<2>::Create(points[i])));
        }

        boost::shared_ptr<CaVesselSegment<2> > p_segment(CaVesselSegment<2>::Create(nodes[0], nodes[1]));

        double initial_vessel_radius = 10.0e-6;
        p_segment->SetRadius(initial_vessel_radius);
        double haematocrit = 0.45;
        p_segment->GetFlowProperties()->SetHaematocrit(haematocrit);
        double viscosity = 2e-3;
        p_segment->GetFlowProperties()->SetViscosity(viscosity);
        p_network->SetSegmentProperties(p_segment);

        std::vector<std::pair<double, double> > network_extents = p_network->GetExtents();
        double y_middle = (network_extents[1].first + network_extents[1].second) /2.0;
        double x_middle = (network_extents[0].first + network_extents[0].second) /2.0;

        std::vector<boost::shared_ptr<CaVessel<2> > >::iterator vessel_iterator;

        std::vector<boost::shared_ptr<CaVessel<2> > > vessels = p_network->GetVessels();

        for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
        {
            if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetStartNode()->GetLocation()[1] >  y_middle)
                {
                    if((*vessel_iterator)->GetStartNode()->GetLocation()[0] >  x_middle)
                    {
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3320);
                    }
                }
            }
            if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetEndNode()->GetLocation()[1] >  y_middle)
                {
                    if((*vessel_iterator)->GetStartNode()->GetLocation()[0] >  x_middle)
                    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3320);
                    }
                }
            }
            if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetStartNode()->GetLocation()[1] <=  y_middle)
                {
                    if((*vessel_iterator)->GetStartNode()->GetLocation()[0] <  x_middle)
                    {
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2090);
                    }
                }
            }
            if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetEndNode()->GetLocation()[1] <=  y_middle)
                {
                    if((*vessel_iterator)->GetStartNode()->GetLocation()[0] <  x_middle)
                    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2090);
                    }
                }
            }

        }

        return p_network;
                }

public:

    void TestOnLattice2dVascularTumourGrowth() throw (Exception)
    {
        // todo the simulation clearly goes wrong when domain is not square ... likely to be a problem with how interpolation is being done and indices being mixed up
        // Set up simulation domain
        double domain_x = 42.0;
        double domain_y = 42.0;
        boost::shared_ptr<Part<2> > p_domain = Part<2>::Create();
        p_domain->AddRectangle(domain_x, domain_y);

        // Create the initial vessel network: hexagonally tesselated vascular network
        boost::shared_ptr<CaVascularNetwork<2> > p_network =  GetHexagonalNetwork(domain_x,domain_y);

        // Create a lattice for the cell population
        double spacing = 1.0;
        unsigned num_x = unsigned(p_domain->GetBoundingBox()[1] / spacing);
        unsigned num_y = unsigned(p_domain->GetBoundingBox()[3] / spacing);
        PottsMeshGenerator<2> generator(num_x, 0, 0, num_y, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();
        p_mesh->Scale(spacing, spacing);

        // Create cell population
        CaBasedCellPopulationWithVesselsGenerator<2> cellPopulationGenerator;
        cellPopulationGenerator.SetIncludeNormalCellPopulation(true);
        boost::shared_ptr<CaBasedCellPopulationWithVessels<2> > cell_population =
                cellPopulationGenerator.CreateCellPopulation(*p_mesh, p_network);

        // Get initial tumour cell region
        double radius = 5.0;
        c_vector<double, 2> origin;
        origin[0] = round(domain_x/2);
        origin[1] = round(domain_y/2);
        boost::shared_ptr<Part<2> > p_sub_domain = Part<2>::Create();
        boost::shared_ptr<Polygon> circle = p_sub_domain->AddCircle(radius, origin);
        std::vector<unsigned> location_indices;
        for (unsigned ind = 0; ind < p_mesh->GetNumNodes(); ind++)
        {
            if (p_sub_domain->IsPointInPart(p_mesh->GetNode(ind)->rGetLocation()))
            {
                location_indices.push_back(ind);
            }
        }

        // mutate all cells inside circle to make them cancerous
        MAKE_PTR(CancerCellMutationState, p_cancerous_state);

        double oxygen_concentration = 30.0;
        for (unsigned i = 0; i < location_indices.size(); i++)
        {
            std::set<CellPtr> cells = cell_population-> GetCellsUsingLocationIndex(location_indices[i]);
            std::set<CellPtr>::iterator it;
            for (it = cells.begin(); it != cells.end(); it++)
            {
                if ((*it)->GetMutationState()->IsType<WildTypeCellMutationState>())
                {
                    (*it)->SetMutationState(p_cancerous_state);
                }
            }
        }

        // update all cells in population to prescribe an initial oxygen concentration and apoptosis time
        std::list<CellPtr> cells = cell_population->rGetCells();
        std::list<CellPtr>::iterator it;
        for (it = cells.begin(); it != cells.end(); ++it)
        {
            (*it)->GetCellData()->SetItem("oxygen", oxygen_concentration);
            (*it)->GetCellData()->SetItem("VEGF", 0.0);
            (*it)->GetCellData()->SetItem("p53", 0.0);
            (*it)->GetCellData()->SetItem("Number_of_cancerous_neighbours", 0.0);
            (*it)->GetCellData()->SetItem("Number_of_normal_neighbours", 0.0);
            (*it)->SetApoptosisTime(3);
        }

        // inform cell_population what to write out to files
        cell_population->SetOutputResultsForChasteVisualizer(false);
        cell_population->AddCellWriter<CellLabelWriter>();
        cell_population->AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        cell_population->AddCellWriter<CellMutationStatesWriter>();
        cell_population->AddCellWriter<CellProliferativeTypesWriter>();
        cell_population->AddCellWriter<CellProliferativePhasesWriter>();

        // set volume fractions occupied by each cell
        boost::shared_ptr<WildTypeCellMutationState> wild_mutation_state(new WildTypeCellMutationState);
        boost::shared_ptr<CancerCellMutationState> cancer_mutation_state(new CancerCellMutationState);
        boost::shared_ptr<QuiescentCancerCellMutationState> quiescent_cancer_mutation_state(new QuiescentCancerCellMutationState);
        boost::shared_ptr<StalkCellMutationState> stalk_mutation_state(new StalkCellMutationState);
        boost::shared_ptr<TipCellMutationState> tip_mutation_state(new TipCellMutationState);
        cell_population->SetVolumeFraction(wild_mutation_state,0.6);
        cell_population->SetVolumeFraction(cancer_mutation_state,0.6);
        cell_population->SetVolumeFraction(quiescent_cancer_mutation_state,0.6);
        cell_population->SetVolumeFraction(stalk_mutation_state,0.4);
        cell_population->SetVolumeFraction(tip_mutation_state,0.4);

        // Create a grid to solve PDEs on
        boost::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
        p_grid->SetSpacing(1.0);
        std::vector<unsigned> extents;
        extents.push_back(domain_x+1); // num_x
        extents.push_back(domain_y+1); // num_y
        extents.push_back(1); // num_z
        p_grid->SetExtents(extents);

        // Create the oxygen pde solver
        boost::shared_ptr<FiniteDifferenceSolver<2> > p_oxygen_solver = FiniteDifferenceSolver<2>::Create();
        p_oxygen_solver->SetGrid(p_grid);
        p_oxygen_solver->SetPde(GetOxygenPde());
        p_oxygen_solver->SetLabel("oxygen");

        // Create the vegf pde solver
        boost::shared_ptr<FiniteDifferenceSolver<2> > p_vegf_solver = FiniteDifferenceSolver<2>::Create();
        p_vegf_solver->SetGrid(p_grid);
        p_vegf_solver->SetPde(GetVegfPde());
        p_vegf_solver->SetLabel("VEGF");
        cell_population->AddPdeHandler(p_vegf_solver);

        // Create the vascular tumour solver, which manages all pde solves
        boost::shared_ptr<VascularTumourSolver<2> > p_vascular_tumour_solver = VascularTumourSolver<2>::Create();
        p_vascular_tumour_solver->SetVesselNetwork(p_network);
        p_vascular_tumour_solver->AddHybridSolver(p_oxygen_solver);
        p_vascular_tumour_solver->AddHybridSolver(p_vegf_solver);
        p_vascular_tumour_solver->SetOutputFrequency(1);
        // add angiogenesis solver to vascular tumour solver
        boost::shared_ptr<AngiogenesisSolverUsingCellPopulationWithVessels<2> > p_angiogenesis_solver(new AngiogenesisSolverUsingCellPopulationWithVessels<2>(cell_population));
        p_vascular_tumour_solver->SetAngiogenesisSolver(p_angiogenesis_solver);
        // todo currently there is an issue with the flow calculation - some blunt-ended vessels contain flow (they shouldn't)
        boost::shared_ptr<FlowSolver<2> > flow_solver(new FlowSolver<2>());
        p_vascular_tumour_solver->SetFlowSolver(flow_solver);

        OnLatticeSimulation<2> simulator(*(cell_population.get()));

        // Create the vascular tumour modifier which integrates with cell based Chaste
        boost::shared_ptr<VascularTumourModifier<2> > p_vascular_tumour_modifier = VascularTumourModifier<2>::Create();
        p_vascular_tumour_modifier->SetVascularTumourSolver(p_vascular_tumour_solver);

        simulator.AddSimulationModifier(p_vascular_tumour_modifier);

        // Create a Cell Concentration tracking modifier and add it to the simulation
        MAKE_PTR(Owen2011TrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        /*
         * Create cell killer to remove apoptotic cell from simulation
         */
        boost::shared_ptr<ApoptoticCellKiller<2> > apoptotic_cell_killer(new ApoptoticCellKiller<2>(cell_population.get()));
        simulator.AddCellKiller(apoptotic_cell_killer);

        std::string resultsDirectoryName = "Test2dVascularTumourGrowth/OnLattice";
        simulator.SetOutputDirectory(resultsDirectoryName);
        simulator.SetSamplingTimestepMultiple(1);
        // todo this seems to break simulations if dt is set to 1 - causes a CVode error:
        //          *CVODE Error -27 in module CVODE function CVode: tout too close to t0 to start integration.
        //          CVODE failed to solve system: CV_TOO_CLOSE
        simulator.SetDt(1.5);
        simulator.SetEndTime(200);

        simulator.Solve();

    }
};

#endif /*TESTOWEN2011TUMOURSPHEROIDSIMULATIONS_HPP_*/