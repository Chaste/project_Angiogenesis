//
//  TestPerfahl2011Model.hpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 03/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#ifndef TESTCABASEDCELLPOPULATIONWITHVESSELS_HPP
#define TESTCABASEDCELLPOPULATIONWITHVESSELS_HPP

#include <cxxtest/TestSuite.h>
#include <boost/lexical_cast.hpp>
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VasculatureGenerator.hpp"
#include "VascularNetwork.hpp"
#include "CaBasedCellPopulation.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VesselNetworkCellPopulationInteractor.hpp"
#include "UniformlyDistributedCellCycleModel.hpp"
#include "Owen2011OxygenBasedCellCycleModel.hpp"
#include "NodeLocationWriter.hpp"
#include "CellLabelWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellsGenerator.hpp"
#include "StalkCellMutationState.hpp"
#include "TipCellMutationState.hpp"
#include "SimulationTime.hpp"
#include "PottsMesh.hpp"
#include "PottsMeshGenerator.hpp"
#include "CaPopulationMigrationRule.hpp"
#include "AngiogenesisSolver.hpp"

#include "FakePetscSetup.hpp"

class TestAngiogenesisWithCaPopulation : public AbstractCellBasedWithTimingsTestSuite
{

public:

    void TestGrowSingleVessel() throw (Exception)
    {
        // Set up the cell population
        PottsMeshGenerator<2> generator(20, 0, 0, 20, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create the vessel network: single vessel in middle of domain
        c_vector<double, 2> start_position;
        start_position[0] = 10;
        start_position[1] = 0;
        VasculatureGenerator<2> network_generator;
        boost::shared_ptr<VascularNetwork<2> > p_network = network_generator.GenerateSingleVessel(10, start_position);

        // Write the initial network to file
        std::string output_directory = "TestAngiogenesisWithCaPopulation";
        MAKE_PTR_ARGS(OutputFileHandler, p_file_handler, (output_directory, false));
        std::string output_filename = p_file_handler->GetOutputDirectoryFullPath().append("InitialVesselNetwork.vtp");

        c_vector<double, 2> tip_position;
        tip_position[0] = 10;
        tip_position[1] = 10;
        p_network->GetNearestNode(tip_position)->SetIsMigrating(true);
        p_network->Write(output_filename);

        // Create endothelial cell population. First fill the domain, then kill all non-vessel cells, then label vessel cells
        std::vector<unsigned> location_indices;
        for (unsigned index = 0; index < p_mesh->GetNumNodes(); index++)
        {
            location_indices.push_back(index);
        }

        std::vector<CellPtr> cells;
        MAKE_PTR(DefaultCellProliferativeType, p_diff_type);
        MAKE_PTR(StalkCellMutationState, p_EC_state);
        MAKE_PTR(TipCellMutationState, p_EC_Tip_state);
        CellsGenerator<Owen2011OxygenBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_diff_type);
        MAKE_PTR_ARGS(CaBasedCellPopulation<2>, p_cell_population, (*p_mesh, cells, location_indices));

        VesselNetworkCellPopulationInteractor<2> interactor = VesselNetworkCellPopulationInteractor<2>();
        interactor.SetVesselNetwork(p_network);
        interactor.PartitionNetworkOverCells(*p_cell_population);
        interactor.KillNonVesselOverlappingCells(*p_cell_population);
        interactor.LabelVesselsInCellPopulation(*p_cell_population, p_EC_Tip_state, p_EC_state);

        std::string output_filename2 = p_file_handler->GetOutputDirectoryFullPath().append("AssociatedVesselNetwork.vtp");
        p_network->Write(output_filename2);
        p_cell_population->AddCellWriter<CellLabelWriter>();
        p_cell_population->AddCellWriter<CellMutationStatesWriter>();
        p_cell_population->AddPopulationWriter<NodeLocationWriter>();

        boost::shared_ptr<CaPopulationMigrationRule<2> > p_migration_rule = CaPopulationMigrationRule<2>::Create();
        p_migration_rule->SetCellPopulation(p_cell_population);

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10.0, 10);
        AngiogenesisSolver<2> angiogenesis_solver;
        angiogenesis_solver.SetMigrationRule(p_migration_rule);
        angiogenesis_solver.SetCellPopulation(p_cell_population);
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetOutputFileHandler(p_file_handler);
        angiogenesis_solver.Run(true);
    }

//    void TestVolumeFractionMap() throw (Exception)
//    {
//        // Create the mesh
//        PottsMeshGenerator<3> generator(20, 0, 0, 20, 0, 0, 21, 0, 0);
//        PottsMesh<3>* p_mesh = generator.GetMesh();
//
//        // Create the vessel network: single vessel in middle of domain
//        c_vector<double, 3> start_position;
//        start_position[0] = 10;
//        start_position[1] = 10;
//        start_position[2] = 0;
//        VasculatureGenerator<3> network_generator;
//        boost::shared_ptr<VascularNetwork<3> > p_network = network_generator.GenerateSingleVessel(20, start_position);
//
//        // Write the initial network to file
//        std::string output_directory = "TestCaBasedCellPopulationWithVesselsSelectTipCell";
//        OutputFileHandler output_file_handler(output_directory, false);
//        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append(
//                "InitialVesselNetwork.vtp");
//        p_network->Write(output_filename);
//
//        // create endothelial cell population
//        std::vector<unsigned> location_indices;
//        for (unsigned index = 0; index < p_mesh->GetNumNodes(); index++)
//        {
//            location_indices.push_back(index);
//        }
//
//        std::vector<CellPtr> cells;
//        MAKE_PTR(DefaultCellProliferativeType, p_diff_type);
//        MAKE_PTR(StalkCellMutationState, p_EC_state);
//        CellsGenerator<Owen2011OxygenBasedCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_diff_type);
//        CaBasedCellPopulation<3> cell_population(*p_mesh, cells, location_indices);
//
//        VesselNetworkCellPopulationInteractor<3> interactor = VesselNetworkCellPopulationInteractor<3>();
//        interactor.SetVesselNetwork(p_network);
//        interactor.PartitionNetworkOverCells(cell_population);
//        interactor.LabelVesselsInCellPopulation(cell_population, p_EC_state, p_EC_state);
//
//        boost::shared_ptr<WildTypeCellMutationState> wild_mutation_state(new WildTypeCellMutationState);
//        boost::shared_ptr<WildTypeCellMutationState> wild_mutation_state2(new WildTypeCellMutationState);
//        boost::shared_ptr<CancerCellMutationState> cancer_mutation_state(new CancerCellMutationState);
//        boost::shared_ptr<CancerCellMutationState> cancer_mutation_state2(new CancerCellMutationState);
//        boost::shared_ptr<MacrophageMutationState> mac_mutation_state(new MacrophageMutationState);
//        boost::shared_ptr<StalkCellMutationState> stalk_mutation_state(new StalkCellMutationState);
//        boost::shared_ptr<StalkCellMutationState> stalk_mutation_state2(new StalkCellMutationState);
//
//        cell_population.SetVolumeFraction(wild_mutation_state, 0.6);
//        cell_population.SetVolumeFraction(cancer_mutation_state2, 0.6);
//        cell_population.SetVolumeFraction(cancer_mutation_state, 0.8);
//        TS_ASSERT_THROWS_ANYTHING(cell_population.SetVolumeFraction(mac_mutation_state, 1.6));
//
//        cell_population.SetVolumeFraction(mac_mutation_state, 0.5);
//        cell_population.SetVolumeFraction(stalk_mutation_state, 0.9);
//
//        TS_ASSERT_DELTA(cell_population.GetOccupyingVolumeFraction(stalk_mutation_state2), 0.9, 1e-3);
//        TS_ASSERT_DELTA(cell_population.GetOccupyingVolumeFraction(cancer_mutation_state), 0.8, 1e-3);
//        TS_ASSERT_DELTA(cell_population.GetOccupyingVolumeFraction(mac_mutation_state), 0.5, 1e-3);
//        TS_ASSERT_DELTA(cell_population.GetOccupyingVolumeFraction(wild_mutation_state2), 0.6, 1e-3);
//
//        TS_ASSERT_EQUALS(cell_population.IsSiteAvailable(1, cell_population->GetCellUsingLocationIndex(1)), false);
//    }

};

#endif /*TESTCABASEDCELLPOPULATIONWITHVESSELS_HPP*/
