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
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VasculatureGenerator.hpp"
#include "CaVascularNetwork.hpp"
#include "AngiogenesisSolver.hpp"
#include "CaBasedCellPopulationWithVessels.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "CellLabelWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "NodeLocationWriter.hpp"
#include "CellsGenerator.hpp"
#include "StalkCellMutationState.hpp"
#include "TipCellMutationState.hpp"
#include "CancerCellMutationState.hpp"
#include "SimulationTime.hpp"
#include "PottsMesh.hpp"
#include "PottsMeshGenerator.hpp"
#include "FakePetscSetup.hpp"

#include "Debug.hpp"

class TestCaBasedCellPopulationWithVessels : public AbstractCellBasedWithTimingsTestSuite
{

public:

    void TestSetUpSingleVessel() throw(Exception)
    {
        // Create the mesh
        PottsMeshGenerator<3> generator(20, 0, 0, 20, 0, 0, 21, 0, 0);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        // Create the vessel network: single vessel in middle of domain
        c_vector<double, 3> start_position;
        start_position[0] = 10;
        start_position[1] = 10;
        start_position[2] = 0;
        VasculatureGenerator<3> network_generator;
        boost::shared_ptr<CaVascularNetwork<3> > p_network = network_generator.GenerateSingleVessel(20 ,start_position);

        // Write the initial network to file
        std::string output_directory = "TestCaBasedCellPopulationWithVessels";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("InitialVesselNetwork.vtp");
        p_network->Write(output_filename);

        // Create some cells to fill the domain
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_transit_type);
        CellsGenerator<StochasticDurationCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_transit_type);

        MAKE_PTR(CancerCellMutationState, p_cancer_cell_mutation_state);
        std::vector<unsigned> location_indices;
        for(unsigned idx = 0; idx< p_mesh->GetNumNodes(); idx++)
        {
            location_indices.push_back(idx);
            cells[idx]->SetMutationState(p_cancer_cell_mutation_state);
        }

        // Create a combined cell and vessel population
        MAKE_PTR_ARGS(CaBasedCellPopulationWithVessels<3>, p_combined_population, (*p_mesh, cells, location_indices));

        // Define mutation states for tips and stalks
        MAKE_PTR(StalkCellMutationState, p_stalk_cell_mutation_state);
        MAKE_PTR(TipCellMutationState, p_tip_cell_mutation_state);

        // Add the vessels
        p_combined_population->SetVesselNetwork(p_network);
        p_combined_population->AsscoiateVesselNetworkWithCells(p_stalk_cell_mutation_state,
                                                               p_tip_cell_mutation_state);

        // Test that the network is correctly associated with the cells
        // are 20 nodes in the network, 2 vessel nodes, and 1 vessel
        TS_ASSERT_EQUALS(p_network->GetNumberOfNodes(),21u);
        TS_ASSERT_EQUALS(p_network->GetNumberOfVessels(),1u);
        TS_ASSERT_EQUALS(p_network->GetNumberOfVesselNodes(),2u);
        // Test that each node has a cell associated with it
        for (unsigned idx = 0; idx < p_network->GetNumberOfNodes(); idx++)
        {
            if(!p_network->GetNode(idx)->HasCell())
            {
                std::cout << p_network->GetNode(idx)->GetLocationVector() << idx << std::endl;
            }
            TS_ASSERT(p_network->GetNode(idx)->HasCell());
        }

        std::string output_filename2 = output_file_handler.GetOutputDirectoryFullPath().append("AssociatedVesselNetwork.vtp");
        p_network->Write(output_filename2);

        // VTK writing needs a simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
        p_combined_population->AddCellWriter<CellLabelWriter>();
        p_combined_population->AddCellWriter<CellMutationStatesWriter>();
        p_combined_population->AddPopulationWriter<NodeLocationWriter>();
        p_combined_population->OpenWritersFiles(output_file_handler);
        p_combined_population->WriteResultsToFiles(output_directory);
        p_combined_population->CloseWritersFiles();
    }
};

#endif /*TESTCABASEDCELLPOPULATIONWITHVESSELS_HPP*/
