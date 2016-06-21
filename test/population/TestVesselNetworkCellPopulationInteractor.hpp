//
//  TestPerfahl2011Model.hpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 03/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#ifndef TestVesselNetworkCellPopulationInteractor_HPP
#define TestVesselNetworkCellPopulationInteractor_HPP

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
#include "Owen2011OxygenBasedCellCycleModel.hpp"
#include "CellLabelWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "NodeLocationWriter.hpp"
#include "CellsGenerator.hpp"
#include "StalkCellMutationState.hpp"
#include "TipCellMutationState.hpp"
#include "CancerCellMutationState.hpp"
#include "MacrophageMutationState.hpp"
#include "SimulationTime.hpp"
#include "PottsMesh.hpp"
#include "PottsMeshGenerator.hpp"
#include "FakePetscSetup.hpp"
#include "CellwiseSourcePde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "PdeAndBoundaryConditions.hpp"
#include "CellBasedPdeHandler.hpp"
#include "AveragedSourcePde.hpp"
#include "VesselNetworkCellPopulationInteractor.hpp"

class TestVesselNetworkCellPopulationInteractor : public AbstractCellBasedWithTimingsTestSuite
{

public:

    void TestSetUpLatticeBasedVessel() throw (Exception)
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
        boost::shared_ptr<VascularNetwork<3> > p_network = network_generator.GenerateSingleVessel(20, start_position);

        // Write the initial network to file
        std::string output_directory = "TestVesselNetworkCellPopulationInteractor";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("InitialVesselNetwork.vtp");
        p_network->Write(output_filename);

        // Create cell population
        // use OnLatticeVascularTumourCellPopulationGenerator to generate cells and associate cells
        // with vessels

        // create endothelial cell population
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index < p_mesh->GetNumNodes(); index++)
        {
            location_indices.push_back(index);
        }

        std::vector<CellPtr> cells;
        MAKE_PTR(DefaultCellProliferativeType, p_diff_type);
        MAKE_PTR(StalkCellMutationState, p_EC_state);
        CellsGenerator<Owen2011OxygenBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_diff_type);
        CaBasedCellPopulation<3> cell_population(*p_mesh, cells, location_indices);

        VesselNetworkCellPopulationInteractor<3> interactor = VesselNetworkCellPopulationInteractor<3>();
        interactor.SetVesselNetwork(p_network);
        interactor.PartitionNetworkOverCells(cell_population);
        interactor.LabelVesselsInCellPopulation(cell_population, p_EC_state, p_EC_state);

        TS_ASSERT_EQUALS(p_network->GetNumberOfNodes(), 21u);
        TS_ASSERT_EQUALS(p_network->GetNumberOfVessels(), 1u);
        TS_ASSERT_EQUALS(p_network->GetNumberOfVesselNodes(), 2u);
        // Test that each node has a cell associated with it
        for (unsigned idx = 0; idx < p_network->GetNumberOfNodes(); idx++)
        {
            if (!p_network->GetNode(idx)->HasCell())
            {
                std::cout << p_network->GetNode(idx)->GetLocationVector() << idx << std::endl;
            }
            TS_ASSERT(p_network->GetNode(idx)->GetNumberOfSegments() > 0);
        }

        std::string output_filename2 = output_file_handler.GetOutputDirectoryFullPath().append("AssociatedVesselNetwork.vtp");
        p_network->Write(output_filename2);

        // VTK writing needs a simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
        cell_population.AddCellWriter<CellLabelWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddPopulationWriter<NodeLocationWriter>();
        cell_population.OpenWritersFiles(output_file_handler);
        cell_population.WriteResultsToFiles(output_directory);
        cell_population.CloseWritersFiles();
    }
};

#endif /*TestVesselNetworkCellPopulationInteractor_HPP*/
