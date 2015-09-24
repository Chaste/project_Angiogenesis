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
#include "CaVascularNetwork.hpp"
#include "CaBasedCellPopulationWithVessels.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "OnLatticeVascularTumourCellPopulationGenerator.hpp"
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
#include "CellwiseSourcePde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "PdeAndBoundaryConditions.hpp"
#include "CellBasedPdeHandler.hpp"
#include "AveragedSourcePde.hpp"

#include "Debug.hpp"

class TestCaBasedCellPopulationWithVessels : public AbstractCellBasedWithTimingsTestSuite
{

public:

    void TestSetUpSingleVessel() throw (Exception)
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
        boost::shared_ptr<CaVascularNetwork<3> > p_network = network_generator.GenerateSingleVessel(20, start_position);

        // Write the initial network to file
        std::string output_directory = "TestCaBasedCellPopulationWithVessels";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append(
                "InitialVesselNetwork.vtp");
        p_network->Write(output_filename);

        // Create cell population
        // use OnLatticeVascularTumourCellPopulationGenerator to generate cells and associate cells
        // with vessels
        OnLatticeVascularTumourCellPopulationGenerator<3> cellPopulationGenerator;
        cellPopulationGenerator.SetIncludeNormalCellPopulation(false);
        boost::shared_ptr<CaBasedCellPopulationWithVessels<3> > cell_population =
                cellPopulationGenerator.CreateCellPopulation(*p_mesh, p_network);

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
            TS_ASSERT(p_network->GetNode(idx)->HasCell());
            TS_ASSERT(p_network->GetNode(idx)->GetNumberOfSegments() > 0);
        }

        std::string output_filename2 = output_file_handler.GetOutputDirectoryFullPath().append(
                "AssociatedVesselNetwork.vtp");
        p_network->Write(output_filename2);

        // VTK writing needs a simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
        cell_population->AddCellWriter<CellLabelWriter>();
        cell_population->AddCellWriter<CellMutationStatesWriter>();
        cell_population->AddPopulationWriter<NodeLocationWriter>();
        cell_population->OpenWritersFiles(output_file_handler);
        cell_population->WriteResultsToFiles(output_directory);
        cell_population->CloseWritersFiles();
    }

    void TestSelectTipCell() throw (Exception)
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
        boost::shared_ptr<CaVascularNetwork<3> > p_network = network_generator.GenerateSingleVessel(20, start_position);

        // Write the initial network to file
        std::string output_directory = "TestCaBasedCellPopulationWithVesselsSelectTipCell";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append(
                "InitialVesselNetwork.vtp");
        p_network->Write(output_filename);

        // Create cell population
        OnLatticeVascularTumourCellPopulationGenerator<3> cellPopulationGenerator;
        cellPopulationGenerator.SetIncludeNormalCellPopulation(false);
        boost::shared_ptr<CaBasedCellPopulationWithVessels<3> > cell_population =
                cellPopulationGenerator.CreateCellPopulation(*p_mesh, p_network);

        // Test that the network is correctly associated with the cells
        // are 20 nodes in the network, 2 vessel nodes, and 1 vessel
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
            TS_ASSERT(p_network->GetNode(idx)->HasCell());
            if (idx == 10)
            {
                cell_population->SelectTipCell(p_network->GetNode(idx)->GetCell());
            }
        }

        TS_ASSERT_EQUALS(cell_population->GetNumberOfTipCells(), 1u);
        TS_ASSERT_EQUALS(cell_population->GetTipCells()[0], p_network->GetNode(10)->GetCell());

        std::string output_filename2 = output_file_handler.GetOutputDirectoryFullPath().append(
                "AssociatedVesselNetwork.vtp");
        p_network->Write(output_filename2);

        // VTK writing needs a simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
        cell_population->AddCellWriter<CellLabelWriter>();
        cell_population->AddCellWriter<CellMutationStatesWriter>();
        cell_population->AddPopulationWriter<NodeLocationWriter>();
        cell_population->OpenWritersFiles(output_file_handler);
        cell_population->WriteResultsToFiles(output_directory);
        cell_population->CloseWritersFiles();

        cell_population->DeselectTipCell(p_network->GetNode(10)->GetCell());
        TS_ASSERT_EQUALS(cell_population->GetNumberOfTipCells(), 0u);

    }

    void TestAngiogenesis() throw (Exception)
    {
        // Create the mesh
        PottsMeshGenerator<3> generator(10, 0, 0, 50, 0, 0, 50, 0, 0);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        // Create the vessel network: single vessel in middle of domain
        c_vector<double, 3> start_position;
        start_position[0] = 5;
        start_position[1] = 10;
        start_position[2] = 0;
        VasculatureGenerator<3> network_generator;
        boost::shared_ptr<CaVascularNetwork<3> > p_network = network_generator.GenerateSingleVessel(20, start_position);

        // Write the initial network to file
        std::string output_directory = "TestCaBasedCellPopulationWithVesselsAngiogenesis";
        OutputFileHandler output_file_handler(output_directory, true);

        OnLatticeVascularTumourCellPopulationGenerator<3> cellPopulationGenerator;
        cellPopulationGenerator.SetIncludeNormalCellPopulation(false);
        boost::shared_ptr<CaBasedCellPopulationWithVessels<3> > cell_population =
                cellPopulationGenerator.CreateCellPopulation(*p_mesh, p_network);

        // add writers for outputting
        cell_population->AddCellWriter<CellLabelWriter>();
        cell_population->AddCellWriter<CellMutationStatesWriter>();
        cell_population->AddPopulationWriter<NodeLocationWriter>();

        // set up PDE
        AveragedSourcePde<3> pde(*(cell_population.get()), -0.30);
        ConstBoundaryCondition<3> p_zero_boundary_condition(1.0);
        PdeAndBoundaryConditions<3> pde_and_bc(&pde, &p_zero_boundary_condition, false);
        pde_and_bc.SetDependentVariableName("VEGF");

        boost::shared_ptr<CellBasedPdeHandler<3> > pde_handler(new CellBasedPdeHandler<3>(cell_population.get()));
        pde_handler->AddPdeAndBc(&pde_and_bc);
        ChastePoint<3> lower(0.0, 0.0, 0.0);
        ChastePoint<3> upper(10.0, 50.0, 50.0);
        ChasteCuboid<3> cuboid(lower, upper);
        pde_handler->UseCoarsePdeMesh(2.0, cuboid, true);
        pde_handler->SetImposeBcsOnCoarseBoundary(true);

        cell_population->AddPdeHandler(pde_handler);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(8, 8);

        while (!(SimulationTime::Instance()->IsFinished()))
        {
            pde_handler->OpenResultsFiles(output_directory);
            pde_handler->SolvePdeAndWriteResultsToFile(1);
            pde_handler->CloseResultsFiles();

            cell_population->UpdateVascularCellPopulation();

            std::string output_filename_vessels = output_file_handler.GetOutputDirectoryFullPath().append("VesselNetwork_");
            p_network->Write(output_filename_vessels + boost::lexical_cast<std::string>(SimulationTime::Instance()->GetTimeStepsElapsed()) + ".vtp");

            cell_population->OpenWritersFiles(output_file_handler);
            cell_population->WriteResultsToFiles(output_directory);
            cell_population->CloseWritersFiles();

            SimulationTime::Instance()->IncrementTimeOneStep();
        }
    }
};

#endif /*TESTCABASEDCELLPOPULATIONWITHVESSELS_HPP*/
