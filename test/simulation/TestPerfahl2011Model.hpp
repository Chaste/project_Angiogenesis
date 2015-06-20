//
//  TestPerfahl2011Model.hpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 03/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#ifndef TestPerfahl2011Model_hpp
#define TestPerfahl2011Model_hpp

#include <cxxtest/TestSuite.h>
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VasculatureGenerator.hpp"
#include "VascularNode.hpp"
#include "CaVesselSegment.hpp"
#include "CaVessel.hpp"
#include "CaVascularNetwork.hpp"
#include "AngiogenesisSolver.hpp"
#include "OnLatticeVascularTumourCellPopulationGenerator.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellLabelWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "NodeLocationWriter.hpp"

#include "SimulationTime.hpp"
#include "PottsMesh.hpp"
#include "PottsMeshGenerator.hpp"
#include "FakePetscSetup.hpp"

#include "Debug.hpp"

class TestPerfahl2011Model : public AbstractCellBasedWithTimingsTestSuite
{

public:


    void TestSetUpSingleVessel() throw(Exception)
    {

        // initialise domain
        // _________________

        unsigned numNodesAcross = 20;
        unsigned numElementsAcross = 0;
        unsigned elementWidth = 0;
        unsigned numNodesUp = 20;
        unsigned numElementsUp=0;
        unsigned elementHeight=0;
        unsigned numNodesDeep = 20;
        unsigned numElementsDeep=0;
        unsigned elementDepth=0;
        bool startAtBottomLeft = false;
        bool isPeriodicInX = false;
        bool isPeriodicInY = false;
        bool isPeriodicInZ = false;

        const unsigned dimensionality = 3;
        //        double dx = 2e-5; // metres

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<dimensionality> generator( numNodesAcross, numElementsAcross, elementWidth,
                numNodesUp, numElementsUp, elementHeight,
                numNodesDeep, numElementsDeep, elementDepth,
                startAtBottomLeft, isPeriodicInX, isPeriodicInY , isPeriodicInZ);
        PottsMesh<dimensionality>* p_mesh = generator.GetMesh();

        // Create vascular network
        // _______________________

        // initially we will just consider a single vessel through the middle of the domain, extending in the z-direction
        c_vector<double, dimensionality> startPosition;
        startPosition[0] = 10;
        startPosition[1] = 10;
        startPosition[2] = 0;
        VasculatureGenerator<dimensionality> vascular_network_generator;
        boost::shared_ptr<CaVascularNetwork<dimensionality> > vascular_network =
                vascular_network_generator.GenerateSingleVessel(numNodesUp-1,startPosition);

        // Write the network to file
        std::string output_directory = "TestPerfahl2011Model";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("InitialVesselNetwork.vtp");
        vascular_network->Write(output_filename);

        // Create cell population
        // ______________________

        // use OnLatticeVascularTumourCellPopulationGenerator to generate cells and associate cells
        // with vessels
        OnLatticeVascularTumourCellPopulationGenerator<dimensionality> cellPopulationGenerator;
        boost::shared_ptr<CaBasedCellPopulation<dimensionality> > cell_population =
                cellPopulationGenerator.CreateCellPopulation(*p_mesh, vascular_network);

        // Test that there are 20 nodes in the network, 2 vessel nodes, and 1 vessel
        TS_ASSERT_EQUALS(vascular_network->GetNumberOfNodes(),20u);
        TS_ASSERT_EQUALS(vascular_network->GetNumberOfVessels(),1u);
        TS_ASSERT_EQUALS(vascular_network->GetNumberOfVesselNodes(),2u);
        // Test that each node has a cell associated with it
        for (unsigned i = 0; i < vascular_network->GetNumberOfNodes(); i++)
        {
            TS_ASSERT(vascular_network->GetNode(i)->HasCell());
        }

        // VTK writing needs a simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        cell_population->AddCellWriter<CellLabelWriter>();
        cell_population->AddCellWriter<CellMutationStatesWriter>();
        cell_population->AddPopulationWriter<NodeLocationWriter>();
        cell_population->OpenWritersFiles(output_file_handler);
        cell_population->WriteResultsToFiles(output_directory);
        cell_population->CloseWritersFiles();

    }

    // to be completed
    void TestSetUpAndSolveOxygenPDE() throw(Exception)
        {

            // initialise domain
            // _________________

            unsigned numNodesAcross = 20;
            unsigned numElementsAcross = 0;
            unsigned elementWidth = 0;
            unsigned numNodesUp = 20;
            unsigned numElementsUp=0;
            unsigned elementHeight=0;
            unsigned numNodesDeep = 20;
            unsigned numElementsDeep=0;
            unsigned elementDepth=0;
            bool startAtBottomLeft = false;
            bool isPeriodicInX = false;
            bool isPeriodicInY = false;
            bool isPeriodicInZ = false;

            const unsigned dimensionality = 3;
            //        double dx = 2e-5; // metres

            // Create a simple 2D PottsMesh
            PottsMeshGenerator<dimensionality> generator( numNodesAcross, numElementsAcross, elementWidth,
                    numNodesUp, numElementsUp, elementHeight,
                    numNodesDeep, numElementsDeep, elementDepth,
                    startAtBottomLeft, isPeriodicInX, isPeriodicInY , isPeriodicInZ);
            PottsMesh<dimensionality>* p_mesh = generator.GetMesh();

            // Create vascular network
            // _______________________

            // initially we will just consider a single vessel through the middle of the domain, extending in the z-direction
            c_vector<double, dimensionality> startPosition;
            startPosition[0] = 10;
            startPosition[1] = 10;
            startPosition[2] = 0;
            VasculatureGenerator<dimensionality> vascular_network_generator;
            boost::shared_ptr<CaVascularNetwork<dimensionality> > vascular_network =
                    vascular_network_generator.GenerateSingleVessel(numNodesUp-1,startPosition);

            // Write the network to file
            std::string output_directory = "TestPerfahl2011Model";
            OutputFileHandler output_file_handler(output_directory, false);
            std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("InitialVesselNetwork.vtp");
            vascular_network->Write(output_filename);

            // Create cell population
            // ______________________

            // use OnLatticeVascularTumourCellPopulationGenerator to generate cells and associate cells
            // with vessels
            OnLatticeVascularTumourCellPopulationGenerator<dimensionality> cellPopulationGenerator;
            boost::shared_ptr<CaBasedCellPopulation<dimensionality> > cell_population =
                    cellPopulationGenerator.CreateCellPopulation(*p_mesh, vascular_network);

            // Test that there are 20 nodes in the network, 2 vessel nodes, and 1 vessel
            TS_ASSERT_EQUALS(vascular_network->GetNumberOfNodes(),20u);
            TS_ASSERT_EQUALS(vascular_network->GetNumberOfVessels(),1u);
            TS_ASSERT_EQUALS(vascular_network->GetNumberOfVesselNodes(),2u);
            // Test that each node has a cell associated with it
            for (unsigned i = 0; i < vascular_network->GetNumberOfNodes(); i++)
            {
                TS_ASSERT(vascular_network->GetNode(i)->HasCell());
            }

            // VTK writing needs a simulation time
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            cell_population->AddCellWriter<CellLabelWriter>();
            cell_population->AddCellWriter<CellMutationStatesWriter>();
            cell_population->AddPopulationWriter<NodeLocationWriter>();
            cell_population->OpenWritersFiles(output_file_handler);
            cell_population->WriteResultsToFiles(output_directory);
            cell_population->CloseWritersFiles();


        }

//    void TestSetUpSingleVessel() throw(Exception)
//        {
//
//            // initialise domain
//            // _________________
//
//            unsigned numNodesAcross = 20;
//            unsigned numElementsAcross = 0;
//            unsigned elementWidth = 0;
//            unsigned numNodesUp = 20;
//            unsigned numElementsUp=0;
//            unsigned elementHeight=0;
//            unsigned numNodesDeep = 20;
//            unsigned numElementsDeep=0;
//            unsigned elementDepth=0;
//            bool startAtBottomLeft = false;
//            bool isPeriodicInX = false;
//            bool isPeriodicInY = false;
//            bool isPeriodicInZ = false;
//
//            const unsigned dimensionality = 3;
//            //        double dx = 2e-5; // metres
//
//            // Create a simple 2D PottsMesh
//            PottsMeshGenerator<dimensionality> generator( numNodesAcross, numElementsAcross, elementWidth,
//                    numNodesUp, numElementsUp, elementHeight,
//                    numNodesDeep, numElementsDeep, elementDepth,
//                    startAtBottomLeft, isPeriodicInX, isPeriodicInY , isPeriodicInZ);
//            PottsMesh<dimensionality>* p_mesh = generator.GetMesh();
//
//            // Create vascular network
//            // _______________________
//
//            // initially we will just consider a single vessel through the middle of the domain, extending in the z-direction
//            c_vector<double, dimensionality> startPosition;
//            startPosition[0] = 10;
//            startPosition[1] = 10;
//            startPosition[2] = 0;
//            VasculatureGenerator<dimensionality> vascular_network_generator;
//            boost::shared_ptr<CaVascularNetwork<dimensionality> > vascular_network =
//                    vascular_network_generator.GenerateSingleVessel(numNodesUp-1,startPosition);
//
//            // Write the network to file
//            std::string output_directory = "TestPerfahl2011Model";
//            OutputFileHandler output_file_handler(output_directory, false);
//            std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("InitialVesselNetwork.vtp");
//            vascular_network->Write(output_filename);
//
//            // Create cell population
//            // ______________________
//
//            // use OnLatticeVascularTumourCellPopulationGenerator to generate cells and associate cells
//            // with vessels
//            OnLatticeVascularTumourCellPopulationGenerator<dimensionality> cellPopulationGenerator;
//            boost::shared_ptr<CaBasedCellPopulation<dimensionality> > cell_population =
//                    cellPopulationGenerator.CreateCellPopulation(*p_mesh, vascular_network);
//
//            // Test that there are 20 nodes in the network, 2 vessel nodes, and 1 vessel
//            TS_ASSERT_EQUALS(vascular_network->GetNumberOfNodes(),20u);
//            TS_ASSERT_EQUALS(vascular_network->GetNumberOfVessels(),1u);
//            TS_ASSERT_EQUALS(vascular_network->GetNumberOfVesselNodes(),2u);
//            // Test that each node has a cell associated with it
//            for (unsigned i = 0; i < vascular_network->GetNumberOfNodes(); i++)
//            {
//                TS_ASSERT(vascular_network->GetNode(i)->HasCell());
//            }
//
//            // VTK writing needs a simulation time
//            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
//
//            cell_population->AddCellWriter<CellLabelWriter>();
//            cell_population->AddCellWriter<CellMutationStatesWriter>();
//            cell_population->AddPopulationWriter<NodeLocationWriter>();
//            cell_population->OpenWritersFiles(output_file_handler);
//            cell_population->WriteResultsToFiles(output_directory);
//            cell_population->CloseWritersFiles();
//
//
//            //  // initialize structural adaptation algorithm
//            //        // __________________________________________
//            //
//            //        boost::shared_ptr<SimpleStructuralAdaptationAlgorithm<dimensionality,int> > StructuralAdaptationAlgorithm(new SimpleStructuralAdaptationAlgorithm<dimensionality,int>());
//            //        // To prevent the vessel radii blowing up or shrinking too far we set maximum and minimum vessel radii here
//            //        // Although this was not done originally in the 2003 paper, it was done in later papers (e.g. Owen et al. (2009), Angiogenesis and
//            //        // vascular remodelling in normal and cancerous tissues
//            //        boost::shared_ptr<Alarcon03RadiusCalculation <dimensionality,int> > radiusCalculation(new Alarcon03RadiusCalculation <dimensionality,int>());
//            //        radiusCalculation->SetMaxRadius(5e-5);
//            //        radiusCalculation->SetMinRadius(1e-6);
//            //        StructuralAdaptationAlgorithm->SetRadiusCalculation(radiusCalculation);
//            //        boost::shared_ptr<ConstantHaematocritCalculation<dimensionality,int> > haematocritCalculator(new ConstantHaematocritCalculation<dimensionality,int>());
//            //        StructuralAdaptationAlgorithm->SetHaematocritCalculation(haematocritCalculator);
//            //
//            //        StructuralAdaptationAlgorithm->SetMaximumNumberOfIterations(1000);
//            //
//            //
//            //        cout << "Structural adaptation algorithm instantiated and customised\n";
//            //
//            //
//            //
//            //        // initialise cellpopulation
//            //        // _________________________
//            //
//            //        boost::shared_ptr<CellPopulation<dimensionality,int> > cellpopulation(new CellPopulation<dimensionality,int>(spatialmesh));
//            //
//            //        for (int i = 0; i < DomainSize_X; i++)
//            //        {
//            //            for (int j = 0; j < DomainSize_Y; j++)
//            //            {
//            //                for (int k = 0; k < DomainSize_Z; k++)
//            //                {
//            //
//            //                    boost::shared_ptr<NormalCell<dimensionality,int> > newcell(new NormalCell<dimensionality,int>(SpatialCoordinate<dimensionality,int>(i,j,k)));
//            //
//            //                    newcell->AddCellularChemical("p53",Concentration(0.0,"unitless"),0.0,0.0);
//            //                    newcell->AddCellularChemical("VEGF",Concentration(0.0,"unitless"),0.0,0.0);
//            //                    newcell->AddCellularChemical("Oxygen",Concentration(0.0,"mmHg"),0.13,0.0);
//            //
//            //                    if (cellpopulation->CellMayMoveInToLocation(SpatialCoordinate<dimensionality,int>(i,j,k),newcell))
//            //                    {
//            //                        cellpopulation->AddCell(newcell);
//            //                    }
//            //
//            //                }
//            //
//            //            }
//            //        }
//            //
//            //        // input small square tumour
//            //
//            //        for (int i = 10; i < 15; i++)
//            //        {
//            //            for (int j = 20; j < 25; j++)
//            //            {
//            //                for (int k = 0; k < 4; k++)
//            //                {
//            //
//            //                    for (int l = 0; l < cellpopulation->GetNumberOfCellsAtLocation(SpatialCoordinate<dimensionality,int>(i,j,k)); l++)
//            //                    {
//            //                        cellpopulation->CellDeath(cellpopulation->GetCell(SpatialCoordinate<dimensionality,int>(i,j,k),l));
//            //                    }
//            //
//            //                    boost::shared_ptr<CancerousCell<dimensionality,int> > newcell(new CancerousCell<dimensionality,int>(SpatialCoordinate<dimensionality,int>(i,j,k)));
//            //
//            //                    newcell->AddCellularChemical("p53",Concentration(0.0,"unitless"),0.0,0.0);
//            //                    newcell->AddCellularChemical("VEGF",Concentration(0.0,"unitless"),0.0,0.0);
//            //                    newcell->AddCellularChemical("Oxygen",Concentration(0.0,"mmHg"),1.3,0.0);
//            //
//            //                    newcell->SetMaximumTimeInQuiescentState(6000);
//            //
//            //                    if (cellpopulation->CellMayMoveInToLocation(SpatialCoordinate<dimensionality,int>(i,j,k),newcell))
//            //                    {
//            //                        cellpopulation->AddCell(newcell);
//            //                    }
//            //
//            //                }
//            //            }
//            //        }
//            //
//            //        // instantiate oxygen calculator
//            //        // _____________________________
//            //
//            //
//            //        boost::shared_ptr<Owen2011OxygenCalculator<dimensionality,int> > oxygenCalculator(new Owen2011OxygenCalculator<dimensionality,int>(cellpopulation,vascularnetwork,"mmHg"));
//            //
//            //
//            //
//            //
//            //        // instantiate VEGF calculator
//            //        // _____________________________
//            //
//            //        boost::shared_ptr<Owen2011VEGFCalculator<dimensionality,int> > vegfCalculator(new Owen2011VEGFCalculator<dimensionality,int>(cellpopulation,vascularnetwork,"unitless"));
//            //
//            //
//            //
//            //
//            //
//            //        // initialize vessel regression model            oxygenCalculator->GetPde()->SetHaematocritOxygenConcentrationConversionFactor(20.0/vascularnetwork->GetArterialHaematocritLevel());
//            //        // \todo need to decide on consistent locations for defining uptake rates for cells and permeabilities for vessels
//            //        oxygenCalculator->GetPde()->SetConstantOxygenUptakeRateForNormalCell(20);
//            //        oxygenCalculator->GetPde()->SetConstantOxygenUptakeRateForCancerousCell(25);
//            //        // __________________________________
//            //
//            //        boost::shared_ptr<Owen09WallShearStressDependentVesselRegressionModel<dimensionality,int> > vesselRegressionModel(new Owen09WallShearStressDependentVesselRegressionModel<dimensionality,int>());
//            //
//            //
//            //        // initialize angiogenesis model
//            //        // _____________________________
//            //
//            //        boost::shared_ptr<Owen09AngiogenesisModel<dimensionality,int> > angiogenesisModel(new Owen09AngiogenesisModel<dimensionality,int>(diffusibles));
//            //        angiogenesisModel->SetRadiusOfExclusion(20e-6);
//            //
//            //
//            //        cout << "Angiogenesis model instantiated\n";
//            //
//            //
//            //
//            //
//            //
//            //        // instantiate fixed duration sub cellular model
//            //        // _____________________________________________
//            //
//            //        boost::shared_ptr<Owen2011SubCellularModel<dimensionality,int> > subCellularModel(new Owen2011SubCellularModel<dimensionality,int>(cellpopulation));
//            ////            boost::shared_ptr<Owen2011QuiescentStateEntryStrategy<dimensionality,int> > quiescentStateEntryStrategy(new Owen2011QuiescentStateEntryStrategy<dimensionality,int>());
//            ////            quiescentStateEntryStrategy->SetEnterQuiescenceO2Concentration(1.8);
//            ////            boost::shared_ptr<Owen2011ProliferativeStateEntryStrategy<dimensionality,int> > proliferativeStateEntryStrategy(new Owen2011ProliferativeStateEntryStrategy<dimensionality,int>());
//            ////            proliferativeStateEntryStrategy->SetEnterProliferativeStateO2Concentration(2);
//            ////
//            ////            subCellularModel->SetQuiescentStateEntryStrategy(quiescentStateEntryStrategy);
//            ////            subCellularModel->SetProliferativeStateEntryStrategy(proliferativeStateEntryStrategy);
//            //
//            //        // instantiate oxygen based cell proliferator
//            //        // __________________________________________
//            //
//            //        boost::shared_ptr<OxygenBasedCellProliferator<dimensionality,int> > proliferator(new OxygenBasedCellProliferator<dimensionality,int>(cellpopulation,diffusibles));
//            //
//            //
//            //        // set up and run simulation
//            //        // _________________________
//            //
//            //
//            //        std::string resultsDirectoryName = "CancerInformaticsPaper/MockSimulationResults";
//            //
//            //        bool overwritePreviousSimulation = false;
//            //        VascularTumourGrowthSimulation_AlarconFamily<dimensionality,int> simulation(cellpopulation,vascularnetwork,diffusibles,resultsDirectoryName,overwritePreviousSimulation);
//            //
//            //        simulation.AddOxygenCalculator(oxygenCalculator);
//            //        simulation.AddVEGFCalculator(vegfCalculator);
//            //        simulation.AddStructuralAdaptationAlgorithm(StructuralAdaptationAlgorithm);
//            //        simulation.AddSubCellularModel(subCellularModel);
//            //        simulation.AddCellProliferator(proliferator);
//            //        simulation.AddAngiogenesisModel(angiogenesisModel);
//            ////            simulation.AddVesselRegressionModel(vesselRegressionModel);
//            //        simulation.SetRandomlyInitialiseCellCycleTimes(true);
//            //        simulation.SetMaximumInitialCellCycleTime(3000);
//            //        simulation.SetDuration(40000); // 7 days ~ 10000 minutes
//            //        simulation.SetApproximateTimeStep(30);
//            //        simulation.SetVerbose(true);
//            //
//            //        simulation.Run();
//
//        }

};

#endif
