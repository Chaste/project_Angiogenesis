//
//  TestSimpleStructuralAdaptationSolver.hpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 03/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#ifndef TestAbstractAngiogenesisSolver_hpp
#define TestAbstractAngiogenesisSolver_hpp

#include <cxxtest/TestSuite.h>
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VasculatureGenerator.hpp"
#include "VascularNode.hpp"
#include "CaVesselSegment.hpp"
#include "CaVessel.hpp"
#include "CaVascularNetwork.hpp"
#include "Part.hpp"
#include "FiniteDifferenceSolver.hpp"
#include "AbstractAngiogenesisSolver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CaVesselSegment.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractGrowthDirectionModifier.hpp"

class TestAbstractAngiogenesisSolver : public AbstractCellBasedTestSuite
{

public:

    void TestSingleVesselGrowth() throw(Exception)
    {
        // Make a network
        boost::shared_ptr<VascularNode<3> > p_node1 = VascularNode<3>::Create(0.0, 0.0, 0.0);
        boost::shared_ptr<VascularNode<3> > p_node2 = VascularNode<3>::Create(100.0, 0.0, 0.0);
        p_node2->SetIsMigrating(true);
        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(CaVesselSegment<3>::Create(p_node1, p_node2));

        boost::shared_ptr<CaVascularNetwork<3> > p_network = boost::shared_ptr<CaVascularNetwork<3> >(new CaVascularNetwork<3>());
        p_network->AddVessel(p_vessel1);
        // Set the one of the nodes to migrate
        p_network->GetVessel(0)->GetEndNode()->SetIsMigrating(true);

        OutputFileHandler output_file_handler("TestAbstractAngiogenesisSolver/SingleVesselGrowth/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        boost::shared_ptr<AbstractGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<AbstractGrowthDirectionModifier<3> >(new AbstractGrowthDirectionModifier<3>());

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AbstractAngiogenesisSolver<3> angiogenesis_solver(p_network);
        angiogenesis_solver.SetOutputDirectory(output_directory);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.Run();

        // Test that the vessel has grown
        TS_ASSERT(p_network->GetNumberOfNodes() == 12u);
        TS_ASSERT_DELTA(p_network->GetVessel(0)->GetEndNode()->GetLocationVector()[0], 200.0, 1.e-6);
    }

    void TestSingleVesselGrowthFromStart() throw(Exception)
    {
        // Make a network
        boost::shared_ptr<VascularNode<3> > p_node1 = VascularNode<3>::Create(0.0, 0.0, 0.0);
        boost::shared_ptr<VascularNode<3> > p_node2 = VascularNode<3>::Create(100.0, 0.0, 0.0);
        p_node1->SetIsMigrating(true);
        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(CaVesselSegment<3>::Create(p_node1, p_node2));

        boost::shared_ptr<CaVascularNetwork<3> > p_network = boost::shared_ptr<CaVascularNetwork<3> >(new CaVascularNetwork<3>());
        p_network->AddVessel(p_vessel1);
        // Set the one of the nodes to migrate
        p_network->GetVessel(0)->GetStartNode()->SetIsMigrating(true);

        OutputFileHandler output_file_handler("TestAbstractAngiogenesisSolver/SingleVesselGrowthFromStart/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        boost::shared_ptr<AbstractGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<AbstractGrowthDirectionModifier<3> >(new AbstractGrowthDirectionModifier<3>());

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AbstractAngiogenesisSolver<3> angiogenesis_solver(p_network);
        angiogenesis_solver.SetOutputDirectory(output_directory);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.Run();

        // Test that the vessel has grown
        TS_ASSERT(p_network->GetNumberOfNodes() == 12u);
        TS_ASSERT_DELTA(p_network->GetVessel(0)->GetStartNode()->GetLocationVector()[0], -100.0, 1.e-6);
    }

    void TestTipTipAnastamosisEvent() throw(Exception)
    {
        // Make a network: two vessels on a collision course
        boost::shared_ptr<VascularNode<3> > p_node1 = VascularNode<3>::Create(0.0, 0.0, 0.0);
        boost::shared_ptr<VascularNode<3> > p_node2 = VascularNode<3>::Create(100.0, 0.0, 0.0);
        p_node2->SetIsMigrating(true);
        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(CaVesselSegment<3>::Create(p_node1, p_node2));

        boost::shared_ptr<VascularNode<3> > p_node3 = VascularNode<3>::Create(200.0, 0.0, 0.0);
        boost::shared_ptr<VascularNode<3> > p_node4 = VascularNode<3>::Create(300.0, 0.0, 0.0);
        p_node3->SetIsMigrating(true);
        boost::shared_ptr<CaVessel<3> > p_vessel2 = CaVessel<3>::Create(CaVesselSegment<3>::Create(p_node3, p_node4));

        boost::shared_ptr<CaVascularNetwork<3> > p_network = boost::shared_ptr<CaVascularNetwork<3> >(new CaVascularNetwork<3>());
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);

        boost::shared_ptr<AbstractGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<AbstractGrowthDirectionModifier<3> >(new AbstractGrowthDirectionModifier<3>());

        OutputFileHandler output_file_handler("TestAbstractAngiogenesisSolver/TipTipAnastamosisEvent/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AbstractAngiogenesisSolver<3> angiogenesis_solver(p_network);
        angiogenesis_solver.SetOutputDirectory(output_directory);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.Run();

        // Test that the vessel has grown
        TS_ASSERT_DELTA(p_network->GetVessel(0)->GetEndNode()->GetLocationVector()[0], 150.0, 1.e-6);
        TS_ASSERT_DELTA(p_network->GetVessel(1)->GetStartNode()->GetLocationVector()[0], 150.0, 1.e-6);
    }

    void TestTipSproutAnastamosisEvent() throw(Exception)
    {
        // Make a network: two vessels on a collision course
        boost::shared_ptr<VascularNode<3> > p_node1 = VascularNode<3>::Create(0.0, 0.0, 0.0);
        boost::shared_ptr<VascularNode<3> > p_node2 = VascularNode<3>::Create(100.0, 0.0, 0.0);
        p_node2->SetIsMigrating(true);
        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(CaVesselSegment<3>::Create(p_node1, p_node2));

        boost::shared_ptr<VascularNode<3> > p_node3 = VascularNode<3>::Create(150.0, 50.0, 0.0);
        boost::shared_ptr<VascularNode<3> > p_node4 = VascularNode<3>::Create(150.0, -50.0, 0.0);
        boost::shared_ptr<CaVessel<3> > p_vessel2 = CaVessel<3>::Create(CaVesselSegment<3>::Create(p_node3, p_node4));

        boost::shared_ptr<CaVascularNetwork<3> > p_network = boost::shared_ptr<CaVascularNetwork<3> >(new CaVascularNetwork<3>());
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);

        boost::shared_ptr<AbstractGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<AbstractGrowthDirectionModifier<3> >(new AbstractGrowthDirectionModifier<3>());

        OutputFileHandler output_file_handler("TestAbstractAngiogenesisSolver/TipSproutAnastamosisEvent/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AbstractAngiogenesisSolver<3> angiogenesis_solver(p_network);
        angiogenesis_solver.SetOutputDirectory(output_directory);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.Run();

        // Test that the vessel has grown
        TS_ASSERT_DELTA(p_network->GetVessel(0)->GetEndNode()->GetLocationVector()[0], 150.0, 1.e-6);
        TS_ASSERT_DELTA(p_network->GetVessel(1)->GetStartNode()->GetLocationVector()[0], 150.0, 1.e-6);
    }

    void TestMultiVessel() throw(Exception)
    {
        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        for(unsigned idx=0; idx<10; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 0.0, 0.0));
        }

        std::vector<boost::shared_ptr<VascularNode<3> > > top_nodes;
        for(unsigned idx=0; idx<10; idx++)
        {
            top_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 100.0, 0.0));
        }

        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVessel<3> > p_vessel2 = CaVessel<3>::Create(top_nodes);

        boost::shared_ptr<AbstractGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<AbstractGrowthDirectionModifier<3> >(new AbstractGrowthDirectionModifier<3>());

        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);
        p_network->FormSprout(ChastePoint<3>(20.0, 0.0, 0.0), ChastePoint<3>(20.0, 10.0, 0.0));
        p_network->FormSprout(ChastePoint<3>(70.0, 100.0, 0.0), ChastePoint<3>(70.0, 90.0, 0.0));

        OutputFileHandler output_file_handler("TestAbstractAngiogenesisSolver/Multisegment/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AbstractAngiogenesisSolver<3> angiogenesis_solver(p_network);
        angiogenesis_solver.SetOutputDirectory(output_directory);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.Run();
    }

    void TestMultiSprout() throw(Exception)
    {
        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        for(unsigned idx=0; idx<10; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 0.0, 0.0));
        }

        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);

        boost::shared_ptr<AbstractGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<AbstractGrowthDirectionModifier<3> >(new AbstractGrowthDirectionModifier<3>());

        std::vector<boost::shared_ptr<VascularNode<3> > > top_nodes;
        for(unsigned idx=1; idx<9; idx++)
        {
            p_network->FormSprout(ChastePoint<3>(double(idx)*10, 0.0, 0.0), ChastePoint<3>(double(idx)*10, 10.0, 0.0));
        }

        OutputFileHandler output_file_handler("TestAbstractAngiogenesisSolver/MultiSprout/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AbstractAngiogenesisSolver<3> angiogenesis_solver(p_network);
        angiogenesis_solver.SetOutputDirectory(output_directory);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.Run();
    }

    void TestMultiSproutWithPde() throw(Exception)
    {
        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        for(unsigned idx=0; idx<11; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 10.0, 0.0));
        }

        std::vector<boost::shared_ptr<VascularNode<3> > > top_nodes;
        for(unsigned idx=0; idx<11; idx++)
        {
            top_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 90.0, 0.0));
        }

        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVessel<3> > p_vessel2 = CaVessel<3>::Create(top_nodes);
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);

        for(unsigned idx=1; idx<10; idx++)
        {
            p_network->FormSprout(ChastePoint<3>(double(idx)*10, 10.0, 0.0), ChastePoint<3>(double(idx)*10, 20.0, 0.0));
        }

        boost::shared_ptr<AbstractGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<AbstractGrowthDirectionModifier<3> >(new AbstractGrowthDirectionModifier<3>());

        // Set up the PDE domain
        boost::shared_ptr<Part<3> > p_domain = Part<3>::Create();
        p_domain->AddCuboid(100, 100, 10);

        // Choose the PDE
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_pde = HybridLinearEllipticPde<3>::Create();
        p_pde->SetDiffusionConstant(0.0033);
        p_pde->SetLinearInUTerm(-2.e-7);

        // Set up and run the solver
        boost::shared_ptr<FiniteDifferenceSolver<3> > p_solver = FiniteDifferenceSolver<3>::Create();
        p_solver->SetExtents(p_domain, 10.0);
        p_solver->SetPde(p_pde);

        OutputFileHandler output_file_handler("TestAbstractAngiogenesisSolver/MultiSproutPde/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AbstractAngiogenesisSolver<3> angiogenesis_solver(p_network);
        angiogenesis_solver.SetOutputDirectory(output_directory);
        angiogenesis_solver.AddPdeSolver(p_solver);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.Run();
    }

    void TestMultiSproutWithFlow() throw(Exception)
    {
        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        for(unsigned idx=0; idx<11; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 10.0, 0.0));
        }
        bottom_nodes[0]->GetFlowProperties()->SetIsInputNode(true);
        bottom_nodes[0]->GetFlowProperties()->SetPressure(3000);
        bottom_nodes[10]->GetFlowProperties()->SetIsOutputNode(true);
        bottom_nodes[10]->GetFlowProperties()->SetPressure(1000);

        std::vector<boost::shared_ptr<VascularNode<3> > > top_nodes;
        for(unsigned idx=0; idx<11; idx++)
        {
            top_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 90.0, 0.0));
        }
        top_nodes[10]->GetFlowProperties()->SetIsInputNode(true);
        top_nodes[10]->GetFlowProperties()->SetPressure(3000);
        top_nodes[0]->GetFlowProperties()->SetIsOutputNode(true);
        top_nodes[0]->GetFlowProperties()->SetPressure(1000);

        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVessel<3> > p_vessel2 = CaVessel<3>::Create(top_nodes);
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);
        p_network->SetSegmentRadii(10.0);

        for(unsigned idx=1; idx<10; idx++)
        {
            p_network->FormSprout(ChastePoint<3>(double(idx)*10, 10.0, 0.0), ChastePoint<3>(double(idx)*10, 20.0, 0.0));
        }

        boost::shared_ptr<AbstractGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<AbstractGrowthDirectionModifier<3> >(new AbstractGrowthDirectionModifier<3>());

        p_network->UpdateSegments();
        std::vector<boost::shared_ptr<CaVesselSegment<3> > > segments = p_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(1.e-3);
        }

        OutputFileHandler output_file_handler("TestAbstractAngiogenesisSolver/MultiSproutFlow/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AbstractAngiogenesisSolver<3> angiogenesis_solver(p_network);
        angiogenesis_solver.SetOutputDirectory(output_directory);
        angiogenesis_solver.SetSolveFlow();
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.Run();
    }

    void TestSproutingWithFlow() throw(Exception)
    {
        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        for(unsigned idx=0; idx<9; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 10.0, 0.0));
        }
        bottom_nodes[0]->GetFlowProperties()->SetIsInputNode(true);
        bottom_nodes[0]->GetFlowProperties()->SetPressure(3000);
        bottom_nodes[8]->GetFlowProperties()->SetIsOutputNode(true);
        bottom_nodes[8]->GetFlowProperties()->SetPressure(1000);

        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);
        p_network->SetSegmentRadii(10.0);

        for(unsigned idx=1; idx<6; idx+=2)
        {
            p_network->FormSprout(ChastePoint<3>(double(idx)*10, 10.0, 0.0), ChastePoint<3>(double(idx)*10, 20.0, 0.0));
        }

        p_network->UpdateSegments();
        std::vector<boost::shared_ptr<CaVesselSegment<3> > > segments = p_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(1.e-3);
        }

        boost::shared_ptr<AbstractGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<AbstractGrowthDirectionModifier<3> >(new AbstractGrowthDirectionModifier<3>());

        OutputFileHandler output_file_handler("TestAbstractAngiogenesisSolver/SproutingFlow/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AbstractAngiogenesisSolver<3> angiogenesis_solver(p_network);
        angiogenesis_solver.SetOutputDirectory(output_directory);
        angiogenesis_solver.SetSolveFlow();
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.SetSproutingProbability(0.5);
        angiogenesis_solver.Run();
    }
};

#endif
