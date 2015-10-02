//
//  TestSimpleStructuralAdaptationSolver.hpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 03/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#ifndef TestOffLatticePrwGrowthDirectionModifier_hpp
#define TestOffLatticePrwGrowthDirectionModifier_hpp

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
#include "SimpleFlowSolver.hpp"
#include "OffLatticePrwGrowthDirectionModifier.hpp"

class TestOffLatticePrwGrowthDirectionModifier : public AbstractCellBasedTestSuite
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

        boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> >(new OffLatticePrwGrowthDirectionModifier<3>());

        OutputFileHandler output_file_handler("TestOffLatticePrwGrowthDirectionModifier/SingleVesselGrowth/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AbstractAngiogenesisSolver<3> angiogenesis_solver(p_network);
        angiogenesis_solver.SetOutputDirectory(output_directory);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.Run();
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

        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);
        p_network->FormSprout(ChastePoint<3>(20.0, 0.0, 0.0), ChastePoint<3>(20.0, 10.0, 0.0));
        p_network->FormSprout(ChastePoint<3>(70.0, 100.0, 0.0), ChastePoint<3>(70.0, 90.0, 0.0));

        boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> >(new OffLatticePrwGrowthDirectionModifier<3>());

        OutputFileHandler output_file_handler("TestOffLatticePrwGrowthDirectionModifier/Multisegment/", false);
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

        std::vector<boost::shared_ptr<VascularNode<3> > > top_nodes;
        for(unsigned idx=1; idx<9; idx++)
        {
            p_network->FormSprout(ChastePoint<3>(double(idx)*10, 0.0, 0.0), ChastePoint<3>(double(idx)*10, 10.0, 0.0));
        }

        boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> >(new OffLatticePrwGrowthDirectionModifier<3>());

        OutputFileHandler output_file_handler("TestOffLatticePrwGrowthDirectionModifier/MultiSprout/", false);
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

        boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> >(new OffLatticePrwGrowthDirectionModifier<3>());

        OutputFileHandler output_file_handler("TestOffLatticePrwGrowthDirectionModifier/MultiSproutPde/", false);
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

        p_network->UpdateSegments();
        std::vector<boost::shared_ptr<CaVesselSegment<3> > > segments = p_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(1.e-3);
        }

        boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> >(new OffLatticePrwGrowthDirectionModifier<3>());

        OutputFileHandler output_file_handler("TestOffLatticePrwGrowthDirectionModifier/MultiSproutFlow/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AbstractAngiogenesisSolver<3> angiogenesis_solver(p_network);
        angiogenesis_solver.SetOutputDirectory(output_directory);

        boost::shared_ptr<SimpleFlowSolver<3> > p_simple_flow_solver =
                boost::shared_ptr<SimpleFlowSolver<3> >(new SimpleFlowSolver<3>());

        angiogenesis_solver.SetFlowSolver(p_simple_flow_solver);
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

        boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> > p_grow_direction_modifier =
                boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> >(new OffLatticePrwGrowthDirectionModifier<3>());

        OutputFileHandler output_file_handler("TestOffLatticePrwGrowthDirectionModifier/SproutingFlow/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10, 10);
        AbstractAngiogenesisSolver<3> angiogenesis_solver(p_network);
        angiogenesis_solver.SetOutputDirectory(output_directory);
        boost::shared_ptr<SimpleFlowSolver<3> > p_simple_flow_solver =SimpleFlowSolver<3>::Create();

        angiogenesis_solver.SetFlowSolver(p_simple_flow_solver);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.Run();
    }
};

#endif
