//
//  TestSimpleStructuralAdaptationSolver.hpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 03/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#ifndef TestAngiogenesisSolver_hpp
#define TestAngiogenesisSolver_hpp

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

class TestAngiogenesisSolver : public CxxTest::TestSuite
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

        OutputFileHandler output_file_handler("TestAngiogenesisSolver/SingleVesselGrowth/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        AngiogenesisSolver<3> angiogenesis_solver(p_network, output_directory);
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

        OutputFileHandler output_file_handler("TestAngiogenesisSolver/SingleVesselGrowthFromStart/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        AngiogenesisSolver<3> angiogenesis_solver(p_network, output_directory);
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

        OutputFileHandler output_file_handler("TestAngiogenesisSolver/TipTipAnastamosisEvent/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        AngiogenesisSolver<3> angiogenesis_solver(p_network, output_directory);
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

        OutputFileHandler output_file_handler("TestAngiogenesisSolver/TipSproutAnastamosisEvent/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        AngiogenesisSolver<3> angiogenesis_solver(p_network, output_directory);
        angiogenesis_solver.Run();

        // Test that the vessel has grown
        TS_ASSERT_DELTA(p_network->GetVessel(0)->GetEndNode()->GetLocationVector()[0], 150.0, 1.e-6);
        TS_ASSERT_DELTA(p_network->GetVessel(1)->GetStartNode()->GetLocationVector()[0], 150.0, 1.e-6);
    }

    void TestOffLatticeAngiogenesisSimulation() throw(Exception)
    {

    }
};

#endif
