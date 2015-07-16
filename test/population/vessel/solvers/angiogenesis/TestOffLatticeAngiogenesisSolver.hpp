//
//  TestSimpleStructuralAdaptationSolver.hpp
//  VascularTumourGrowthModellingFramework
//
//  Created by Anthony Connor on 03/12/2012.
//  Copyright (c) 2012 Anthony Connor. All rights reserved.
//

#ifndef TestOffLatticeAngiogenesisSolver_hpp
#define TestOffLatticeAngiogenesisSolver_hpp

#include <cxxtest/TestSuite.h>
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VasculatureGenerator.hpp"
#include "VascularNode.hpp"
#include "CaVesselSegment.hpp"
#include "CaVessel.hpp"
#include "CaVascularNetwork.hpp"
#include "OffLatticeAngiogenesisSolver.hpp"

class TestOffLatticeAngiogenesisSolver : public CxxTest::TestSuite
{

public:

    void TestSingleVesselGrowth() throw(Exception)
    {
        // Make a network
        boost::shared_ptr<VascularNode<3> > p_node1 = VascularNode<3>::Create(0.0, 0.0, 0.0);
        boost::shared_ptr<VascularNode<3> > p_node2 = VascularNode<3>::Create(10.0, 0.0, 0.0);
        p_node2->SetIsMigrating(true);
        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(CaVesselSegment<3>::Create(p_node1, p_node2));

        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);

        OutputFileHandler output_file_handler("TestOffLatticeAngiogenesisSolver/SingleVesselGrowth/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        OffLatticeAngiogenesisSolver<3> angiogenesis_solver(p_network, output_directory);
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

        OutputFileHandler output_file_handler("TestOffLatticeAngiogenesisSolver/Multivessel/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        OffLatticeAngiogenesisSolver<3> angiogenesis_solver(p_network, output_directory);
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

        OutputFileHandler output_file_handler("TestOffLatticeAngiogenesisSolver/MultiSprout/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        OffLatticeAngiogenesisSolver<3> angiogenesis_solver(p_network, output_directory);
        angiogenesis_solver.Run();
    }
};

#endif
