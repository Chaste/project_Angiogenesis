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

class TestMigrationWithPde : public CxxTest::TestSuite
{

public:

    void DontTestSingleVesselGrowth() throw(Exception)
    {
        // Make a network
        boost::shared_ptr<VascularNode<3> > p_node1 = VascularNode<3>::Create(0.0, 20.0, 4.0);
        boost::shared_ptr<VascularNode<3> > p_node2 = VascularNode<3>::Create(10.0, 20.0, 4.0);
        p_node2->SetIsMigrating(true);
        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(CaVesselSegment<3>::Create(p_node1, p_node2));

        boost::shared_ptr<CaVascularNetwork<3> > p_network = boost::shared_ptr<CaVascularNetwork<3> >(new CaVascularNetwork<3>());
        p_network->AddVessel(p_vessel1);
        // Set the one of the nodes to migrate
        p_network->GetVessel(0)->GetEndNode()->SetIsMigrating(true);

        // Set up the PDE domain
        boost::shared_ptr<Part<3> > p_domain = Part<3>::Create();
        p_domain->AddCuboid(100, 40, 8);

        // Choose the PDE
        boost::shared_ptr<HybridLinearEllipticPde<3> > p_pde = HybridLinearEllipticPde<3>::Create();
        p_pde->SetDiffusionConstant(0.0033);
        p_pde->SetLinearInUTerm(-2.e-7);

        // Set up and run the solver
        boost::shared_ptr<FiniteDifferenceSolver<3> > p_solver = FiniteDifferenceSolver<3>::Create();
        p_solver->SetExtents(p_domain, 4.0);
        p_solver->SetPde(p_pde);

        OutputFileHandler output_file_handler("TestMigrationWithPde/SingleVessel/", false);
        std::string output_directory = output_file_handler.GetOutputDirectoryFullPath();

        // Grow the vessel
        AbstractAngiogenesisSolver<3> angiogenesis_solver(p_network, output_directory);
        angiogenesis_solver.SetPdeSolver(p_solver);
        angiogenesis_solver.Run();
    }
};

#endif
