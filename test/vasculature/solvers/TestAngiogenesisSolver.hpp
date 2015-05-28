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
#include "CaVascularNetwork.hpp"
#include "AngiogenesisSolver.hpp"

class TestAngiogenesisSolver : public CxxTest::TestSuite
{

public:

    void TestSingleVesselGrowth() throw(Exception)
    {
        // Make a network
        VasculatureGenerator<3> network_generator;
        boost::shared_ptr<CaVascularNetwork<3> > p_network = network_generator.GenerateSingleVessel();

        // Set the one of the nodes to migrate
        p_network->GetVessel(0)->GetEndNode()->SetIsMigrating(true);

        // Grow the vessel
        AngiogenesisSolver<3> angiogenesis_solver(p_network);
        angiogenesis_solver.Run();

        // Test that the vessel has grown
        TS_ASSERT(p_network->GetNumberOfNodes() == 3u);
        TS_ASSERT_DELTA(p_network->GetVessel(1)->GetEndNode()->GetLocationVector[0], 200.0, 1.e-6);
    }
};

#endif
