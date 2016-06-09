/*

Copyright (c) 2005-2014, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTREGRESSIONSSOLVER_HPP
#define TESTREGRESSIONSSOLVER_HPP

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VascularNetwork.hpp"
#include "WallShearStressBasedRegressionSolver.hpp"
#include "AngiogenesisSolver.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestRegressionSolver : public AbstractCellBasedWithTimingsTestSuite
{

public:

    void TestSingleVesselRegression() throw(Exception)
    {
        // Make a vessel
        boost::shared_ptr<VascularNode<2> > p_node1 = VascularNode<2>::Create(0.0, 0.0);
        boost::shared_ptr<VascularNode<2> > p_node2 = VascularNode<2>::Create(0.0,100.0);
        boost::shared_ptr<Vessel<2> > p_vessel = Vessel<2>::Create(p_node1, p_node2);
        boost::shared_ptr<VascularNetwork<2> > p_network = VascularNetwork<2>::Create();
        p_network->AddVessel(p_vessel);

        // Set a wall shear stress below the threshold
        double wss_threshold = 10.0;
        double vessel_wss = 5.0;
        p_vessel->GetSegment(0)->GetFlowProperties()->SetWallShearStress(vessel_wss);

        // Set up the regression solver
        WallShearStressBasedRegressionSolver<2> regression_solver = WallShearStressBasedRegressionSolver<2>();
        regression_solver.SetVesselNetwork(p_network);
        regression_solver.SetLowWallShearStressThreshold(wss_threshold);
        regression_solver.SetMaximumTimeWithLowWallShearStress(5);

        // Run the solver for six 'increments'
        for(unsigned idx=0 ; idx<6; idx++)
        {
            regression_solver.Increment();
            TS_ASSERT(p_network->GetVessel(0)->HasRegressionTimerStarted());
        }

        TS_ASSERT_EQUALS(p_network->GetNumberOfVessels(), 0);
    }
};

#endif // TESTANGIOGENESISSOLVER_HPP
