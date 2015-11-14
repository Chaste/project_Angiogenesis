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

#ifndef TESTANASTA_HPP
#define TESTANASTA_HPP

#include <cxxtest/TestSuite.h>
#include "OffLatticePrwGrowthDirectionModifier.hpp"
#include "OffLatticeTipAttractionGrowthDirectionModifier.hpp"
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
#include "AngiogenesisSolver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CaVesselSegment.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FlowSolver.hpp"
#include "AbstractSproutingRule.hpp"
#include "OffLatticeRandomNormalSproutingRule.hpp"
#include "RandomNumberGenerator.hpp"

class TestAnastamosis : public AbstractCellBasedTestSuite
{

public:

    void TestParallelConterflow() throw(Exception)
    {
        RandomNumberGenerator::Instance()->Reseed(123456);

        // Make a network
        std::vector<boost::shared_ptr<VascularNode<3> > > bottom_nodes;
        for(unsigned idx=0; idx<41; idx++)
        {
            bottom_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 10.0, 0.0));
        }

        std::vector<boost::shared_ptr<VascularNode<3> > > top_nodes;
        for(unsigned idx=0; idx<41; idx++)
        {
            top_nodes.push_back(VascularNode<3>::Create(double(idx)*10, 450.0, 0.0));
        }

        boost::shared_ptr<CaVessel<3> > p_vessel1 = CaVessel<3>::Create(bottom_nodes);
        boost::shared_ptr<CaVessel<3> > p_vessel2 = CaVessel<3>::Create(top_nodes);
        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);
        p_network->SetSegmentRadii(10.0);

        for(unsigned idx=1; idx<40; idx+=2)
        {
            p_network->FormSprout(ChastePoint<3>(double(idx)*10, 10.0, 0.0), ChastePoint<3>(double(idx)*10, 20.0, 0.0));
            p_network->FormSprout(ChastePoint<3>(double(idx)*10, 450.0, 0.0), ChastePoint<3>(double(idx)*10, 440.0, 0.0));
        }

        boost::shared_ptr<OffLatticePrwGrowthDirectionModifier<3> > p_grow_direction_modifier = OffLatticePrwGrowthDirectionModifier<3>::Create();
        boost::shared_ptr<OffLatticeTipAttractionGrowthDirectionModifier<3> > p_grow_direction_modifier2 = OffLatticeTipAttractionGrowthDirectionModifier<3>::Create();
        p_grow_direction_modifier2->SetNetwork(p_network);
        boost::shared_ptr<OffLatticeRandomNormalSproutingRule<3> > p_sprouting_rule = OffLatticeRandomNormalSproutingRule<3>::Create();
        p_sprouting_rule->SetSproutingProbability(0.005);

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(40, 40);
        AngiogenesisSolver<3> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier);
        angiogenesis_solver.AddGrowthDirectionModifier(p_grow_direction_modifier2);
        angiogenesis_solver.SetSproutingRule(p_sprouting_rule);
        angiogenesis_solver.SetAnastamosisRadius(4.0);
        angiogenesis_solver.SetEndTime(40.0);

        MAKE_PTR_ARGS(OutputFileHandler, p_handler, ("TestAnastamosis/"));
        angiogenesis_solver.SetFileHandler(p_handler);
        angiogenesis_solver.Run();
    }

};

#endif
