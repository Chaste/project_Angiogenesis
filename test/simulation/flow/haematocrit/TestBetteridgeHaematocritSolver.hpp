/*

 Copyright (c) 2005-2015, University of Oxford.
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

#ifndef TestSimpleStructuralAdaptationSolver_hpp
#define TestSimpleStructuralAdaptationSolver_hpp

#include <cxxtest/TestSuite.h>
#include "StructuralAdaptationSolver.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VasculatureGenerator.hpp"
#include "FlowSolver.hpp"
#include "VasculatureData.hpp"
#include "SimulationTime.hpp"
#include "FakePetscSetup.hpp"
#include "BetteridgeHaematocritSolver.hpp"
#include "PoiseuilleImpedanceCalculator.hpp"

class TestBetteridgeSolver : public CxxTest::TestSuite
{

public:

    void TestHexagonalNetworkBetteridgeHaematocrit() throw(Exception)
    {
        // Specify the network dimensions
        double vessel_length = 80.0;

        // Generate the network
        VasculatureGenerator<2> vascular_network_generator;
        boost::shared_ptr<VascularNetwork<2> > vascular_network = vascular_network_generator.GenerateHexagonalNetwork(800.0,
                                                                                                                        1000.0,
                                                                                                                        vessel_length);

        std::vector<ChastePoint<2> > points;
        points.push_back(ChastePoint<2>(0, 0));
        points.push_back(ChastePoint<2>(5, 0));
        std::vector<boost::shared_ptr<VascularNode<2> > > nodes;
        for(unsigned i=0; i < points.size(); i++)
        {
            nodes.push_back(boost::shared_ptr<VascularNode<2> > (VascularNode<2>::Create(points[i])));
        }

        boost::shared_ptr<VesselSegment<2> > p_segment(VesselSegment<2>::Create(nodes[0], nodes[1]));

        double radius = 10.0;
        p_segment->SetRadius(radius);
        double haematocrit = 0.45;
        p_segment->GetFlowProperties()->SetHaematocrit(haematocrit);
        vascular_network->SetSegmentProperties(p_segment);

        std::vector<std::pair<double, double> > extents = vascular_network->GetExtents();
        double y_middle = (extents[1].first + extents[1].second) /2.0;
        double x_middle = (extents[0].first + extents[0].second) /2.0;

        std::vector<boost::shared_ptr<Vessel<2> > >::iterator vessel_iterator;

        std::vector<boost::shared_ptr<Vessel<2> > > vessels = vascular_network->GetVessels();

        for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
        {
            if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetStartNode()->GetLocation()[1] >  y_middle)
                {
                    if((*vessel_iterator)->GetStartNode()->GetLocation()[0] >  x_middle)
                    {
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3320);
                    }
                }
            }
            if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetEndNode()->GetLocation()[1] >  y_middle)
                {
                    if((*vessel_iterator)->GetStartNode()->GetLocation()[0] >  x_middle)
                    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3320);
                    }
                }
            }
            if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetStartNode()->GetLocation()[1] <=  y_middle)
                {
                    if((*vessel_iterator)->GetStartNode()->GetLocation()[0] <  x_middle)
                    {
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2090);
                    }
                }
            }
            if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetEndNode()->GetLocation()[1] <=  y_middle)
                {
                    if((*vessel_iterator)->GetStartNode()->GetLocation()[0] <  x_middle)
                    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2090);
                    }
                }
            }
        }

        std::vector<boost::shared_ptr<VesselSegment<2> > > segments = vascular_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(1.e-3);
        }


        PoiseuilleImpedanceCalculator<2> impedance_calculator;
        impedance_calculator.Calculate(vascular_network);
        FlowSolver<2> solver;
        solver.SetVesselNetwork(vascular_network);
        solver.SetUp();
        solver.Solve();

        OutputFileHandler output_file_handler("TestBetteridgeHaematocritSolver", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexNet.vtp");
        vascular_network->Write(output_filename);

        boost::shared_ptr<BetteridgeHaematocritSolver<2> > p_haematocrit_calculator(new BetteridgeHaematocritSolver<2>());
        p_haematocrit_calculator->Calculate(vascular_network);

        std::string output_filename2 = output_file_handler.GetOutputDirectoryFullPath().append("HexNetHemo.vtp");
        vascular_network->Write(output_filename2);
    }
};

#endif
