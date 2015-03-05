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

#ifndef TESTSIMPLEFLOWSOLVER_HPP_
#define TESTSIMPLEFLOWSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VasculatureGenerator.hpp"
#include "FakePetscSetup.hpp"
#include "SimpleFlowSolver.hpp"
#include "VasculatureData.hpp"

class TestSimpleFlowSolver : public CxxTest::TestSuite
{
public:

    void TestFlowThroughHexagonalNetwork() throw(Exception)
    {
        // Specify the network dimensions
        double vessel_length = 80.0;

        // Generate the network
        VasculatureGenerator<2> vascular_network_generator;
        boost::shared_ptr<CaVascularNetwork<2> > vascular_network = vascular_network_generator.GenerateHexagonalNetwork(1000,1000,vessel_length);

        VasculatureData data;
        double impedance = 10.0;
        double flow_rate = 0.0;
        data.SetData("Impedance", impedance);
        data.SetData("Flow Rate", flow_rate);
        vascular_network->SetVesselData(data);
        vascular_network->SetSegmentData(data);

        VasculatureData node_data;
        double pressure = 10.0;
        node_data.SetData("Pressure", pressure);
        bool is_input = false;
        node_data.SetData("Is Input", is_input);
        bool is_output = false;
        node_data.SetData("Is Output", is_output);
        vascular_network->SetNodeData(node_data);

        std::vector<std::pair<double, double> > extents = vascular_network->GetExtents();
        double y_middle = (extents[1].first + extents[1].second) /2.0;

        typename std::vector<boost::shared_ptr<CaVessel<2> > >::iterator vessel_iterator;

        std::vector<boost::shared_ptr<CaVessel<2> > > vessels = vascular_network->GetVessels();

        for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
        {
            if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetStartNode()->GetLocation()[1] >  y_middle)
                {
                    (*vessel_iterator)->GetStartNode()->GetDataContainer().SetData("Is Input", true);
                    (*vessel_iterator)->GetStartNode()->GetDataContainer().SetData("Pressure", 3393);
                }
            }
            if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetEndNode()->GetLocation()[1] >  y_middle)
                {
                    (*vessel_iterator)->GetEndNode()->GetDataContainer().SetData("Is Input", true);
                    (*vessel_iterator)->GetEndNode()->GetDataContainer().SetData("Pressure", 3393);
                }
            }
            if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetStartNode()->GetLocation()[1] <=  y_middle)
                {
                    (*vessel_iterator)->GetStartNode()->GetDataContainer().SetData("Is Output", true);
                    (*vessel_iterator)->GetStartNode()->GetDataContainer().SetData("Pressure", 1993);
                }
            }
            if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
            {
                if((*vessel_iterator)->GetEndNode()->GetLocation()[1] <=  y_middle)
                {
                    (*vessel_iterator)->GetEndNode()->GetDataContainer().SetData("Is Output", true);
                    (*vessel_iterator)->GetEndNode()->GetDataContainer().SetData("Pressure", 1993);
                }
            }

        }

        SimpleFlowSolver<2> solver;
        solver.Implement(vascular_network);

        // Write the network to file
        OutputFileHandler output_file_handler("TestSimpleFlowSolver", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexagonalVesselNetwork.vtp");
        vascular_network->WriteToFile(output_filename);
    }

};

#endif /*TESTSIMPLEFLOWSOLVER_HPP_*/
