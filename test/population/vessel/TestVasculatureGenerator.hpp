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

#ifndef TESTVASCULARNETWORKGENERATOR_HPP_
#define TESTVASCULARNETWORKGENERATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VasculatureGenerator.hpp"
#include "FakePetscSetup.hpp"

class TestVasculatureGenerator : public CxxTest::TestSuite
{
public:

    void TestGenerateAndWriteHexagonalNetwork() throw (Exception)
    {
        // Specify the network dimensions
        double vessel_length = 5.0;

        // Generate the network
        VasculatureGenerator<2> vascular_network_generator;
        boost::shared_ptr<CaVascularNetwork<2> > vascular_network = vascular_network_generator.GenerateHexagonalUnit(vessel_length);

        // Pattern the unit
        std::vector<unsigned> num_units;
        num_units.push_back(3);
        num_units.push_back(3);
        vascular_network_generator.PatternUnitByTranslation(vascular_network, num_units);

        // Write the network to file
        OutputFileHandler output_file_handler("TestVasculatureGenerator", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexagonalVesselNetwork.vtp");
        vascular_network->Write(output_filename);
    }

    void TestGenerate3dHexagonalNetwork() throw (Exception)
    {
        // Specify the network dimensions
        double vessel_length = 40.0;

        // Generate the network
        VasculatureGenerator<3> vascular_network_generator;
        boost::shared_ptr<CaVascularNetwork<3> > vascular_network = vascular_network_generator.GenerateHexagonalUnit(vessel_length);

        // Pattern the unit
        std::vector<unsigned> num_units;
        num_units.push_back(3);
        num_units.push_back(3);
        num_units.push_back(3);
        vascular_network_generator.PatternUnitByTranslation(vascular_network, num_units);

        // Write the network to file
        OutputFileHandler output_file_handler("TestVasculatureGenerator", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexagonalVesselNetwork3d.vtp");
        vascular_network->Write(output_filename);
    }

    void TestGenerateSimpleDivergeAndConvergeNetwork() throw (Exception)
    {
        // Specify the network dimensions
        double segment_length = 20.0;
        c_vector<double, 3> start_location;
        c_vector<double, 3> end_location;

        start_location[0] = 0.0;
        start_location[1] = 0.0;
        start_location[2] = 0.0;

        end_location[0] = 1000.0;
        end_location[1] = 0.0;
        end_location[2] = 0.0;

        // Generate the network
        VasculatureGenerator<3> vascular_network_generator;
        OutputFileHandler output_file_handler("TestVasculatureGenerator", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("DivergeAndConvergeNetwork");
        boost::shared_ptr<CaVascularNetwork<3> > vascular_network = vascular_network_generator.GenerateSimpleDivergeAndConvergeNetwork(start_location,
                                                                                                                                       end_location,
                                                                                                                                       segment_length,
                                                                                                                                       output_filename);
    }

#ifdef CHASTE_VTK
    void TestGeneratorWithVtkInput() throw(Exception)
    {
        // Locate the input file
        FileFinder fileFinder("projects/Angiogenesis/test/data/tapmeier_network.vtp", RelativeTo::ChasteSourceRoot);
        TS_ASSERT(fileFinder.Exists());
        TS_ASSERT(fileFinder.IsFile());
        std::string input_filename = fileFinder.GetAbsolutePath();

        // Generate the network
        VasculatureGenerator<3> vascular_network_generator;
        boost::shared_ptr<CaVascularNetwork<3> > vascular_network = vascular_network_generator.GenerateNetworkFromVtkFile(input_filename);

        // Write the network to file
        OutputFileHandler output_file_handler("TestVasculatureGenerator", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("VtkVesselNetwork.vtp");
        vascular_network->MergeCoincidentNodes();
        vascular_network->Write(output_filename);

        std::string output_filename1 = output_file_handler.GetOutputDirectoryFullPath().append("VesselNetworkConnectivityGraph.gv");
        vascular_network->WriteConnectivity(output_filename1);

    }
#endif // CHASTE_VTK
};

#endif /*TESTVASCULARNETWORKGENERATOR_HPP_*/