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
        std::vector<unsigned> num_units(2,3);
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
        std::vector<unsigned> num_units(3, 3);
        vascular_network_generator.PatternUnitByTranslation(vascular_network, num_units);

        // Write the network to file
        OutputFileHandler output_file_handler("TestVasculatureGenerator", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexagonalVesselNetwork3d.vtp");
        vascular_network->Write(output_filename);
    }

    void TestGenerateDivergeAndConvergeNetwork() throw (Exception)
    {
        // Specify the network dimensions
        double segment_length = 20.0;
        c_vector<double, 3> start_location = zero_vector<double>(3);
        c_vector<double, 3> end_location = unit_vector<double>(3,1) * 1000.0;

        VasculatureGenerator<3> vascular_network_generator;

        boost::shared_ptr<CaVascularNetwork<3> > p_network = vascular_network_generator.GenerateDivergeAndConvergeNetwork(start_location,
                                                                                                                                       end_location,
                                                                                                                                       segment_length);
        OutputFileHandler output_file_handler("TestVasculatureGenerator", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("DivergeAndConvergeNetwork");
        p_network->Write(output_filename + ".vtp");
    }

    void TestVoronoiNetwork() throw (Exception)
    {
        // Generate the network
        VasculatureGenerator<3> vascular_network_generator;
        OutputFileHandler output_file_handler("TestVasculatureGenerator", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("VoronoiNetwork.vtp");
        boost::shared_ptr<CaVascularNetwork<3> > p_network = vascular_network_generator.GenerateVoronoiNetwork(100, 100, 100, 100);
        p_network->Write(output_filename);

    }

    void TestParallelNetworks() throw (Exception)
    {
        boost::shared_ptr<Part<3> > p_part = Part<3>::Create();
        p_part->AddCuboid(1000.0, 1000.0, 50.0);
        VasculatureGenerator<3> network_generator;
        boost::shared_ptr<CaVascularNetwork<3> > p_network = network_generator.GenerateParrallelNetwork(p_part,
                                                                                                        1.e-4,
                                                                                                        VesselDistribution::REGULAR);
        OutputFileHandler output_file_handler("TestVasculatureGenerator/Parallel", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("RegularNetwork.vtp");
        p_network->Write(output_filename);

        boost::shared_ptr<CaVascularNetwork<3> > p_network2 = network_generator.GenerateParrallelNetwork(p_part,
                                                                                                        1.e-4,
                                                                                                        VesselDistribution::UNIFORM);

        std::string output_filename2 = output_file_handler.GetOutputDirectoryFullPath().append("UniformNetwork.vtp");
        p_network2->Write(output_filename2);

        boost::shared_ptr<CaVascularNetwork<3> > p_network3 = network_generator.GenerateParrallelNetwork(p_part,
                                                                                                        1.e-4,
                                                                                                        VesselDistribution::UNIFORM,
                                                                                                        20.0);

        std::string output_filename3 = output_file_handler.GetOutputDirectoryFullPath().append("UniformExclusionNetwork.vtp");
        p_network3->Write(output_filename3);

        boost::shared_ptr<CaVascularNetwork<3> > p_network4 = network_generator.GenerateParrallelNetwork(p_part,
                                                                                                        1.e-4,
                                                                                                        VesselDistribution::TWO_LAYER);

        std::string output_filename4 = output_file_handler.GetOutputDirectoryFullPath().append("TwoLayerNetwork.vtp");
        p_network3->Write(output_filename4);

        boost::shared_ptr<CaVascularNetwork<3> > p_network5 = network_generator.GenerateParrallelNetwork(p_part,
                                                                                                        1.e-4,
                                                                                                        VesselDistribution::TWO_LAYER,
                                                                                                        20.0);

        std::string output_filename5 = output_file_handler.GetOutputDirectoryFullPath().append("TwoLayerExclusionNetwork.vtp");
        p_network3->Write(output_filename5);
    }

    void Test3dNetworks() throw (Exception)
    {
        boost::shared_ptr<Part<3> > p_part = Part<3>::Create();
        p_part->AddCuboid(1000.0, 1000.0, 1000.0);
        std::vector<double> density;
        density.push_back(8.e-5);
        density.push_back(8.e-5);
        density.push_back(1.e-5);
        VasculatureGenerator<3> network_generator;
        boost::shared_ptr<CaVascularNetwork<3> > p_network = network_generator.Generate3dNetwork(p_part,
                                                                                                        density,
                                                                                                        VesselDistribution::REGULAR);
        OutputFileHandler output_file_handler("TestVasculatureGenerator/3d", false);
        std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("RegularNetwork.vtp");
        p_network->Write(output_filename);

        std::vector<double> density2;
        density2.push_back(4.e-5);
        density2.push_back(4.e-5);
        density2.push_back(4.e-5);
        boost::shared_ptr<CaVascularNetwork<3> > p_network2 = network_generator.Generate3dNetwork(p_part,
                                                                                                  density2,
                                                                                                        VesselDistribution::UNIFORM);

        std::string output_filename2 = output_file_handler.GetOutputDirectoryFullPath().append("UniformNetwork.vtp");
        p_network2->Write(output_filename2);

        boost::shared_ptr<CaVascularNetwork<3> > p_network3 = network_generator.Generate3dNetwork(p_part,
                                                                                                         density2,
                                                                                                        VesselDistribution::UNIFORM,
                                                                                                        20.0);

        std::string output_filename3 = output_file_handler.GetOutputDirectoryFullPath().append("UniformExclusionNetwork.vtp");
        p_network3->Write(output_filename3);

        boost::shared_ptr<CaVascularNetwork<3> > p_network4 = network_generator.Generate3dNetwork(p_part,
                                                                                                         density2,
                                                                                                        VesselDistribution::TWO_LAYER);

        std::string output_filename4 = output_file_handler.GetOutputDirectoryFullPath().append("TwoLayerNetwork.vtp");
        p_network3->Write(output_filename4);

        boost::shared_ptr<CaVascularNetwork<3> > p_network5 = network_generator.Generate3dNetwork(p_part,
                                                                                                         density2,
                                                                                                        VesselDistribution::TWO_LAYER,
                                                                                                        20.0);

        std::string output_filename5 = output_file_handler.GetOutputDirectoryFullPath().append("TwoLayerExclusionNetwork.vtp");
        p_network3->Write(output_filename5);
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
