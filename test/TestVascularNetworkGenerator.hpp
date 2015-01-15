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
#include "VascularNetworkGenerator.hpp"
#include "FakePetscSetup.hpp"

class TestVascularNetworkGenerator : public CxxTest::TestSuite
{
public:

	void TestGenerateAndWriteHexagonalNetwork() throw(Exception)
	{
		// Specify the network dimensions
		unsigned width = 50u;
		unsigned height = 25u;
		unsigned vessel_length = 5u;

		// Generate the network
		VascularNetworkGenerator<2> vascular_network_generator;
		boost::shared_ptr<CaVascularNetwork<2> > vascular_network = vascular_network_generator.GenerateHexagonalNetwork(width, height, vessel_length);

		///\todo Add some checks for the generated network, e.g. node positions, number of nodes etc.

		// Write the network to file
		OutputFileHandler output_file_handler("TestVascularNetworkGenerator", false);
		std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("HexagonalVesselNetwork.vtk");
		vascular_network->SaveVasculatureDataToFile(output_filename);
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
		VascularNetworkGenerator<3> vascular_network_generator;
		boost::shared_ptr<CaVascularNetwork<3> > vascular_network = vascular_network_generator.GenerateNetworkFromVtkFile(input_filename);

		///\todo Add some checks for the generated network, e.g. node positions, number of nodes etc.

		// Write the network to file
		OutputFileHandler output_file_handler("TestVascularNetworkGenerator", false);
		std::string output_filename = output_file_handler.GetOutputDirectoryFullPath().append("VtkVesselNetwork.vtk");
		vascular_network->SaveVasculatureDataToFile(output_filename);
     }
	#endif // CHASTE_VTK
};

#endif /*TESTVASCULARNETWORKGENERATOR_HPP_*/
