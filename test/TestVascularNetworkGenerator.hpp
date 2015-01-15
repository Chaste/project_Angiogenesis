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

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <algorithm>

#include "FileFinder.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "VascularNetworkGenerator.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestVascularNetworkGenerator : public CxxTest::TestSuite
{
public:

	void TestVascularNetworkGeneratorWithHexagonallyTesselatedVesselNetwork() throw(Exception)
	{

		try
		{

			// Create a simple 2D PottsMesh
			unsigned width = 50;
			unsigned height = 25;
			unsigned vessel_length = 5;

			boost::shared_ptr<CaVessel<2> > prototypeVessel(new CaVessel<2>());
			prototypeVessel->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(0.0,"unitless"),3.0*pow(10.0,(-6)));
			VascularNetworkGenerator<2> vascularNetworkGenerator(prototypeVessel);

			boost::shared_ptr<CaVascularNetwork<2> > vascularNetwork = vascularNetworkGenerator.GenerateHexagonallyTesselatedVascularNetwork(width,height,vessel_length,"North East");

			// write vascular cell data out to file

			std::string output_directory = "TestVascularNetworkGenerator";
			std::string filename;
			std::stringstream sstm;
			OutputFileHandler output_file_handler(output_directory, false);

			sstm << output_file_handler.GetOutputDirectoryFullPath() << "HexagonallyTessellatedVesselNetwork.vtk";

			filename = sstm.str();

			vascularNetwork->SaveVasculatureDataToFile(filename);

		}
		catch(Exception &e)
		{
			std::cout << e.GetMessage() << std::endl;
		}

	}

	void TestVascularNetworkGeneratorWithTapmeierNetwork() throw(Exception)
        		{

		try
		{

			boost::shared_ptr<CaVessel<3> > prototypeVessel(new CaVessel<3>());
			prototypeVessel->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(0.0,"unitless"),3.0*pow(10.0,(-6)));
			VascularNetworkGenerator<3> vascularNetworkGenerator(prototypeVessel);

			FileFinder fileFinder("projects/Angiogenesis/test/data/tapmeier_network.vtp",RelativeTo::ChasteSourceRoot);

			std::string input_filename = fileFinder.GetAbsolutePath();
			TS_ASSERT(fileFinder.IsFile());
			TS_ASSERT(fileFinder.Exists());

			boost::shared_ptr<CaVascularNetwork<3> > vascularNetwork = vascularNetworkGenerator.GenerateVascularNetworkFromVtkFile(input_filename);

			// write vascular cell data out to file

			std::string output_directory = "TestVascularNetworkGenerator";
			std::string output_filename;
			std::stringstream sstm;
			OutputFileHandler output_file_handler(output_directory, false);

			sstm << output_file_handler.GetOutputDirectoryFullPath() << "TapmeierNetwork.vtk";

			output_filename = sstm.str();

			vascularNetwork->SaveVasculatureDataToFile(output_filename);

		}
		catch(Exception &e)
		{
			std::cout << e.GetMessage() << std::endl;
		}

        		}

};

#endif /*TESTVASCULARNETWORKGENERATOR_HPP_*/
